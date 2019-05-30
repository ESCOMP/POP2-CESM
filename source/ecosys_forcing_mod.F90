module ecosys_forcing_mod

  ! !DESCRIPTION:
  !  This module sets up data types to keep track of the forcing data
  !  being sent to MARBL. Data may be passed through to MARBL from the flux
  !  coupler, read from a file, or set to a constant value based on namelist
  !  variables.
  !
  !  This module mostly handles surface flux forcing, but also handles some
  !  interior tendency forcing -- namely tracer restoring.

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
  use ecosys_tracers_and_saved_state_mod, only : dic_ind, alk_ind, dic_alt_co2_ind, alk_alt_co2_ind
  use ecosys_tracers_and_saved_state_mod, only : di13c_ind, di14c_ind
  use ecosys_tracers_and_saved_state_mod, only : no3_ind, po4_ind, don_ind, donr_ind, dop_ind, dopr_ind
  use ecosys_tracers_and_saved_state_mod, only : sio3_ind, fe_ind, doc_ind, docr_ind, do13ctot_ind, do14ctot_ind

  use forcing_timeseries_mod, only : forcing_timeseries_dataset

  implicit none
  private

  public :: ecosys_forcing_init
  public :: ecosys_forcing_set_surface_time_varying_forcing_data
  public :: ecosys_forcing_comp_stf_riv
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
     integer   (kind=int_kind) :: year_first
     integer   (kind=int_kind) :: year_last
     integer   (kind=int_kind) :: year_align
     character (char_len)      :: tintalgo
     character (char_len)      :: taxMode
     integer   (kind=int_kind) :: strdata_inputlist_ind
     integer   (kind=int_kind) :: strdata_var_ind
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
     real(kind=r8)                        :: unit_conv_factor ! unit conversion factor, incorporates scale_factor
     logical(kind=log_kind)               :: ltime_varying    ! does this field vary in time
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
     integer (int_kind)                 :: rank
     logical(kind=log_kind)             :: ldim3_is_depth     ! is 3rd dimension depth
     real(r8), allocatable              :: field_0d(:,:,:)    ! (nx, ny, iblock)
     real(r8), allocatable              :: field_1d(:,:,:,:)  ! (nx, ny, [any dim], iblock)
  contains
     procedure, public :: add_forcing_field => forcing_fields_add
  end type forcing_fields_type

  type(forcing_fields_type), dimension(:), allocatable, public :: riv_flux_forcing_fields
  type(forcing_fields_type), dimension(:), allocatable, public :: surface_flux_forcings
  type(forcing_fields_type), dimension(:), allocatable, public :: interior_tendency_forcings

  !---------------------------------------------------------------------
  !  Variables read in via &ecosys_forcing_data_nml
  !---------------------------------------------------------------------

  character(char_len) :: dust_flux_source             ! option for atmospheric dust deposition
  type(tracer_read)   :: dust_flux_input              ! namelist input for dust_flux
  character(char_len) :: iron_flux_source             ! option for atmospheric iron deposition
  real(r8)            :: dust_ratio_thres             ! coarse/fine dust ratio threshold, used in iron_flux_source=='driver-derived' computation
  real(r8)            :: fe_bioavail_frac_offset
  real(r8)            :: dust_ratio_to_fe_bioavail_frac_r
  type(tracer_read)   :: iron_flux_input              ! namelist input for iron_flux
  type(tracer_read)   :: fesedflux_input              ! namelist input for fesedflux
  type(tracer_read)   :: feventflux_input             ! namelist input for feventflux
  character(char_len) :: o2_consumption_scalef_opt    ! option for specification of o2_consumption_scalef
  real(r8)            :: o2_consumption_scalef_const  ! constant for o2_consumption_scalef_opt=const
  type(tracer_read)   :: o2_consumption_scalef_input  ! file info for o2_consumption_scalef_opt=file_time_invariant
  character(char_len) :: p_remin_scalef_opt           ! option for specification of p_remin_scalef
  real(r8)            :: p_remin_scalef_const         ! constant for p_remin_scalef_opt=const
  type(tracer_read)   :: p_remin_scalef_input         ! file info for p_remin_scalef_opt=file_time_invariant
  logical(log_kind)   :: lignore_driver_ndep
  character(char_len) :: ndep_data_type               ! type of ndep forcing
  type(tracer_read)   :: nox_flux_monthly_input       ! namelist input for nox_flux_monthly
  type(tracer_read)   :: nhy_flux_monthly_input       ! namelist input for nhy_flux_monthly
  integer(int_kind)   :: ndep_shr_stream_year_first   ! first year in stream to use
  integer(int_kind)   :: ndep_shr_stream_year_last    ! last year in stream to use
  integer(int_kind)   :: ndep_shr_stream_year_align   ! align ndep_shr_stream_year_first with this model year
  character(char_len) :: ndep_shr_stream_file         ! file containing domain and input data
  real(r8)            :: ndep_shr_stream_scale_factor ! unit conversion factor
  character(char_len) :: gas_flux_forcing_opt         ! option for forcing gas fluxes
  character(char_len) :: gas_flux_forcing_file        ! file containing gas flux forcing fields
  type(tracer_read)   :: gas_flux_fice                ! ice fraction for gas fluxes
  type(tracer_read)   :: gas_flux_ws                  ! wind speed for gas fluxes
  type(tracer_read)   :: gas_flux_ap                  ! atmospheric pressure for gas fluxes
  character(char_len) :: atm_co2_opt                  ! option for atmospheric co2 concentration
  real(r8)            :: atm_co2_const                ! value of atmospheric co2 (ppm, dry-air, 1 atm)
  character(char_len) :: atm_alt_co2_opt              ! option for atmospheric alternative CO2
  real(r8)            :: atm_alt_co2_const            ! value of atmospheric alternative co2 (ppm, dry-air, 1 atm)
  logical(log_kind)   :: liron_patch                  ! flag for iron patch fertilization
  character(char_len) :: iron_patch_flux_filename     ! file containing name of iron patch file
  integer(int_kind)   :: iron_patch_month             ! integer month to add patch flux
  integer(int_kind)   :: ciso_atm_model_year          ! arbitrary model year
  integer(int_kind)   :: ciso_atm_data_year           ! year in atmospheric ciso data that corresponds to ciso_atm_model_year
  real(r8)            :: ciso_atm_d13c_const          ! atmospheric d13C constant [permil]
  real(r8)            :: ciso_atm_d14c_const          ! atmospheric D14C constant [permil]
  real(r8),dimension(3)::ciso_atm_d14c_lat_band_vals  ! atmospheric D14C constant [permil]
  character(char_len) :: ciso_atm_d13c_opt            ! option for d13C (varying or constant forcing)
  character(char_len) :: ciso_atm_d13c_filename       ! filename for varying atm d13C
  character(char_len) :: ciso_atm_d14c_opt            ! option for d14C (varying or constant forcing)
  character(char_len) :: ciso_atm_d14c_filename       ! filename for varying atm D14C

  type (forcing_timeseries_dataset) :: &
    ciso_atm_d13c_forcing_dataset                     ! data structure for atm d13C timeseries

  !-----------------------------------------------------------------------
  !  tracer restoring related variables
  !-----------------------------------------------------------------------

  character(char_len), dimension(marbl_tracer_cnt) :: restorable_tracer_names    ! list of tracers that we can apply restoring to
                                                                                 !   (i.e., for which restoring fields are available)
  character(char_len), dimension(marbl_tracer_cnt) :: restore_data_filenames     ! list of files containing the restoring fields
  character(char_len), dimension(marbl_tracer_cnt) :: restore_data_file_varnames ! file varnames in restore_data_filenames corresponding to tracers names in restorable_tracer_names
  integer(int_kind), dimension(marbl_tracer_cnt)   :: restore_year_first         ! first year in stream to use
  integer(int_kind), dimension(marbl_tracer_cnt)   :: restore_year_last          ! last year in stream to use
  integer(int_kind), dimension(marbl_tracer_cnt)   :: restore_year_align         ! align restore_year_first with this model year
  real(r8), dimension(marbl_tracer_cnt)            :: restore_scale_factor       ! unit conversion factor

  ! Also need to know what time scale to restore on
  character(char_len)                              :: restore_inv_tau_opt
  real(r8)                                         :: restore_inv_tau_const
  type(tracer_read), target                        :: restore_inv_tau_input

  ! Is NDEP available from the driver?
  logical(log_kind), public :: ldriver_has_ndep
  logical(log_kind), public :: ldriver_has_atm_co2_diag
  logical(log_kind), public :: ldriver_has_atm_co2_prog

  ! Data type for reading interior tendency forcing from shr_stream
  type (strdata_input_type), pointer :: interior_strdata_inputlist_ptr(:)

  !-----------------------------------------------------------------------
  !  riv_flux forcing
  !-----------------------------------------------------------------------

  type (strdata_input_type), pointer :: riv_flux_strdata_inputlist_ptr(:)

  !-----------------------------------------------------------------------
  !  surface flux forcing
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

  type (strdata_input_type), pointer :: surface_strdata_inputlist_ptr(:)

  ! Surface flux and interior tendency forcing fields
  real(r8), allocatable, target :: iron_patch_flux(:,:,:)
  real(r8)                      :: dust_flux_in(nx_block, ny_block, max_blocks_clinic)

  integer (int_kind) :: ecosys_pre_sflux_timer
  integer (int_kind) :: ecosys_riv_flux_strdata_create_timer
  integer (int_kind) :: ecosys_surface_strdata_create_timer
  integer (int_kind) :: ecosys_interior_strdata_create_timer
  integer (int_kind) :: ecosys_riv_flux_strdata_advance_timer
  integer (int_kind) :: ecosys_surface_strdata_advance_timer
  integer (int_kind) :: ecosys_interior_strdata_advance_timer

  ! Some surface flux forcing fields need special treatment, so we store indices
  integer(int_kind) :: dust_dep_ind = 0, &
                       Fe_dep_ind   = 0, &
                       bc_dep_ind   = 0, &
                       box_atm_co2_ind     = 0, &
                       box_atm_co2_dup_ind = 0, &
                       ifrac_ind    = 0, &
                       ap_ind       = 0, &
                       sst_ind      = 0, &
                       sss_ind      = 0, &
                       ext_C_flux_ind = 0, &
                       ext_P_flux_ind = 0, &
                       ext_Si_flux_ind = 0, &
                       u10sqr_ind   = 0, &
                       d13c_ind     = 0, &
                       d14c_ind     = 0

  integer(int_kind) :: din_riv_flux_ind = 0, &
                       dip_riv_flux_ind = 0, &
                       don_riv_flux_ind = 0, &
                       dop_riv_flux_ind = 0, &
                       dsi_riv_flux_ind = 0, &
                       dfe_riv_flux_ind = 0, &
                       dic_riv_flux_ind = 0, &
                       alk_riv_flux_ind = 0, &
                       doc_riv_flux_ind = 0

  ! We also need to track the indices of all the interior tendency forcing fields
  integer(int_kind), public :: dustflux_ind       = 0, &
                       PAR_col_frac_ind   = 0, &
                       surf_shortwave_ind = 0, &
                       potemp_ind         = 0, &
                       salinity_ind       = 0, &
                       pressure_ind       = 0, &
                       fesedflux_ind      = 0

  !-----------------------------------------------------------------------
  ! Other private variables
  !-----------------------------------------------------------------------

  ! Virtual fluxes
  real(r8), dimension(marbl_tracer_cnt) :: surf_avg                      ! average surface tracer values

  real(r8) :: iron_frac_in_atm_fine_dust
  real(r8) :: iron_frac_in_atm_coarse_dust
  real(r8) :: iron_frac_in_seaice_dust
  real(r8) :: iron_frac_in_atm_bc
  real(r8) :: iron_frac_in_seaice_bc

  real(r8) :: d14c_glo_avg       ! global average D14C over the ocean, computed from current D14C field

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine ecosys_forcing_init(ciso_on, land_mask,                     &
                                 marbl_req_surface_flux_forcings,        &
                                 marbl_req_interior_tendency_forcings,   &
                                 forcing_nml,                            &
                                 lhas_riv_flux)

    use ecosys_tracers_and_saved_state_mod, only : set_defaults_tracer_read

    use marbl_interface_public_types, only : marbl_forcing_fields_type

    use constants, only : delim_fmt, ndelim_fmt

    use domain, only : distrb_clinic

    use mcog, only : mcog_nbins

    use ecosys_forcing_saved_state_mod, only : lbox_atm_co2, box_atm_co2_init_val

    logical,                         intent(in)    :: ciso_on
    logical,                         intent(in)    :: land_mask(:,:,:)
    type(marbl_forcing_fields_type), intent(in)    :: marbl_req_surface_flux_forcings(:)
    type(marbl_forcing_fields_type), intent(in)    :: marbl_req_interior_tendency_forcings(:)
    character(len=*),                intent(in)    :: forcing_nml
    logical,                         intent(out)   :: lhas_riv_flux(:) ! true if a tracer has a riverine flux

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter  :: subname = 'ecosys_forcing_mod:ecosys_forcing_init'
    character(len=char_len)      :: err_msg

    character(char_len)      :: marbl_varname, tracer_name, units
    integer (int_kind)       :: nml_error                  ! error flag for nml read
    character(char_len_long) :: ioerror_msg
    integer                  :: m, n
    type(forcing_monthly_every_ts), pointer :: file_details
    logical                  :: var_processed

    character(char_len) :: riv_flux_shr_stream_file
    integer(int_kind)   :: riv_flux_shr_stream_year_first
    integer(int_kind)   :: riv_flux_shr_stream_year_last
    integer(int_kind)   :: riv_flux_shr_stream_year_align
    character(char_len) :: riv_flux_din_file_varname
    real(r8)            :: riv_flux_din_scale_factor
    character(char_len) :: riv_flux_dip_file_varname
    real(r8)            :: riv_flux_dip_scale_factor
    character(char_len) :: riv_flux_don_file_varname
    real(r8)            :: riv_flux_don_scale_factor
    character(char_len) :: riv_flux_dop_file_varname
    real(r8)            :: riv_flux_dop_scale_factor
    character(char_len) :: riv_flux_dsi_file_varname
    real(r8)            :: riv_flux_dsi_scale_factor
    character(char_len) :: riv_flux_dfe_file_varname
    real(r8)            :: riv_flux_dfe_scale_factor
    character(char_len) :: riv_flux_dic_file_varname
    real(r8)            :: riv_flux_dic_scale_factor
    character(char_len) :: riv_flux_alk_file_varname
    real(r8)            :: riv_flux_alk_scale_factor
    character(char_len) :: riv_flux_doc_file_varname
    real(r8)            :: riv_flux_doc_scale_factor

    ! Virtual fluxes set in namelist
    real(r8) :: surf_avg_dic_const
    real(r8) :: surf_avg_alk_const
    real(r8) :: surf_avg_di13c_const
    real(r8) :: surf_avg_di14c_const

    namelist /ecosys_forcing_data_nml/                                        &
         dust_flux_source, dust_flux_input, iron_flux_source,                 &
         dust_ratio_thres, fe_bioavail_frac_offset, dust_ratio_to_fe_bioavail_frac_r, &
         iron_flux_input, fesedflux_input, feventflux_input,                  &
         o2_consumption_scalef_opt, o2_consumption_scalef_const,              &
         o2_consumption_scalef_input,                                         &
         p_remin_scalef_opt, p_remin_scalef_const, p_remin_scalef_input,      &
         lignore_driver_ndep, ndep_data_type,                                 &
         nox_flux_monthly_input, nhy_flux_monthly_input,                      &
         ndep_shr_stream_year_first, ndep_shr_stream_year_last,               &
         ndep_shr_stream_year_align, ndep_shr_stream_file,                    &
         ndep_shr_stream_scale_factor,                                        &
         riv_flux_shr_stream_file, riv_flux_shr_stream_year_first,            &
         riv_flux_shr_stream_year_last, riv_flux_shr_stream_year_align,       &
         riv_flux_din_file_varname, riv_flux_din_scale_factor,                &
         riv_flux_dip_file_varname, riv_flux_dip_scale_factor,                &
         riv_flux_don_file_varname, riv_flux_don_scale_factor,                &
         riv_flux_dop_file_varname, riv_flux_dop_scale_factor,                &
         riv_flux_dsi_file_varname, riv_flux_dsi_scale_factor,                &
         riv_flux_dfe_file_varname, riv_flux_dfe_scale_factor,                &
         riv_flux_dic_file_varname, riv_flux_dic_scale_factor,                &
         riv_flux_alk_file_varname, riv_flux_alk_scale_factor,                &
         riv_flux_doc_file_varname, riv_flux_doc_scale_factor,                &
         gas_flux_forcing_opt,                                                &
         gas_flux_forcing_file, gas_flux_fice, gas_flux_ws, gas_flux_ap,      &
         atm_co2_opt, atm_co2_const, box_atm_co2_init_val,                    &
         atm_alt_co2_opt, atm_alt_co2_const,                                  &
         liron_patch, iron_patch_flux_filename, iron_patch_month,             &
         ciso_atm_d13c_opt, ciso_atm_d13c_const, ciso_atm_d13c_filename,      &
         ciso_atm_d14c_opt, ciso_atm_d14c_const, ciso_atm_d14c_lat_band_vals, &
         ciso_atm_d14c_filename,                                              &
         ciso_atm_model_year, ciso_atm_data_year, restorable_tracer_names,    &
         restore_data_filenames, restore_data_file_varnames,                  &
         restore_year_first, restore_year_last, restore_year_align,           &
         restore_scale_factor,                                                &
         restore_inv_tau_opt, restore_inv_tau_const, restore_inv_tau_input,   &
         surf_avg_alk_const, surf_avg_dic_const,                              &
         surf_avg_di13c_const, surf_avg_di14c_const,                          &
         iron_frac_in_atm_fine_dust, iron_frac_in_atm_coarse_dust,            &
         iron_frac_in_seaice_dust, iron_frac_in_atm_bc, iron_frac_in_seaice_bc

    !-----------------------------------------------------------------------
    !  &ecosys_forcing_data_nml
    !-----------------------------------------------------------------------

    gas_flux_forcing_opt  = 'drv'
    gas_flux_forcing_file = 'unknown'
    call set_defaults_tracer_read(gas_flux_fice, file_varname='FICE')
    call set_defaults_tracer_read(gas_flux_ws, file_varname='XKW')
    call set_defaults_tracer_read(gas_flux_ap, file_varname='P')
    dust_flux_source             = 'driver'
    call set_defaults_tracer_read(dust_flux_input, file_varname='dust_flux')
    iron_flux_source             = 'driver-derived'
    dust_ratio_thres             = 55.0_r8
    fe_bioavail_frac_offset      = 0.01_r8
    dust_ratio_to_fe_bioavail_frac_r = 170.0_r8
    call set_defaults_tracer_read(iron_flux_input, file_varname='iron_flux')
    call set_defaults_tracer_read(fesedflux_input, file_varname='FESEDFLUXIN')
    call set_defaults_tracer_read(feventflux_input, file_varname='FESEDFLUXIN')
    o2_consumption_scalef_opt   = 'const'
    o2_consumption_scalef_const = c1
    call set_defaults_tracer_read(o2_consumption_scalef_input, file_varname='o2_consumption_scalef')
    p_remin_scalef_opt   = 'const'
    p_remin_scalef_const = c1
    call set_defaults_tracer_read(p_remin_scalef_input, file_varname='p_remin_scalef')
    lignore_driver_ndep = .false.
    ndep_data_type = 'monthly-calendar'
    call set_defaults_tracer_read(nox_flux_monthly_input, file_varname='nox_flux')
    call set_defaults_tracer_read(nhy_flux_monthly_input, file_varname='nhy_flux')
    ndep_shr_stream_year_first = 1
    ndep_shr_stream_year_last  = 1
    ndep_shr_stream_year_align = 1
    ndep_shr_stream_file       = 'unknown'
    ndep_shr_stream_scale_factor = c1
    riv_flux_shr_stream_file         = 'unknown'
    riv_flux_shr_stream_year_first   = 1900
    riv_flux_shr_stream_year_last    = 1900
    riv_flux_shr_stream_year_align   = 1900
    riv_flux_din_file_varname        = 'din_riv_flux'
    riv_flux_din_scale_factor        = c1
    riv_flux_dip_file_varname        = 'dip_riv_flux'
    riv_flux_dip_scale_factor        = c1
    riv_flux_don_file_varname        = 'don_riv_flux'
    riv_flux_don_scale_factor        = c1
    riv_flux_dop_file_varname        = 'dop_riv_flux'
    riv_flux_dop_scale_factor        = c1
    riv_flux_dsi_file_varname        = 'dsi_riv_flux'
    riv_flux_dsi_scale_factor        = c1
    riv_flux_dfe_file_varname        = 'dfe_riv_flux'
    riv_flux_dfe_scale_factor        = c1
    riv_flux_dic_file_varname        = 'dic_riv_flux'
    riv_flux_dic_scale_factor        = c1
    riv_flux_alk_file_varname        = 'alk_riv_flux'
    riv_flux_alk_scale_factor        = c1
    riv_flux_doc_file_varname        = 'doc_riv_flux'
    riv_flux_doc_scale_factor        = c1
    liron_patch              = .false.
    iron_patch_flux_filename = 'unknown_iron_patch_filename'
    iron_patch_month         = 1
    atm_co2_opt          = 'const'
    atm_co2_const        = 280.0_r8
    box_atm_co2_init_val = 280.0_r8
    atm_alt_co2_opt      = 'const'
    atm_alt_co2_const    = 280.0_r8
    ciso_atm_d13c_opt                       = 'const'
    ciso_atm_d13c_const                     = -6.610_r8
    ciso_atm_d13c_filename                  = 'unknown'
    ciso_atm_d14c_opt                       = 'lat_bands'
    ciso_atm_d14c_const                     = c0
    ciso_atm_d14c_lat_band_vals(:)          = (/ -2.3_r8, -4.0_r8, -5.8_r8 /)
    ciso_atm_d14c_filename                  = 'unknown'
    ciso_atm_model_year                     = 1
    ciso_atm_data_year                      = 1

    restorable_tracer_names(:)    = ''
    restore_data_filenames(:)     = ''
    restore_data_file_varnames(:) = ''
    restore_year_first(:)         = 1
    restore_year_last(:)          = 1
    restore_year_align(:)         = 1
    restore_scale_factor(:)       = c1

    restore_inv_tau_opt   = 'const'
    restore_inv_tau_const = c0

    call set_defaults_tracer_read(restore_inv_tau_input, file_varname='RESTORE_INV_TAU_MARGINAL_SEA_ONLY')

    surf_avg_alk_const       = 2225.0_r8
    surf_avg_dic_const       = 1944.0_r8
    surf_avg_di13c_const     = 1944.0_r8
    surf_avg_di14c_const     = 1944.0_r8

    iron_frac_in_atm_fine_dust   = 0.035_r8
    iron_frac_in_atm_coarse_dust = 0.035_r8
    iron_frac_in_seaice_dust     = 0.035_r8
    iron_frac_in_atm_bc          = 0.06_r8
    iron_frac_in_seaice_bc       = 0.06_r8

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

    ! Consistency checks
    ! Diagnostic ATM_CO2 requested from driver but driver doesn't provide it?
    if ((trim(atm_co2_opt) .eq. 'drv_diag') .and. (.not. ldriver_has_atm_co2_diag)) then
        call document(subname, "ERROR: atm_co2_opt is requesting diagnostic CO2 from coupler but coupler is not providing it!")
        call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    ! Prognostic ATM_CO2 requested from driver but driver doesn't provide it?
    if ((trim(atm_co2_opt) .eq. 'drv_prog') .and. (.not. ldriver_has_atm_co2_prog)) then
        call document(subname, "ERROR: atm_co2_opt is requesting prognostic CO2 from coupler but coupler is not providing it!")
        call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    ! NDEP requested from driver but driver doesn't provide it?
    if ((ndep_data_type.eq.'driver') .and. (.not. ldriver_has_ndep)) then
        call document(subname, "ERROR: ndep_data_type is 'driver' but coupler is not providing nitrogen deposition!")
        call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    ! NDEP provided by driver but not read and not intentionally
    ! ignored?
    if ((.not. lignore_driver_ndep) .and. (ldriver_has_ndep) .and. &
        (ndep_data_type.ne.'driver')) then
        ! If either named field index is greater than 0, then either
        ! ndep_data_type must be driver or lignore_driver_ndep must be .true.
        call document(subname, "ERROR: coupler is providing nitrogen deposition but ndep_data_type is not 'driver'")
        call document(subname, "To use the coupler values, set ndep_data_type = 'driver'")
        call document(subname, 'Otherwise set lignore_driver_ndep = .true.')
        call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    if (ciso_on) then
       call ciso_init_atm_D13_D14()
    end if

    ! Set surf_avg for all tracers
    surf_avg(:) = c0
    if (any((/dic_ind, dic_alt_co2_ind, alk_ind, alk_alt_co2_ind/).eq.0)) then
      call document(subname, 'dic_ind, alk_ind, dic_alt_co2_ind, and alk_alt_co2_ind must be non-zero')
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    surf_avg(dic_ind)         = surf_avg_dic_const
    surf_avg(dic_alt_co2_ind) = surf_avg_dic_const
    surf_avg(alk_ind)         = surf_avg_alk_const
    surf_avg(alk_alt_co2_ind) = surf_avg_alk_const

    if (ciso_on) then
       if (any((/di13c_ind, di14c_ind/).eq.0)) then
         call document(subname, 'di13c_ind and di14c_ind must be non-zero')
         call exit_POP(sigAbort, 'Stopping in ' // subname)
       end if
       surf_avg(di13c_ind) = surf_avg_di13c_const
       surf_avg(di14c_ind) = surf_avg_di14c_const
    end if

    call get_timer(ecosys_pre_sflux_timer               , 'ECOSYS_PRE_SFLUX'               , 1, distrb_clinic%nprocs)
    call get_timer(ecosys_riv_flux_strdata_create_timer , 'ecosys_riv_flux_strdata_create' , 1, distrb_clinic%nprocs)
    call get_timer(ecosys_surface_strdata_create_timer  , 'ecosys_surface_strdata_create'  , 1, distrb_clinic%nprocs)
    call get_timer(ecosys_interior_strdata_create_timer , 'ecosys_interior_strdata_create' , 1, distrb_clinic%nprocs)
    call get_timer(ecosys_riv_flux_strdata_advance_timer, 'ecosys_riv_flux_strdata_advance', 1, distrb_clinic%nprocs)
    call get_timer(ecosys_surface_strdata_advance_timer , 'ecosys_surface_strdata_advance' , 1, distrb_clinic%nprocs)
    call get_timer(ecosys_interior_strdata_advance_timer, 'ecosys_interior_strdata_advance', 1, distrb_clinic%nprocs)

    !--------------------------------------------------------------------------
    !  River flux forcing
    !--------------------------------------------------------------------------

    call init_riv_flux_forcing_fields(ciso_on, &
         riv_flux_shr_stream_file, riv_flux_shr_stream_year_first, &
         riv_flux_shr_stream_year_last, riv_flux_shr_stream_year_align, &
         riv_flux_din_file_varname, riv_flux_din_scale_factor, &
         riv_flux_dip_file_varname, riv_flux_dip_scale_factor, &
         riv_flux_don_file_varname, riv_flux_don_scale_factor, &
         riv_flux_dop_file_varname, riv_flux_dop_scale_factor, &
         riv_flux_dsi_file_varname, riv_flux_dsi_scale_factor, &
         riv_flux_dfe_file_varname, riv_flux_dfe_scale_factor, &
         riv_flux_dic_file_varname, riv_flux_dic_scale_factor, &
         riv_flux_alk_file_varname, riv_flux_alk_scale_factor, &
         riv_flux_doc_file_varname, riv_flux_doc_scale_factor, &
         lhas_riv_flux)

    !--------------------------------------------------------------------------
    !  Surface flux forcing
    !--------------------------------------------------------------------------

    ! Set up where surface flux forcing data will come from
    fice_file_loc%input             = gas_flux_fice
    xkw_file_loc%input              = gas_flux_ws
    ap_file_loc%input               = gas_flux_ap
    dust_flux_file_loc%input        = dust_flux_input
    iron_flux_file_loc%input        = iron_flux_input
    nox_flux_monthly_file_loc%input = nox_flux_monthly_input
    nhy_flux_monthly_file_loc%input = nhy_flux_monthly_input

    allocate(surface_flux_forcings(size(marbl_req_surface_flux_forcings)))

    ! allocating to zero size eases the implementation of growing the array
    allocate(surface_strdata_inputlist_ptr(0))

    do n=1,size(surface_flux_forcings)
      marbl_varname = marbl_req_surface_flux_forcings(n)%metadata%varname
      units         = marbl_req_surface_flux_forcings(n)%metadata%field_units
      select case (trim(marbl_req_surface_flux_forcings(n)%metadata%varname))
        case ('d13c')
          d13c_ind = n
          call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
               marbl_varname=marbl_varname, field_units=units,                           &
               driver_varname='d13C', rank=2, id=n)

        case ('d14c')
          d14c_ind = n
          call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
               marbl_varname=marbl_varname, field_units=units,                           &
               driver_varname='D14C', rank=2, id=n)

        case ('u10_sqr')
          u10sqr_ind = n
          call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
               marbl_varname=marbl_varname, field_units=units,                           &
               driver_varname='U10_SQR', rank=2, id=n)

        case ('sst')
          sst_ind = n
          call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
               marbl_varname=marbl_varname, field_units=units,                           &
               driver_varname='SST', rank=2, id=n)

        case ('sss')
          sss_ind = n
          call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
               marbl_varname=marbl_varname, field_units=units,                           &
               driver_varname='SSS', rank=2, id=n)

        case ('xco2')
          if (trim(atm_co2_opt).eq.'const') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='const', &
                 marbl_varname=marbl_varname, field_units=units,                        &
                 field_constant=atm_co2_const, rank=2, id=n)
          else if (trim(atm_co2_opt).eq.'drv_prog') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='named_field', &
                 marbl_varname=marbl_varname, field_units=units,                              &
                 named_field='ATM_CO2_PROG', rank=2, id=n)
          else if (trim(atm_co2_opt).eq.'drv_diag') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='named_field', &
                 marbl_varname=marbl_varname, field_units=units,                              &
                 named_field='ATM_CO2_DIAG', rank=2, id=n)
          else if (trim(atm_co2_opt).eq.'box_atm_co2') then
            box_atm_co2_ind = n
            call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
                 marbl_varname=marbl_varname, field_units=units,                           &
                 driver_varname='box_atm_co2', rank=2, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(atm_co2_opt),                     &
                 'is not a valid option for atm_co2_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('xco2_alt_co2')
          if (trim(atm_alt_co2_opt).eq.'const') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='const', &
                 marbl_varname=marbl_varname, field_units=units,                        &
                 field_constant=atm_alt_co2_const, rank=2, id=n)
          else if (trim(atm_alt_co2_opt).eq.'box_atm_co2') then
            if (trim(atm_co2_opt).eq.'box_atm_co2') then
              box_atm_co2_dup_ind = n
            else
              box_atm_co2_ind = n
            end if
            call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
                 marbl_varname=marbl_varname, field_units=units,                           &
                 driver_varname='box_atm_co2', rank=2, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(atm_alt_co2_opt),                 &
                 'is not a valid option for atm_alt_co2_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('Ice Fraction')
          ifrac_ind = n
          if (trim(gas_flux_forcing_opt).eq.'drv') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
                 marbl_varname=marbl_varname, field_units=units,                           &
                 driver_varname='ICE Fraction', rank=2, id=n)
          else if (trim(gas_flux_forcing_opt).eq.'file') then
            file_details => fice_file_loc
            call init_monthly_surface_flux_forcing_metadata(file_details)
            call surface_flux_forcings(n)%add_forcing_field(                    &
                 field_source='POP monthly calendar',                                 &
                 marbl_varname=marbl_varname, field_units=units,                      &
                 forcing_calendar_name=file_details, rank=2, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(gas_flux_forcing_opt),            &
                 'is not a valid option for gas_flux_forcing_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('Atmospheric Pressure')
          ap_ind = n
          if (trim(gas_flux_forcing_opt).eq.'drv') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
                 marbl_varname=marbl_varname, field_units=units,                           &
                 driver_varname='AP_FILE_INPUT', rank=2, id=n)
          else if (trim(gas_flux_forcing_opt).eq.'file') then
            file_details => ap_file_loc
            call init_monthly_surface_flux_forcing_metadata(file_details)
            call surface_flux_forcings(n)%add_forcing_field(                    &
                 field_source='POP monthly calendar',                                 &
                 marbl_varname=marbl_varname, field_units=units,                      &
                                 forcing_calendar_name=file_details, rank=2, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(gas_flux_forcing_opt),            &
                 'is not a valid option for gas_flux_forcing_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('Dust Flux')
          dust_dep_ind = n
          if (trim(dust_flux_source).eq.'driver') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
                 marbl_varname=marbl_varname, field_units=units,                           &
                 driver_varname='DUST_FLUX', rank=2, id=n)
          else if (trim(dust_flux_source).eq.'monthly-calendar') then
            file_details => dust_flux_file_loc
            call init_monthly_surface_flux_forcing_metadata(file_details)
            call surface_flux_forcings(n)%add_forcing_field(                    &
                 field_source='POP monthly calendar',                                 &
                 marbl_varname=marbl_varname, field_units=units,                      &
                 forcing_calendar_name=file_details, rank=2, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(dust_flux_source),                &
                 'is not a valid option for dust_flux_source'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('Iron Flux')
          if (trim(iron_flux_source).eq.'driver-derived') then
            bc_dep_ind = n
            call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
                 marbl_varname=marbl_varname, field_units=units,                           &
                 driver_varname='BLACK_CARBON_FLUX', rank=2, id=n)
          else if (trim(iron_flux_source).eq.'monthly-calendar') then
            Fe_dep_ind = n
            file_details => iron_flux_file_loc
            call init_monthly_surface_flux_forcing_metadata(file_details)
            call surface_flux_forcings(n)%add_forcing_field(                    &
                 field_source='POP monthly calendar',                                 &
                 marbl_varname=marbl_varname, field_units=units,                      &
                 forcing_calendar_name=file_details, rank=2, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(iron_flux_source),                &
                 'is not a valid option for iron_flux_source'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('NOx Flux')
          if (trim(ndep_data_type).eq.'shr_stream') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='shr_stream', &
                 strdata_inputlist_ptr=surface_strdata_inputlist_ptr,                        &
                 marbl_varname=marbl_varname, field_units=units,                             &
                 unit_conv_factor=ndep_shr_stream_scale_factor,                              &
                 file_varname='NOy_deposition',                                              &
                 year_first = ndep_shr_stream_year_first,                                    &
                 year_last = ndep_shr_stream_year_last,                                      &
                 year_align = ndep_shr_stream_year_align,                                    &
                 filename = ndep_shr_stream_file,                                            &
                 rank = 2, id = n)
          else if (trim(ndep_data_type).eq.'monthly-calendar') then
            file_details => nox_flux_monthly_file_loc
            call init_monthly_surface_flux_forcing_metadata(file_details)
            call surface_flux_forcings(n)%add_forcing_field(                    &
                 field_source='POP monthly calendar',                                 &
                 marbl_varname=marbl_varname, field_units=units,                      &
                 forcing_calendar_name=file_details, rank=2, id=n)
          else if (trim(ndep_data_type).eq.'driver') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='named_field', &
                 marbl_varname=marbl_varname, field_units=units,                              &
                 named_field='ATM_NOy', rank=2, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(ndep_data_type),                  &
                 'is not a valid option for ndep_data_type'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('NHy Flux')
          if (trim(ndep_data_type).eq.'shr_stream') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='shr_stream', &
                 strdata_inputlist_ptr=surface_strdata_inputlist_ptr,                        &
                 marbl_varname=marbl_varname, field_units=units,                             &
                 unit_conv_factor=ndep_shr_stream_scale_factor,                              &
                 file_varname='NHx_deposition',                                              &
                 year_first = ndep_shr_stream_year_first,                                    &
                 year_last = ndep_shr_stream_year_last,                                      &
                 year_align = ndep_shr_stream_year_align,                                    &
                 filename = ndep_shr_stream_file,                                            &
                 rank = 2, id = n)
          else if (trim(ndep_data_type).eq.'monthly-calendar') then
            file_details => nhy_flux_monthly_file_loc
            call init_monthly_surface_flux_forcing_metadata(file_details)
            call surface_flux_forcings(n)%add_forcing_field(                    &
                 field_source='POP monthly calendar',                                 &
                 marbl_varname=marbl_varname, field_units=units,                      &
                 forcing_calendar_name=file_details, rank=2, id=n)
          else if (trim(ndep_data_type).eq.'driver') then
            call surface_flux_forcings(n)%add_forcing_field(field_source='named_field', &
                 marbl_varname=marbl_varname, field_units=units,                              &
                 named_field='ATM_NHx', rank=2, id=n)
          else
            write(err_msg, "(A,1X,A)") trim(ndep_data_type),                  &
                 'is not a valid option for ndep_data_type'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
          end if

        case ('external C Flux')
          ext_C_flux_ind = n
          call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
               marbl_varname=marbl_varname, field_units=units,                           &
               driver_varname='ext_C_flux', rank=2, id=n)

        case ('external P Flux')
          ext_P_flux_ind = n
          call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
               marbl_varname=marbl_varname, field_units=units,                           &
               driver_varname='ext_P_flux', rank=2, id=n)

        case ('external Si Flux')
          ext_Si_flux_ind = n
          call surface_flux_forcings(n)%add_forcing_field(field_source='internal', &
               marbl_varname=marbl_varname, field_units=units,                           &
               driver_varname='ext_Si_flux', rank=2, id=n)

        case DEFAULT
          write(err_msg, "(A,1X,A)") trim(marbl_req_surface_flux_forcings(n)%metadata%varname), &
                         'is not a valid surface flux forcing field name.'
          call document(subname, err_msg)
          call exit_POP(sigAbort, 'Stopping in ' // subname)
      end select
    end do

    !--------------------------------------------------------------------------

    lbox_atm_co2 = box_atm_co2_ind > 0

    !--------------------------------------------------------------------------
    !  Interior tendency forcing
    !--------------------------------------------------------------------------

    allocate(interior_tendency_forcings(size(marbl_req_interior_tendency_forcings)))

    ! allocating to zero size eases the implementation of growing the array
    allocate(interior_strdata_inputlist_ptr(0))

    do n=1,size(interior_tendency_forcings)
      marbl_varname = marbl_req_interior_tendency_forcings(n)%metadata%varname
      units = marbl_req_interior_tendency_forcings(n)%metadata%field_units

      var_processed = .false.
      ! Check to see if this forcing field is tracer restoring
      if (index(marbl_varname,'Restoring Field').gt.0) then
        tracer_name = trim(marbl_varname(1:scan(marbl_varname,' ')))
        do m=1,marbl_tracer_cnt
          if (trim(tracer_name).eq.trim(restorable_tracer_names(m))) then
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
            call interior_tendency_forcings(n)%add_forcing_field(                &
                       field_source='shr_stream',                             &
                       strdata_inputlist_ptr=interior_strdata_inputlist_ptr,  &
                       marbl_varname=marbl_varname, field_units=units,        &
                       filename=restore_data_filenames(m),                    &
                       file_varname=restore_data_file_varnames(m),            &
                       year_first=restore_year_first(m),                      &
                       year_last=restore_year_last(m),                        &
                       year_align=restore_year_align(m),                      &
                       unit_conv_factor=restore_scale_factor(m),              &
                       rank=3, dim3_len=km, id=n)
            var_processed = .true.
            exit
          end if
        end do
      end if

      ! Check to see if this forcing field is a restoring time scale
      if (index(marbl_varname,'Restoring Inverse Timescale').gt.0) then
        select case (trim(restore_inv_tau_opt))
          case('const')
            call interior_tendency_forcings(n)%add_forcing_field(                &
                       field_source='const',                                  &
                       marbl_varname=marbl_varname, field_units=units,        &
                       field_constant=restore_inv_tau_const,                  &
                       rank=3, dim3_len=km, id=n)
          case('file_time_invariant')
            call interior_tendency_forcings(n)%add_forcing_field(                &
                       field_source='file_time_invariant',                    &
                       marbl_varname=marbl_varname, field_units=units,        &
                       filename=restore_inv_tau_input%filename,               &
                       file_varname=restore_inv_tau_input%file_varname,       &
                       unit_conv_factor=restore_inv_tau_input%scale_factor,   &
                       rank=3, dim3_len=km, id=n)
          case DEFAULT
            write(err_msg, "(A,1X,A)") trim(restore_inv_tau_opt),             &
                 'is not a valid option for restore_inv_tau_opt'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
        end select
        var_processed = .true.
      end if

      if (.not.var_processed) then
        select case (trim(marbl_req_interior_tendency_forcings(n)%metadata%varname))
          case ('Dust Flux')
            dustflux_ind = n
            call interior_tendency_forcings(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='dust_flux', rank=2, id=n)
          case ('PAR Column Fraction')
            PAR_col_frac_ind = n
            call interior_tendency_forcings(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='PAR_col_frac', rank=3, dim3_len=mcog_nbins,  &
                          ldim3_is_depth=.false., id=n)
          case ('Surface Shortwave')
            surf_shortwave_ind = n
            call interior_tendency_forcings(n)%add_forcing_field(field_source='internal',  &
                          marbl_varname=marbl_varname, field_units=units,               &
                          driver_varname='surf_shortwave', rank=3, dim3_len=mcog_nbins, &
                          ldim3_is_depth=.false., id=n)
          case ('Potential Temperature')
            potemp_ind = n
            call interior_tendency_forcings(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='temperature', rank=3, dim3_len=km, id=n)
          case ('Salinity')
            salinity_ind = n
            call interior_tendency_forcings(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='salinity', rank=3, dim3_len=km, id=n)
          case ('Pressure')
            pressure_ind = n
            call interior_tendency_forcings(n)%add_forcing_field(field_source='internal', &
                          marbl_varname=marbl_varname, field_units=units,              &
                          driver_varname='pressure', rank=3, dim3_len=km, id=n)
          case ('Iron Sediment Flux')
            fesedflux_ind = n
            call interior_tendency_forcings(n)%add_forcing_field(                &
                          field_source='file_time_invariant',                 &
                          marbl_varname=marbl_varname, field_units=units,     &
                          filename=fesedflux_input%filename,                  &
                          file_varname=fesedflux_input%file_varname,          &
                          unit_conv_factor=fesedflux_input%scale_factor,      &
                          rank=3, dim3_len=km, id=n)
          case ('O2 Consumption Scale Factor')
            select case (trim(o2_consumption_scalef_opt))
            case ('const')
              call interior_tendency_forcings(n)%add_forcing_field(field_source='const', &
                   marbl_varname=marbl_varname, field_units=units, &
                   field_constant=o2_consumption_scalef_const, rank=3, dim3_len=km, id=n)
            case ('file_time_invariant')
              call interior_tendency_forcings(n)%add_forcing_field( &
                   field_source='file_time_invariant', &
                   marbl_varname=marbl_varname, field_units=units, &
                   filename=o2_consumption_scalef_input%filename, &
                   file_varname=o2_consumption_scalef_input%file_varname, &
                   unit_conv_factor=o2_consumption_scalef_input%scale_factor, &
                   rank=3, dim3_len=km, id=n)
            case default
              call document(subname, 'unknown o2_consumption_scalef_opt', o2_consumption_scalef_opt)
              call exit_POP(sigAbort, 'Stopping in ' // subname)
            end select
          case ('Particulate Remin Scale Factor')
            select case (trim(p_remin_scalef_opt))
            case ('const')
              call interior_tendency_forcings(n)%add_forcing_field(field_source='const', &
                   marbl_varname=marbl_varname, field_units=units, &
                   field_constant=p_remin_scalef_const, rank=3, dim3_len=km, id=n)
            case ('file_time_invariant')
              call interior_tendency_forcings(n)%add_forcing_field( &
                   field_source='file_time_invariant', &
                   marbl_varname=marbl_varname, field_units=units, &
                   filename=p_remin_scalef_input%filename, &
                   file_varname=p_remin_scalef_input%file_varname, &
                   unit_conv_factor=p_remin_scalef_input%scale_factor, &
                   rank=3, dim3_len=km, id=n)
            case default
              call document(subname, 'unknown p_remin_scalef_opt', p_remin_scalef_opt)
              call exit_POP(sigAbort, 'Stopping in ' // subname)
            end select
          case DEFAULT
            write(err_msg, "(A,1X,A)") trim(marbl_req_interior_tendency_forcings(n)%metadata%varname), &
                           'is not a valid interior tendency forcing field name.'
            call document(subname, err_msg)
            call exit_POP(sigAbort, 'Stopping in ' // subname)
        end select
      end if

    end do ! do n=1,size(interior_tendency_forcings)

    if ((bc_dep_ind.ne.0).and.(dust_dep_ind.eq.0)) then
      call document(subname, "If deriving iron flux, must provide dust flux!")
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    call set_time_invariant_forcing_data(surface_flux_forcings)
    call set_time_invariant_forcing_data(interior_tendency_forcings)

    call read_monthly_calendar_forcing_data(surface_flux_forcings, land_mask)

    call forcing_init_post_processing(land_mask)

    call mask_time_invariant_forcing_data(surface_flux_forcings, land_mask)
    call mask_time_invariant_forcing_data(interior_tendency_forcings, land_mask)

  end subroutine ecosys_forcing_init

  !*****************************************************************************

  subroutine init_riv_flux_forcing_fields(ciso_on, &
       riv_flux_shr_stream_file, riv_flux_shr_stream_year_first, &
       riv_flux_shr_stream_year_last, riv_flux_shr_stream_year_align, &
       riv_flux_din_file_varname, riv_flux_din_scale_factor, &
       riv_flux_dip_file_varname, riv_flux_dip_scale_factor, &
       riv_flux_don_file_varname, riv_flux_don_scale_factor, &
       riv_flux_dop_file_varname, riv_flux_dop_scale_factor, &
       riv_flux_dsi_file_varname, riv_flux_dsi_scale_factor, &
       riv_flux_dfe_file_varname, riv_flux_dfe_scale_factor, &
       riv_flux_dic_file_varname, riv_flux_dic_scale_factor, &
       riv_flux_alk_file_varname, riv_flux_alk_scale_factor, &
       riv_flux_doc_file_varname, riv_flux_doc_scale_factor, &
       lhas_riv_flux)

    logical          , intent(in)  :: ciso_on
    character(len=*) , intent(in)  :: riv_flux_shr_stream_file
    integer(int_kind), intent(in)  :: riv_flux_shr_stream_year_first
    integer(int_kind), intent(in)  :: riv_flux_shr_stream_year_last
    integer(int_kind), intent(in)  :: riv_flux_shr_stream_year_align
    character(len=*) , intent(in)  :: riv_flux_din_file_varname
    real(r8)         , intent(in)  :: riv_flux_din_scale_factor
    character(len=*) , intent(in)  :: riv_flux_dip_file_varname
    real(r8)         , intent(in)  :: riv_flux_dip_scale_factor
    character(len=*) , intent(in)  :: riv_flux_don_file_varname
    real(r8)         , intent(in)  :: riv_flux_don_scale_factor
    character(len=*) , intent(in)  :: riv_flux_dop_file_varname
    real(r8)         , intent(in)  :: riv_flux_dop_scale_factor
    character(len=*) , intent(in)  :: riv_flux_dsi_file_varname
    real(r8)         , intent(in)  :: riv_flux_dsi_scale_factor
    character(len=*) , intent(in)  :: riv_flux_dfe_file_varname
    real(r8)         , intent(in)  :: riv_flux_dfe_scale_factor
    character(len=*) , intent(in)  :: riv_flux_dic_file_varname
    real(r8)         , intent(in)  :: riv_flux_dic_scale_factor
    character(len=*) , intent(in)  :: riv_flux_alk_file_varname
    real(r8)         , intent(in)  :: riv_flux_alk_scale_factor
    character(len=*) , intent(in)  :: riv_flux_doc_file_varname
    real(r8)         , intent(in)  :: riv_flux_doc_scale_factor
    logical          , intent(out) :: lhas_riv_flux(:) ! true if a tracer has a riverine flux

    !--------------------------------------------------------------------------
    !  local variables
    !--------------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:init_riv_flux_forcing_fields'

    character(len=*), parameter :: riv_flux_components(9) = (/ &
         'din', 'dip', 'don', 'dop', 'dsi', 'dfe', 'dic', 'alk', 'doc' /)

    real(r8)            :: scale_factor
    character(char_len) :: file_varname, marbl_varname, units
    integer             :: n

    !--------------------------------------------------------------------------

    lhas_riv_flux(:) = .false.

    ! allocating to zero size eases the implementation of growing the array
    allocate(riv_flux_strdata_inputlist_ptr(0))

    if (trim(riv_flux_shr_stream_file) == 'unknown') then
      allocate(riv_flux_forcing_fields(0))
      return
    endif

    allocate(riv_flux_forcing_fields(size(riv_flux_components)))

    do n = 1, size(riv_flux_components)
      select case (trim(riv_flux_components(n)))
        case ('din')
          din_riv_flux_ind = n
          file_varname = riv_flux_din_file_varname
          scale_factor = riv_flux_din_scale_factor
          lhas_riv_flux(no3_ind) = .true.
        case ('dip')
          dip_riv_flux_ind = n
          file_varname = riv_flux_dip_file_varname
          scale_factor = riv_flux_dip_scale_factor
          lhas_riv_flux(po4_ind) = .true.
        case ('don')
          don_riv_flux_ind = n
          file_varname = riv_flux_don_file_varname
          scale_factor = riv_flux_don_scale_factor
          lhas_riv_flux(don_ind) = .true.
          lhas_riv_flux(donr_ind) = .true.
        case ('dop')
          dop_riv_flux_ind = n
          file_varname = riv_flux_dop_file_varname
          scale_factor = riv_flux_dop_scale_factor
          lhas_riv_flux(dop_ind) = .true.
          lhas_riv_flux(dopr_ind) = .true.
        case ('dsi')
          dsi_riv_flux_ind = n
          file_varname = riv_flux_dsi_file_varname
          scale_factor = riv_flux_dsi_scale_factor
          lhas_riv_flux(sio3_ind) = .true.
        case ('dfe')
          dfe_riv_flux_ind = n
          file_varname = riv_flux_dfe_file_varname
          scale_factor = riv_flux_dfe_scale_factor
          lhas_riv_flux(fe_ind) = .true.
        case ('dic')
          dic_riv_flux_ind = n
          file_varname = riv_flux_dic_file_varname
          scale_factor = riv_flux_dic_scale_factor
          lhas_riv_flux(dic_ind) = .true.
          lhas_riv_flux(dic_alt_co2_ind) = .true.
          if (ciso_on) then
            lhas_riv_flux(di13c_ind) = .true.
            lhas_riv_flux(di14c_ind) = .true.
          endif
        case ('alk')
          alk_riv_flux_ind = n
          file_varname = riv_flux_alk_file_varname
          scale_factor = riv_flux_alk_scale_factor
          lhas_riv_flux(alk_ind) = .true.
          lhas_riv_flux(alk_alt_co2_ind) = .true.
        case ('doc')
          doc_riv_flux_ind = n
          file_varname = riv_flux_doc_file_varname
          scale_factor = riv_flux_doc_scale_factor
          lhas_riv_flux(doc_ind) = .true.
          lhas_riv_flux(docr_ind) = .true.
          if (ciso_on) then
            lhas_riv_flux(do13ctot_ind) = .true.
            lhas_riv_flux(do14ctot_ind) = .true.
          endif
        case default
          call document(subname, 'unhandled riv_flux file_varname ', file_varname)
          call exit_POP(sigAbort, 'Stopping in ' // subname)
      end select

      marbl_varname = 'none'
      units         = 'nmol/cm^2/s'

      call riv_flux_forcing_fields(n)%add_forcing_field(field_source='shr_stream',    &
                            strdata_inputlist_ptr=riv_flux_strdata_inputlist_ptr,     &
                            marbl_varname=marbl_varname, field_units=units,           &
                            filename=riv_flux_shr_stream_file,                        &
                            file_varname=file_varname, unit_conv_factor=scale_factor, &
                            year_first=riv_flux_shr_stream_year_first,                &
                            year_last=riv_flux_shr_stream_year_last,                  &
                            year_align=riv_flux_shr_stream_year_align,                &
                            taxMode='extend',                                         &
                            rank=2, id=n)
    end do

  end subroutine init_riv_flux_forcing_fields

  !*****************************************************************************

  subroutine set_time_invariant_forcing_data(forcing_fields_in)

    use passive_tracer_tools, only : read_field

    type(forcing_fields_type), intent(inout) :: forcing_fields_in(:)

    !--------------------------------------------------------------------------
    !  local variables
    !--------------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:set_time_invariant_forcing_data'
    integer :: n

    !--------------------------------------------------------------------------
    !  Set any constant forcing fields and
    !  read any forcing fields that are time invariant
    !--------------------------------------------------------------------------

    do n=1, size(forcing_fields_in)
      associate (forcing_field => forcing_fields_in(n), &
                 metadata      => forcing_fields_in(n)%metadata)
        if (.not. metadata%ltime_varying) then
          select case (trim(metadata%field_source))
            case ('const','zero')
              if (forcing_field%rank == 2) then
                forcing_field%field_0d = metadata%field_constant_info%field_constant
              else
                forcing_field%field_1d = metadata%field_constant_info%field_constant
              end if
            case ('file_time_invariant')
              if (forcing_field%rank == 2) then
                call read_field('nc', metadata%field_file_info%filename,   &
                                metadata%field_file_info%file_varname,     &
                                forcing_field%field_0d)
              else
                call read_field('nc', metadata%field_file_info%filename,   &
                                metadata%field_file_info%file_varname,     &
                                forcing_field%field_1d)
              end if
            case default
              call document(subname, 'unhandled time_invariant field_source ', trim(metadata%field_source))
              call exit_POP(sigAbort, 'Stopping in ' // subname)
          end select
        end if
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
    !  For POP monthly calendar surface flux forcing fields, read all
    !  12 months of data
    !--------------------------------------------------------------------------

    do n=1, size(forcing_fields_in)
      associate (metadata => forcing_fields_in(n)%metadata)
        if (trim(metadata%field_source).eq.'POP monthly calendar') then
           file => metadata%field_monthly_calendar_info%forcing_calendar_name

           allocate(work_read(nx_block, ny_block, 12, nblocks_clinic))
           if (trim(file%input%filename) == 'unknown') then
              file%input%filename = gas_flux_forcing_file
           end if
           if (trim(file%input%filename) /= 'none') then
              allocate(file%data(nx_block, ny_block, nblocks_clinic, 1, 12))
              call read_field(file%input%file_fmt, file%input%filename, file%input%file_varname, work_read)
              do iblock=1, nblocks_clinic
                 do m=1, 12
                    file%data(:,:,iblock,1,m) = work_read(:,:,m,iblock)
                    where (.not. land_mask(:,:,iblock)) file%data(:,:,iblock,1,m) = c0
                    file%data(:,:,iblock,1,m) = file%data(:,:,iblock,1,m) * file%input%scale_factor
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
           if (n == Fe_dep_ind .and. liron_patch) then
              allocate(iron_patch_flux(nx_block, ny_block, nblocks_clinic))
              call read_field(file%input%file_fmt, file%input%filename, iron_patch_flux_filename, iron_patch_flux)
              do iblock=1, nblocks_clinic
                 do m=1, 12
                    where (.not. land_mask(:,:,iblock)) iron_patch_flux(:,:,iblock) = c0
                    file%data(:,:,iblock,1,m) = iron_patch_flux(:,:,iblock) * file%input%scale_factor
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

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character (*), parameter :: subname = 'ecosys_forcing_mod:forcing_init_post_processing'

    integer :: iblock, k, n ! loop indices

    real (r8), allocatable, target :: feventflux(:,:,:,:) !  Fe from vents

    do iblock = 1, nblocks_clinic

      !-----------------------------------------------------------------------
      !  apply subsurface adjustment to fesedflux
      !-----------------------------------------------------------------------

      if (fesedflux_ind > 0) then
        associate (fesedflux => interior_tendency_forcings(fesedflux_ind)%field_1d)
          call add_subsurf_to_bottom(land_mask(:,:,iblock), fesedflux(:,:,:,iblock), iblock)
        end associate
      endif

      !-----------------------------------------------------------------------
      ! apply unit_conv_factor to all non-time-varying interior and surface flux forcing fields
      !-----------------------------------------------------------------------

      do n=1,size(surface_flux_forcings)
        if (.not. surface_flux_forcings(n)%metadata%ltime_varying) then
           call apply_unit_conv_factor(land_mask(:,:,iblock), surface_flux_forcings(n), iblock)
        end if
      end do

      do n=1,size(interior_tendency_forcings)
        if (.not. interior_tendency_forcings(n)%metadata%ltime_varying) then
           call apply_unit_conv_factor(land_mask(:,:,iblock), interior_tendency_forcings(n), iblock)
        end if
      end do

    enddo

    !-----------------------------------------------------------------------
    ! add feventflux*scale_factor to fesedflux
    ! abort if feventflux is provided and MARBL is not asking for fesedflux
    ! apply subsurface adjustment after reading in field
    !-----------------------------------------------------------------------

    if (feventflux_input%filename /= 'unknown') then
      if (fesedflux_ind == 0) then
        call document(subname, 'feventflux_input%filename', feventflux_input%filename)
        call exit_POP(sigAbort, 'feventflux is specified, but fesedflux is not requested by MARBL')
      else
        associate (fesedflux => interior_tendency_forcings(fesedflux_ind)%field_1d)
          allocate(feventflux(nx_block, ny_block, km, nblocks_clinic))
          call read_field(feventflux_input%file_fmt, &
               feventflux_input%filename, &
               feventflux_input%file_varname, &
               feventflux)
          do iblock = 1, nblocks_clinic
            call add_subsurf_to_bottom(land_mask(:,:,iblock), feventflux(:,:,:,iblock), iblock)
            do k = 1, km
              where (land_mask(:,:,iblock) .and. (k.le.KMT(:,:,iblock)))
                fesedflux(:,:,k,iblock) = fesedflux(:,:,k,iblock) &
                     + feventflux(:,:,k,iblock) * feventflux_input%scale_factor
              end where
            enddo
          enddo
        end associate
      endif
    endif

  end subroutine forcing_init_post_processing

  !*****************************************************************************

  subroutine add_subsurf_to_bottom(land_mask, field, iblock)

    !  add subsurface values to bottom layer, to accomodate overflow pop-ups

    use grid, only : KMT

    logical               , intent(in)    :: land_mask(:,:)
    real(kind=r8)         , intent(inout) :: field(:,:,:) ! (nx, ny, km)
    integer(kind=int_kind), intent(in)    :: iblock

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer(kind=int_kind) :: i, j ! loop indices

    do j=1, ny_block
      do i=1, nx_block
        if (land_mask(i,j) .and. KMT(i,j,iblock) > 0 .and. KMT(i,j,iblock) < km) then
          field(i,j,KMT(i,j,iblock)) = field(i,j,KMT(i,j,iblock)) + sum(field(i,j,KMT(i,j,iblock)+1:km))
        endif
      enddo
    enddo

  end subroutine add_subsurf_to_bottom

  !*****************************************************************************

  subroutine apply_unit_conv_factor(land_mask, forcing_field, iblock)

    use grid, only : KMT

    logical                  , intent(in)    :: land_mask(:,:)
    type(forcing_fields_type), intent(inout) :: forcing_field
    integer(kind=int_kind)   , intent(in)    :: iblock

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer(kind=int_kind) :: k ! loop index

    associate (metadata => forcing_field%metadata)
      if (metadata%unit_conv_factor == c1) return

      if (forcing_field%rank == 2) then
        where (land_mask)
          forcing_field%field_0d(:,:,iblock) = &
               metadata%unit_conv_factor * forcing_field%field_0d(:,:,iblock)
        endwhere
      else
        if (forcing_field%ldim3_is_depth) then
          do k = 1, km
            where (land_mask .and. (k .le. KMT(:,:,iblock)))
              forcing_field%field_1d(:,:,k,iblock) = &
                   metadata%unit_conv_factor * forcing_field%field_1d(:,:,k,iblock)
            endwhere
          end do
        else
          do k = 1, size(forcing_field%field_1d,3)
            where (land_mask)
              forcing_field%field_1d(:,:,k,iblock) = &
                   metadata%unit_conv_factor * forcing_field%field_1d(:,:,k,iblock)
            endwhere
          end do
        endif
      endif
    end associate

  end subroutine apply_unit_conv_factor

  !*****************************************************************************

  subroutine mask_time_invariant_forcing_data(forcing_fields_in, land_mask)

    use grid, only : KMT

    type(forcing_fields_type), intent(inout) :: forcing_fields_in(:)
    logical,                   intent(in)    :: land_mask(:,:,:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer(kind=int_kind) :: n, iblock, k ! loop indices

    do n=1,size(forcing_fields_in)
      associate (forcing_field => forcing_fields_in(n), &
                 metadata      => forcing_fields_in(n)%metadata)
        if (.not. metadata%ltime_varying) then
          do iblock=1,nblocks_clinic
            if (forcing_field%rank == 2) then
              where (.not.land_mask(:,:,iblock))
                forcing_field%field_0d(:,:,iblock) = c0
              endwhere
            else
              if (forcing_field%ldim3_is_depth) then
                do k = 1, km
                  where (.not.land_mask(:,:,iblock) .or. (k.gt.KMT(:,:,iblock)))
                    forcing_field%field_1d(:,:,k,iblock) = c0
                  end where
                end do
              else
                do k = 1, size(forcing_field%field_1d,3)
                  where (.not.land_mask(:,:,iblock))
                    forcing_field%field_1d(:,:,k,iblock) = c0
                  end where
                end do
              end if
            end if
          end do
        end if
      end associate
    end do

  end subroutine mask_time_invariant_forcing_data

  !*****************************************************************************

  subroutine ecosys_forcing_set_interior_time_varying_forcing_data( &
         FRACR_BIN, QSW_RAW_BIN, QSW_BIN, ecosys_qsw_distrb_const, land_mask)

    use blocks                , only : get_block
    use constants             , only : p5, salt_to_ppt
    use domain                , only : blocks_clinic
    use grid                  , only : KMT
    use mcog                  , only : mcog_nbins
    use prognostic            , only : TRACER, oldtime, curtime
    use state_mod             , only : ref_pressure
    use strdata_interface_mod , only : POP_strdata_create, POP_strdata_advance

    real(r8), dimension(nx_block, ny_block, mcog_nbins, nblocks_clinic), intent(in) :: FRACR_BIN
    real(r8), dimension(nx_block, ny_block, mcog_nbins, nblocks_clinic), intent(in) :: QSW_RAW_BIN
    real(r8), dimension(nx_block, ny_block, mcog_nbins, nblocks_clinic), intent(in) :: QSW_BIN
    logical,                                                             intent(in) :: ecosys_qsw_distrb_const
    logical,                                                             intent(in) :: land_mask(:,:,:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character (*), parameter :: subname = 'ecosys_forcing_mod:ecosys_forcing_set_interior_time_varying_forcing_data'
    logical(log_kind)        :: first_call = .true.
    integer(int_kind)        :: field_index, i, j, k, iblock, n   ! loop indices
    integer(int_kind)        :: n0(nblocks_clinic)
    integer(int_kind)        :: stream_index                ! index into interior_strdata_inputlist_ptr array
    integer(int_kind)        :: var_ind                     ! var index in interior_strdata_inputlist_ptr entry
    type(block)              :: this_block                  ! block info for the current block

    !-----------------------------------------------------------------------
    ! Initialize interior_strdata_inputlist_ptr entries (only once)
    !-----------------------------------------------------------------------

    if (first_call) then
      call timer_start(ecosys_interior_strdata_create_timer)
      do n = 1, size(interior_strdata_inputlist_ptr)
        call POP_strdata_create(interior_strdata_inputlist_ptr(n))
      end do
      call timer_stop(ecosys_interior_strdata_create_timer)
      first_call = .false.
    end if ! first_call

    !-----------------------------------------------------------------------
    ! advance interior_strdata_inputlist_ptr entries
    !-----------------------------------------------------------------------

    call timer_start(ecosys_interior_strdata_advance_timer)
    call POP_strdata_advance(interior_strdata_inputlist_ptr(:))
    call timer_stop(ecosys_interior_strdata_advance_timer)

    if ((nblocks_clinic > 0) .and. (size(interior_strdata_inputlist_ptr) > 0)) then
      n0(1) = 0
      do iblock = 1, nblocks_clinic-1
        this_block = get_block(blocks_clinic(iblock), iblock)
        n0(iblock+1) = n0(iblock) + km*(this_block%je-this_block%jb+1)*(this_block%ie-this_block%ib+1)
      enddo
    end if

    !$OMP PARALLEL DO PRIVATE(iblock,this_block,field_index,k,stream_index,var_ind,n,j,i)
    do iblock = 1, nblocks_clinic
      this_block = get_block(blocks_clinic(iblock), iblock)
      do field_index = 1, size(interior_tendency_forcings)
!!!     the following associate construct seems to be incompatible with the OMP directive
!!!     associate (forcing_field => interior_tendency_forcings(field_index), &
!!!                metadata      => interior_tendency_forcings(field_index)%metadata)
          select case (trim(interior_tendency_forcings(field_index)%metadata%field_source))
            case('internal')
              if (field_index .eq. dustflux_ind) then
                interior_tendency_forcings(field_index)%field_0d(:,:,iblock) = dust_flux_in(:,:,iblock)
              else if (field_index .eq. PAR_col_frac_ind) then
                interior_tendency_forcings(field_index)%field_1d(:,:,:,iblock) = FRACR_BIN(:,:,:,iblock)
              else if (field_index .eq. surf_shortwave_ind) then
                if (ecosys_qsw_distrb_const) then
                  interior_tendency_forcings(field_index)%field_1d(:,:,:,iblock) = QSW_RAW_BIN(:,:,:,iblock)
                else
                  interior_tendency_forcings(field_index)%field_1d(:,:,:,iblock) = QSW_BIN(:,:,:,iblock)
                end if
              else if (field_index .eq. potemp_ind) then
                ! --- average 2 time levels into 1 ---
                interior_tendency_forcings(field_index)%field_1d(:,:,:,iblock) = &
                  p5 * (TRACER(:,:,:,1,oldtime,iblock) + TRACER(:,:,:,1,curtime,iblock))
              else if (field_index .eq. salinity_ind) then
                ! --- average 2 time levels into 1, and convert from msu to psu ---
                interior_tendency_forcings(field_index)%field_1d(:,:,:,iblock) = &
                  p5 * (TRACER(:,:,:,2,oldtime,iblock) + TRACER(:,:,:,2,curtime,iblock)) * salt_to_ppt
              else if (field_index .eq. pressure_ind) then
                do k = 1, km
                  interior_tendency_forcings(field_index)%field_1d(:,:,k,iblock) = ref_pressure(k)
                end do
              end if
            case('shr_stream')
              stream_index = interior_tendency_forcings(field_index)%metadata%field_file_info%strdata_inputlist_ind
              var_ind      = interior_tendency_forcings(field_index)%metadata%field_file_info%strdata_var_ind
              n = n0(iblock)
              do k=1,km
                do j = this_block%jb, this_block%je
                  do i = this_block%ib, this_block%ie
                    n = n + 1
                    if (land_mask(i,j,iblock) .and. k .le. KMT(i,j,iblock)) then
                      interior_tendency_forcings(field_index)%field_1d(i,j,k,iblock) = &
                        interior_strdata_inputlist_ptr(stream_index)%sdat%avs(1)%rAttr(var_ind, n)
                    else
                      interior_tendency_forcings(field_index)%field_1d(i,j,k,iblock) = c0
                    endif
                  enddo
                enddo
              enddo
          end select

          if (interior_tendency_forcings(field_index)%metadata%ltime_varying) then
            call apply_unit_conv_factor(land_mask(:,:,iblock), interior_tendency_forcings(field_index), iblock)
          end if
!!!     end associate
      end do
    end do
    !$OMP END PARALLEL DO

    call adjust_interior_time_varying_data()

  end subroutine ecosys_forcing_set_interior_time_varying_forcing_data

  !***********************************************************************

  subroutine adjust_interior_time_varying_data()
    ! This subroutine is empty because there are no interior tendency forcing
    ! fields that need to be modified before being passed to MARBL
  end subroutine adjust_interior_time_varying_data

  !*****************************************************************************

  subroutine ecosys_forcing_set_surface_time_varying_forcing_data( &
       ciso_on,                               &
       land_mask,                             &
       u10_sqr,                               &
       ifrac,                                 &
       press,                                 &
       atm_fine_dust_flux,                    &
       atm_coarse_dust_flux,                  &
       seaice_dust_flux,                      &
       atm_black_carbon_flux,                 &
       seaice_black_carbon_flux,              &
       sst,                                   &
       sss)

    ! !DESCRIPTION:
    !  Compute surface fluxes for ecosys tracer module.

    use POP_HaloMod           , only : POP_HaloUpdate
    use POP_GridHorzMod       , only : POP_gridHorzLocCenter
    use POP_FieldMod          , only : POP_fieldKindScalar
    use POP_ErrorMod          , only : POP_Success
    use domain                , only : POP_haloClinic
    use domain                , only : blocks_clinic
    use blocks                , only : get_block
    use constants             , only : field_loc_center
    use constants             , only : field_type_scalar
    use forcing_tools         , only : interpolate_forcing
    use forcing_tools         , only : update_forcing_data
    use named_field_mod       , only : named_field_get
    use time_management       , only : thour00
    use strdata_interface_mod , only : POP_strdata_create
    use strdata_interface_mod , only : POP_strdata_advance
    use passive_tracer_tools  , only : read_field
    use marbl_constants_mod   , only : molw_Fe

    use ecosys_forcing_saved_state_mod , only : ecosys_forcing_saved_state_get_var_val
    use ecosys_forcing_saved_state_mod , only : box_atm_co2_forcing_saved_state_id

    implicit none

    logical,   intent(in)  :: ciso_on
    logical,   intent(in)  :: land_mask                (nx_block,ny_block,max_blocks_clinic)
    real (r8), intent(in)  :: u10_sqr                  (nx_block,ny_block,max_blocks_clinic) ! 10m wind speed squared (cm/s)**2
    real (r8), intent(in)  :: ifrac                    (nx_block,ny_block,max_blocks_clinic) ! sea ice fraction (non-dimensional)
    real (r8), intent(in)  :: press                    (nx_block,ny_block,max_blocks_clinic) ! sea level atmospheric pressure (dyne/cm**2)
    real (r8), intent(in)  :: atm_fine_dust_flux       (nx_block,ny_block,max_blocks_clinic) ! fine dust flux from atm (g/cm**2/s)
    real (r8), intent(in)  :: atm_coarse_dust_flux     (nx_block,ny_block,max_blocks_clinic) ! coarse dust flux from atm (g/cm**2/s)
    real (r8), intent(in)  :: seaice_dust_flux         (nx_block,ny_block,max_blocks_clinic) ! dust flux from seaice (g/cm**2/s)
    real (r8), intent(in)  :: atm_black_carbon_flux    (nx_block,ny_block,max_blocks_clinic) ! black carbon flux from atm (g/cm**2/s)
    real (r8), intent(in)  :: seaice_black_carbon_flux (nx_block,ny_block,max_blocks_clinic) ! black carbon flux from seaice (g/cm**2/s)
    real (r8), intent(in)  :: sst                      (nx_block,ny_block,max_blocks_clinic) ! sea surface temperature (c)
    real (r8), intent(in)  :: sss                      (nx_block,ny_block,max_blocks_clinic) ! sea surface salinity (psu)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character (*), parameter       :: subname = 'ecosys_forcing_mod:ecosys_forcing_set_surface_time_varying_forcing_data'
    logical   (log_kind)           :: first_call = .true.
    type      (block)              :: this_block                                            ! block info for the current block
    integer   (int_kind)           :: index                                                 ! field index
    integer   (int_kind)           :: riv_flux_ind                                          ! index into riv_flux_forcing_fields
    integer   (int_kind)           :: i, j, iblock, n                                       ! loop indices
    integer   (int_kind)           :: errorCode                                             ! errorCode from HaloUpdate call
    integer   (int_kind)           :: tracer_bndy_loc(1)                                    ! location   for ghost tracer_bndy_type cell updates
    integer   (int_kind)           :: tracer_bndy_type(1)                                   ! field type for ghost tracer_bndy_type cell updates
    character (char_len)           :: tracer_data_names(1)                                  ! short names for input data fields
    real      (r8)                 :: interp_work(nx_block, ny_block, nblocks_clinic, 1)    ! temp array for interpolate_forcing output
    real      (r8)                 :: shr_stream(nx_block, ny_block, nblocks_clinic)
    real      (r8)                 :: d13c(nx_block, ny_block, max_blocks_clinic)           ! atm 13co2 value
    real      (r8)                 :: d14c(nx_block, ny_block, max_blocks_clinic)           ! atm 14co2 value
    type(forcing_monthly_every_ts), pointer :: file
    integer   (int_kind)           :: stream_index                                          ! index into surface_strdata_inputlist_ptr array
    integer   (int_kind)           :: var_ind                                               ! var index in surface_strdata_inputlist_ptr entry

    real      (r8)                 :: atm_fe_bioavail_frac(nx_block, ny_block)
    real      (r8)                 :: seaice_fe_bioavail_frac(nx_block, ny_block)
    real      (r8)                 :: dust_ratio_to_fe_bioavail_frac

    !-----------------------------------------------------------------------

    call timer_start(ecosys_pre_sflux_timer)

    !-----------------------------------------------------------------------
    ! Update carbon isotope atmosphere deltas if appropriate
    !-----------------------------------------------------------------------

    if (ciso_on) then
       call ciso_update_atm_d13C_D14C(land_mask, d13c, d14c)
    end if

    !-----------------------------------------------------------------------
    ! Initialize {riv_flux,surface}_strdata_inputlist_ptr entries (only once)
    !-----------------------------------------------------------------------

    if (first_call) then
      call timer_start(ecosys_surface_strdata_create_timer)
      do n = 1, size(surface_strdata_inputlist_ptr)
        call POP_strdata_create(surface_strdata_inputlist_ptr(n))
      end do
      call timer_stop(ecosys_surface_strdata_create_timer)

      call timer_start(ecosys_riv_flux_strdata_create_timer)
      do n = 1, size(riv_flux_strdata_inputlist_ptr)
        call POP_strdata_create(riv_flux_strdata_inputlist_ptr(n))
      end do
      call timer_stop(ecosys_riv_flux_strdata_create_timer)

      first_call = .false.
    end if ! first_call

    !-----------------------------------------------------------------------
    ! advance {riv_flux,surface}_strdata_inputlist_ptr entries
    !-----------------------------------------------------------------------

    call timer_start(ecosys_surface_strdata_advance_timer)
    call POP_strdata_advance(surface_strdata_inputlist_ptr(:))
    call timer_stop(ecosys_surface_strdata_advance_timer)

    call timer_start(ecosys_riv_flux_strdata_advance_timer)
    call POP_strdata_advance(riv_flux_strdata_inputlist_ptr(:))
    call timer_stop(ecosys_riv_flux_strdata_advance_timer)

    !-----------------------------------------------------------------------
    !  loop through riv_flux forcing fields
    !-----------------------------------------------------------------------

    do index = 1, size(riv_flux_forcing_fields)
       associate (forcing_field => riv_flux_forcing_fields(index), &
                  metadata      => riv_flux_forcing_fields(index)%metadata)
          select case (trim(metadata%field_source))

          !------------------------------------
          case ('shr_stream')
          !------------------------------------

             stream_index = metadata%field_file_info%strdata_inputlist_ind
             var_ind      = metadata%field_file_info%strdata_var_ind

             n = 0
             do iblock = 1, nblocks_clinic
                this_block = get_block(blocks_clinic(iblock), iblock)
                do j = this_block%jb, this_block%je
                   do i = this_block%ib, this_block%ie
                      n = n + 1
                      shr_stream(i,j,iblock) = riv_flux_strdata_inputlist_ptr(stream_index)%sdat%avs(1)%rAttr(var_ind,n)
                   enddo
                enddo
             enddo

             call POP_HaloUpdate(shr_stream, POP_haloClinic, &
                  POP_gridHorzLocCenter, POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
             if (errorCode /= POP_Success) then
                call document(subname, 'error updating halo for shr_stream field')
                call exit_POP(sigAbort, 'Stopping in ' // subname)
             endif

             do iblock = 1, nblocks_clinic
                where (land_mask(:,:,iblock))
                   forcing_field%field_0d(:,:,iblock) = shr_stream(:,:,iblock)
                endwhere
             enddo

          end select

          if (metadata%ltime_varying) then
             do iblock = 1, nblocks_clinic
                call apply_unit_conv_factor(land_mask(:,:,iblock), forcing_field, iblock)
             enddo
          end if
       end associate
    end do  ! index

    !-----------------------------------------------------------------------
    !  loop through surface flux forcing fields
    !-----------------------------------------------------------------------

    do index = 1, size(surface_flux_forcings)
       associate (forcing_field => surface_flux_forcings(index), &
                  metadata      => surface_flux_forcings(index)%metadata)
          select case (trim(metadata%field_source))

          !------------------------------------
          case ('POP monthly calendar')
          !------------------------------------

             file => metadata%field_monthly_calendar_info%forcing_calendar_name

             if (file%has_data) then
                if (thour00 >= file%data_update) then
                   tracer_data_names(1) = file%input%file_varname
                   tracer_bndy_loc(1)   = field_loc_center
                   tracer_bndy_type(1)  = field_type_scalar
                   call update_forcing_data(                            &
                        forcing_time         = file%data_time,          &
                        forcing_time_min_loc = file%data_time_min_loc,  &
                        forcing_interp_type  = file%interp_type,        &
                        forcing_data_next    = file%data_next,          &
                        forcing_data_update  = file%data_update,        &
                        forcing_data_type    = file%data_type,          &
                        forcing_data_inc     = file%data_inc,           &
                        field                = file%data(:,:,:,:,1:12), &
                        forcing_data_rescale = file%data_renorm,        &
                        forcing_data_label   = metadata%marbl_varname,  &
                        forcing_data_names   = tracer_data_names,       &
                        forcing_bndy_loc     = tracer_bndy_loc,         &
                        forcing_bndy_type    = tracer_bndy_type,        &
                        forcing_infile       = file%filename,           &
                        forcing_infile_fmt   = file%input%file_fmt)
                endif

                call interpolate_forcing(                            &
                     interp               = interp_work,             &
                     field                = file%data(:,:,:,:,1:12), &
                     forcing_time         = file%data_time,          &
                     forcing_interp_type  = file%interp_type,        &
                     forcing_time_min_loc = file%data_time_min_loc,  &
                     forcing_interp_freq  = file%interp_freq,        &
                     forcing_interp_inc   = file%interp_inc,         &
                     forcing_interp_next  = file%interp_next,        &
                     forcing_interp_last  = file%interp_last,        &
                     nsteps_run_check     = 0)

                forcing_field%field_0d = interp_work(:,:,:,1)
             endif

          !------------------------------------
          case ('named_field')
          !------------------------------------

          do iblock = 1,nblocks_clinic
             call named_field_get(metadata%field_named_info%field_ind, iblock, &
                                  forcing_field%field_0d(:,:,iblock))
          end do
          !------------------------------------
          case ('internal')
          !------------------------------------

             do iblock = 1,nblocks_clinic

                if (index == ifrac_ind) then
                   forcing_field%field_0d(:,:,iblock) = ifrac(:,:,iblock)

                else if (index == ap_ind) then
                   !  assume PRESS is in cgs units (dyne/cm**2) since that is what is
                   !    required for pressure forcing in barotropic
                   !  want units to be atmospheres
                   !  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
                   !  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
                   forcing_field%field_0d(:,:,iblock) = press(:,:,iblock) / 101.325e+4_r8

                else if (index == sst_ind) then
                   forcing_field%field_0d(:,:,iblock) = sst(:,:,iblock)

                else if (index == sss_ind) then
                   forcing_field%field_0d(:,:,iblock) = sss(:,:,iblock)

                else if (index == box_atm_co2_ind .or. index == box_atm_co2_dup_ind) then
                   forcing_field%field_0d(:,:,iblock) = &
                        ecosys_forcing_saved_state_get_var_val(box_atm_co2_forcing_saved_state_id)

                else if (index == ext_C_flux_ind) then
                   forcing_field%field_0d(:,:,iblock) = c0

                   riv_flux_ind = dic_riv_flux_ind
                   if (riv_flux_ind > 0) then
                     forcing_field%field_0d(:,:,iblock) = forcing_field%field_0d(:,:,iblock) + &
                          riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
                   end if

                   riv_flux_ind = doc_riv_flux_ind
                   if (riv_flux_ind > 0) then
                     forcing_field%field_0d(:,:,iblock) = forcing_field%field_0d(:,:,iblock) + &
                          riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
                   end if

                else if (index == ext_P_flux_ind) then
                   forcing_field%field_0d(:,:,iblock) = c0

                   riv_flux_ind = dip_riv_flux_ind
                   if (riv_flux_ind > 0) then
                     forcing_field%field_0d(:,:,iblock) = forcing_field%field_0d(:,:,iblock) + &
                          riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
                   end if

                   riv_flux_ind = dop_riv_flux_ind
                   if (riv_flux_ind > 0) then
                     forcing_field%field_0d(:,:,iblock) = forcing_field%field_0d(:,:,iblock) + &
                          riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
                   end if

                else if (index == ext_Si_flux_ind) then
                   riv_flux_ind = dsi_riv_flux_ind
                   if (riv_flux_ind > 0) then
                     forcing_field%field_0d(:,:,iblock) = riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
                   else
                     forcing_field%field_0d(:,:,iblock) = c0
                   end if

                else if (index == bc_dep_ind) then
                   ! compute iron_flux in g/cm^2/s

                   ! compute component from atm

                   dust_ratio_to_fe_bioavail_frac = c1 / dust_ratio_to_fe_bioavail_frac_r

                   where (atm_coarse_dust_flux(:,:,iblock) < dust_ratio_thres * atm_fine_dust_flux(:,:,iblock))
                     atm_fe_bioavail_frac(:,:) = fe_bioavail_frac_offset + dust_ratio_to_fe_bioavail_frac * &
                       (dust_ratio_thres - atm_coarse_dust_flux(:,:,iblock) / atm_fine_dust_flux(:,:,iblock))
                   elsewhere
                     atm_fe_bioavail_frac(:,:) = fe_bioavail_frac_offset
                   end where

                   forcing_field%field_0d(:,:,iblock) = atm_fe_bioavail_frac(:,:) * &
                        (iron_frac_in_atm_fine_dust * atm_fine_dust_flux(:,:,iblock) + &
                         iron_frac_in_atm_coarse_dust * atm_coarse_dust_flux(:,:,iblock) + &
                         iron_frac_in_atm_bc * atm_black_carbon_flux(:,:,iblock))

                   ! add component from seaice

                   seaice_fe_bioavail_frac(:,:) = atm_fe_bioavail_frac(:,:)

                   forcing_field%field_0d(:,:,iblock) = forcing_field%field_0d(:,:,iblock) + seaice_fe_bioavail_frac(:,:) * &
                        (iron_frac_in_seaice_dust * seaice_dust_flux(:,:,iblock) + &
                         iron_frac_in_seaice_bc * seaice_black_carbon_flux(:,:,iblock))

                   ! convert to nmol/cm^2/s
                   forcing_field%field_0d(:,:,iblock) = (1.0e9_r8 / molw_Fe) * forcing_field%field_0d(:,:,iblock)

                else if (index == u10sqr_ind) then
                   forcing_field%field_0d(:,:,iblock) = u10_sqr(:,:,iblock)

                else if (index == dust_dep_ind) then
                   forcing_field%field_0d(:,:,iblock) = atm_fine_dust_flux(:,:,iblock) + atm_coarse_dust_flux(:,:,iblock) + &
                        seaice_dust_flux(:,:,iblock)

                else if (index == d13c_ind) then
                   forcing_field%field_0d(:,:,iblock) = d13c(:,:,iblock)

                else if (index == d14c_ind) then
                   forcing_field%field_0d(:,:,iblock) = d14c(:,:,iblock)

                end if  ! index

             end do

          !------------------------------------
          case ('shr_stream')
          !------------------------------------

             stream_index = metadata%field_file_info%strdata_inputlist_ind
             var_ind      = metadata%field_file_info%strdata_var_ind

             n = 0
             do iblock = 1, nblocks_clinic
                this_block = get_block(blocks_clinic(iblock), iblock)
                do j = this_block%jb, this_block%je
                   do i = this_block%ib, this_block%ie
                      n = n + 1
                      shr_stream(i,j,iblock) = surface_strdata_inputlist_ptr(stream_index)%sdat%avs(1)%rAttr(var_ind,n)
                   enddo
                enddo
             enddo

             call POP_HaloUpdate(shr_stream, POP_haloClinic, &
                  POP_gridHorzLocCenter, POP_fieldKindScalar, errorCode, fillValue = 0.0_r8)
             if (errorCode /= POP_Success) then
                call document(subname, 'error updating halo for shr_stream field')
                call exit_POP(sigAbort, 'Stopping in ' // subname)
             endif

             do iblock = 1, nblocks_clinic
                where (land_mask(:,:,iblock))
                   forcing_field%field_0d(:,:,iblock) = shr_stream(:,:,iblock)
                endwhere
             enddo

          end select ! file, constant, driver, shr_stream

          if (metadata%ltime_varying) then
             do iblock = 1, nblocks_clinic
                call apply_unit_conv_factor(land_mask(:,:,iblock), forcing_field, iblock)
             enddo
          end if
       end associate
    end do  ! index

    call adjust_surface_time_varying_data()

    call timer_stop(ecosys_pre_sflux_timer)

  end subroutine ecosys_forcing_set_surface_time_varying_forcing_data

  !***********************************************************************

  subroutine ecosys_forcing_comp_stf_riv(ciso_on, stf_riv, iblock)

    ! DESCRIPTION:
    ! compute river fluxes

    use marbl_constants_mod,                only : R13C_std, R14C_std
    use constants,                          only : p001

    implicit none

    logical,                intent(in)    :: ciso_on
    real (r8),              intent(inout) :: stf_riv(:,:,:)
    integer(kind=int_kind), intent(in)    :: iblock

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:ecosys_forcing_comp_stf_riv'

    integer (int_kind) :: riv_flux_ind
    integer (int_kind) :: processed_field_cnt
    real (r8)          :: conv_factor

    real (r8), parameter :: DOCriv_refract = 0.2_r8
    real (r8), parameter :: DONriv_refract = 0.1_r8
    real (r8), parameter :: DOPriv_refract = 0.025_r8

    !-------------------------------------------------------------------------
    ! River input of isotopic DIC and DOC.
    ! River input of carbon is currently constant and from file.
    ! So the isotopic carbon input is also done very simplified with one value
    ! globally, even though data shows it should vary from river to river.
    !
    ! Using constant delta values of
    ! d13C=-10 permil for DIC (Mook 1986, Raymond et al 2004)
    ! D14C= atmos_D14C - 50 permil for DIC (based on very few data points and
    !       discussion with N. Gruber)
    ! d13C=-27.6 permil for DOC (Raymond et al 2004)
    ! D14C=-50 permil for DOC (Raymond et al 2004), Gruber et al
    !-------------------------------------------------------------------------

    processed_field_cnt = 0

    ! because alk stf_riv has two terms (alk_riv_flux_ind and din_riv_flux_ind),
    ! initialize it to zero and have both terms be cumulative
    stf_riv(:,:,alk_ind) = c0
    stf_riv(:,:,alk_alt_co2_ind) = c0

    riv_flux_ind = din_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,no3_ind)  = riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)

      ! subtract din term from ALK stf_riv
      stf_riv(:,:,alk_ind)  = stf_riv(:,:,alk_ind) - riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
      stf_riv(:,:,alk_alt_co2_ind)  = stf_riv(:,:,alk_alt_co2_ind) - riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
    endif

    riv_flux_ind = dip_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,po4_ind)  = riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
    endif

    riv_flux_ind = don_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,don_ind)  = (c1 - DONriv_refract) * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
      stf_riv(:,:,donr_ind) =       DONriv_refract  * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
    endif

    riv_flux_ind = dop_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,dop_ind)  = (c1 - DOPriv_refract) * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
      stf_riv(:,:,dopr_ind) =       DOPriv_refract  * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
    endif

    riv_flux_ind = dsi_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,sio3_ind) = riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
    endif

    riv_flux_ind = dfe_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,fe_ind)   = riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
    endif

    riv_flux_ind = dic_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,dic_ind)  = riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
      stf_riv(:,:,dic_alt_co2_ind) = riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)

      if (ciso_on) then
        conv_factor = (-10.0_r8 * p001 + c1) * R13C_std
        stf_riv(:,:,di13c_ind) = conv_factor * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)

        conv_factor = ((d14c_glo_avg - 50.0_r8) * p001 + c1) * R14C_std
        stf_riv(:,:,di14c_ind) = conv_factor * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
      endif
    endif

    riv_flux_ind = alk_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,alk_ind)  = stf_riv(:,:,alk_ind) + riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
      stf_riv(:,:,alk_alt_co2_ind)  = stf_riv(:,:,alk_alt_co2_ind) + riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
    endif

    riv_flux_ind = doc_riv_flux_ind
    if (riv_flux_ind > 0) then
      processed_field_cnt   = processed_field_cnt + 1
      stf_riv(:,:,doc_ind)  = (c1 - DOCriv_refract) * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
      stf_riv(:,:,docr_ind) =       DOCriv_refract  * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)

      if (ciso_on) then
        conv_factor = (-27.6_r8 * p001 + c1) * R13C_std
        stf_riv(:,:,do13ctot_ind) = conv_factor * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)

        conv_factor = (-50.0_r8 * p001 + c1) * R14C_std
        stf_riv(:,:,do14ctot_ind) = conv_factor * riv_flux_forcing_fields(riv_flux_ind)%field_0d(:,:,iblock)
      endif
    endif

    if (processed_field_cnt /= size(riv_flux_forcing_fields(:))) then
      call document(subname, 'processed_field_cnt', processed_field_cnt)
      call document(subname, 'size(riv_flux_forcing_fields(:))', size(riv_flux_forcing_fields(:)))
      call exit_POP(sigAbort, 'mismatch between processed_field_cnt and size(riv_flux_forcing_fields(:))')
    end if

  end subroutine ecosys_forcing_comp_stf_riv

  !***********************************************************************

  subroutine adjust_surface_time_varying_data()

    use time_management,    only : imonth
    use marbl_settings_mod, only : parm_Fe_bioavail

    integer :: index, iblock

    !-----------------------------------------------------------------------
    ! Some surface flux forcing fields need to be modified before being
    ! sent to MARBL
    !-----------------------------------------------------------------------

    do iblock = 1,nblocks_clinic

       ! Reduce surface dust flux due to assumed instant surface dissolution
       index = dust_dep_ind
       if (index.gt.0) then
         surface_flux_forcings(index)%field_0d(:,:,iblock) = &
              surface_flux_forcings(index)%field_0d(:,:,iblock) * 0.98_r8
       end if

       ! Add iron patch (if available)
       ! Apply bioavail scaling
       index = Fe_dep_ind
       if (index.gt.0) then
         if (liron_patch .and. imonth == iron_patch_month) then
           surface_flux_forcings(index)%field_0d(:,:,iblock) =      &
                surface_flux_forcings(index)%field_0d(:,:,iblock) + &
                iron_patch_flux(:,:,iblock)
         endif
         surface_flux_forcings(index)%field_0d(:,:,iblock) = &
              surface_flux_forcings(index)%field_0d(:,:,iblock) * parm_Fe_bioavail
       endif

       ! Add iron patch (if available)
       index = bc_dep_ind
       if (index.gt.0) then
         if (liron_patch .and. imonth == iron_patch_month) then
           surface_flux_forcings(index)%field_0d(:,:,iblock) =      &
                surface_flux_forcings(index)%field_0d(:,:,iblock) + &
                iron_patch_flux(:,:,iblock) * parm_Fe_bioavail
         endif
       endif

       index = ifrac_ind
       if (index.gt.0) then
         if (surface_flux_forcings(index)%metadata%field_source == 'internal') then
           where (surface_flux_forcings(index)%field_0d(:,:,iblock) < c0)    &
                surface_flux_forcings(index)%field_0d(:,:,iblock) = c0
           where (surface_flux_forcings(index)%field_0d(:,:,iblock) > c1)    &
                surface_flux_forcings(index)%field_0d(:,:,iblock) = c1
         else
           ! Apply OCMIP ice fraction mask when input is from a file.
           where (surface_flux_forcings(index)%field_0d(:,:,iblock) < 0.2000_r8) &
                surface_flux_forcings(index)%field_0d(:,:,iblock) = 0.2000_r8
           where (surface_flux_forcings(index)%field_0d(:,:,iblock) > 0.9999_r8) &
                surface_flux_forcings(index)%field_0d(:,:,iblock) = 0.9999_r8
         end if
       end if

    end do

    index = dust_dep_ind
    if (index.gt.0) then
      dust_flux_in(:,:,1:nblocks_clinic) = surface_flux_forcings(index)%field_0d(:,:,1:nblocks_clinic)
    end if

  end subroutine adjust_surface_time_varying_data

  !***********************************************************************

  subroutine ciso_update_atm_d13C_D14C (land_mask, d13C, D14C)

    ! Updates module variables d13C and D14C (for atmospheric ratios)

    use grid,                   only : TAREA
    use domain,                 only : blocks_clinic
    use domain,                 only : distrb_clinic
    use blocks,                 only : get_block
    use global_reductions,      only : global_sum
    use forcing_timeseries_mod, only : forcing_timeseries_dataset_update_data
    use forcing_timeseries_mod, only : forcing_timeseries_dataset_get_var
    use c14_atm_forcing_mod,    only : c14_atm_forcing_update_data, c14_atm_forcing_get_data

    implicit none

    logical,   intent(in)  :: land_mask(nx_block, ny_block, max_blocks_clinic)
    real (r8), intent(out) :: d13C(nx_block, ny_block, max_blocks_clinic)  ! atm 13co2 value
    real (r8), intent(out) :: D14C(nx_block, ny_block, max_blocks_clinic)  ! atm 14co2 value

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:ciso_update_atm_d13C_D14C'
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
         atm_d13C_curr, & ! current value of atm d13C (for file option)
         d14c_sum_tmp,  & ! temp for local sum of D14C
         tarea_sum_tmp    ! temp for local sum of TAREA

    !-----------------------------------------------------------------------

    work1(:,:) = c0
    d14c_local_sums(:)  = c0
    tarea_local_sums(:) = c0

    if (trim(ciso_atm_d13c_opt) == 'file') &
      call forcing_timeseries_dataset_update_data(ciso_atm_d13c_forcing_dataset)

    call c14_atm_forcing_update_data

    !-----------------------------------------------------------------------
    ! Loop over blocks
    !-----------------------------------------------------------------------

    do iblock = 1, nblocks_clinic

       !-----------------------------------------------------------------------
       !  Set d13C (constant or from files read in _init) and put on global grid
       !-----------------------------------------------------------------------

       select case (trim(ciso_atm_d13c_opt))
       case ('const')
          d13C(:,:,:) = ciso_atm_d13c_const
       case ('file')
          call forcing_timeseries_dataset_get_var(ciso_atm_d13c_forcing_dataset, &
            varind=1, data_1d=atm_d13C_curr)
          d13C(:,:,iblock) = atm_d13C_curr
       case default
          call document(subname, 'unknown ciso_atm_d13c_opt')
          call exit_POP(sigAbort, 'Stopping in ' // subname)
       end select

       !-----------------------------------------------------------------------
       !  Set D14C
       !-----------------------------------------------------------------------

       call c14_atm_forcing_get_data(iblock, D14C(:,:,iblock))

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

  end subroutine ciso_update_atm_d13C_D14C

  !***********************************************************************

  subroutine ciso_init_atm_D13_D14

    !---------------------------------------------------------------------
    ! !DESCRIPTION:
    !  Initialize surface flux computations for the ecosys_ciso tracer module.
    !  Includes reading CO2 and d13C and D14C data from file if option file is used
    !---------------------------------------------------------------------

    use forcing_timeseries_mod, only : forcing_timeseries_init_dataset
    use c14_atm_forcing_mod,    only : c14_atm_forcing_init

    character(len=*), parameter :: subname = 'ecosys_forcing_mod:ciso_init_atm_D13_D14'

    !-------------------------------------------------------------------------
    !     Set d13C data source
    !-------------------------------------------------------------------------

    call document(subname, 'ciso_atm_d13c_opt', ciso_atm_d13c_opt)

    select case (trim(ciso_atm_d13c_opt))

    case ('const')

       call document(subname, 'ciso_atm_d13c_const', ciso_atm_d13c_const)

    case ('file')

       call forcing_timeseries_init_dataset(ciso_atm_d13c_filename, &
          varnames      = (/ 'delta13co2_in_air' /), &
          model_year    = ciso_atm_model_year, &
          data_year     = ciso_atm_data_year, &
          taxmode_start = 'endpoint', &
          taxmode_end   = 'extrapolate', &
          dataset       = ciso_atm_d13c_forcing_dataset)

    case default

       call exit_POP(sigAbort, 'unknown ciso_atm_d13c_opt')

    end select

    !-------------------------------------------------------------------------
    !     Set D14C data source
    !-------------------------------------------------------------------------

    call c14_atm_forcing_init('ecosys_forcing', ciso_atm_d14c_opt, &
        ciso_atm_d14c_const, ciso_atm_d14c_lat_band_vals, &
        ciso_atm_d14c_filename, ciso_atm_model_year, ciso_atm_data_year)

  end subroutine ciso_init_atm_D13_D14

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

  subroutine forcing_file_init(this, filename, file_varname, rank, &
       year_first, year_last, year_align, tintalgo, taxMode, strdata_inputlist_ptr)

    use strdata_interface_mod, only : POP_strdata_type_set
    use strdata_interface_mod, only : POP_strdata_type_match
    use strdata_interface_mod, only : POP_strdata_type_append_field
    use strdata_interface_mod, only : POP_strdata_type_cp
    use strdata_interface_mod, only : POP_strdata_type_field_count

    class(forcing_file_type),                    intent(inout) :: this
    character(len=*),                            intent(in)    :: filename
    character(len=*),                            intent(in)    :: file_varname
    integer (kind=int_kind),                     intent(in)    :: rank
    integer(kind=int_kind),            optional, intent(in)    :: year_first
    integer(kind=int_kind),            optional, intent(in)    :: year_last
    integer(kind=int_kind),            optional, intent(in)    :: year_align
    character(len=*),                  optional, intent(in)    :: tintalgo
    character(len=*),                  optional, intent(in)    :: taxMode
    type(strdata_input_type), pointer, optional, intent(inout) :: strdata_inputlist_ptr(:)

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    integer (int_kind) :: strdata_inputlist_size
    integer (int_kind) :: m
    type (strdata_input_type)          :: strdata_input_var
    type (strdata_input_type), pointer :: strdata_inputlist_tmp_ptr(:)

    this%filename     = filename
    this%file_varname = file_varname
    if (present(year_first)) this%year_first   = year_first
    if (present(year_last )) this%year_last    = year_last
    if (present(year_align)) this%year_align   = year_align

    if (present(strdata_inputlist_ptr)) then
      call POP_strdata_type_set(strdata_input_var, &
        file_name   = this%filename,     &
        field       = this%file_varname, &
        timer_label = 'marbl_file',      &
        year_first  = this%year_first,   &
        year_last   = this%year_last,    &
        year_align  = this%year_align,   &
        depth_flag  = (rank == 3),       &
        tintalgo    = tintalgo,          &
        taxMode     = taxMode)

      ! if specified file and year specs have already been specified in an strdata_inputlist_ptr entry
      !   then append file_varname to that strdata_inputlist_ptr entry
      ! otherwise generate a new strdata_inputlist_ptr entry

      associate(n => this%strdata_inputlist_ind)
        strdata_inputlist_size = size(strdata_inputlist_ptr)
        do n = 1, strdata_inputlist_size
          if (POP_strdata_type_match(strdata_input_var, strdata_inputlist_ptr(n))) then
            call POP_strdata_type_append_field(this%file_varname, strdata_inputlist_ptr(n))
            exit
          endif
        end do

        if (n > strdata_inputlist_size) then
          allocate(strdata_inputlist_tmp_ptr(n))
          do m = 1, strdata_inputlist_size
            call POP_strdata_type_cp(strdata_inputlist_ptr(m), strdata_inputlist_tmp_ptr(m))
          end do

          deallocate(strdata_inputlist_ptr)
          strdata_inputlist_ptr => strdata_inputlist_tmp_ptr

          call POP_strdata_type_cp(strdata_input_var, strdata_inputlist_ptr(n))
        endif

        this%strdata_var_ind = POP_strdata_type_field_count(strdata_inputlist_ptr(n))
      end associate
    end if

  end subroutine forcing_file_init

  !*****************************************************************************

  subroutine forcing_monthly_calendar_init(this, forcing_calendar_name)

    class(forcing_monthly_calendar_type)   , intent(inout) :: this
    type (forcing_monthly_every_ts), target, intent(in)    :: forcing_calendar_name

    this%forcing_calendar_name => forcing_calendar_name

  end subroutine forcing_monthly_calendar_init

  !*****************************************************************************

  subroutine forcing_field_metadata_set(this, &
       field_source,                          &
       marbl_varname,                         &
       field_units,                           &
       rank,                                  &
       unit_conv_factor,                      &
       field_constant,                        &
       driver_varname,                        &
       named_field,                           &
       filename,                              &
       file_varname,                          &
       year_first, year_last, year_align,     &
       tintalgo, taxMode,                     &
       strdata_inputlist_ptr,                 &
       forcing_calendar_name)

    use io_tools        , only : document

    class(forcing_fields_metadata_type), intent(inout) :: this
    character (len=*),                   intent(in)    :: field_source  ! must  have valid_field_source value)
    character (len=*),                   intent(in)    :: marbl_varname ! required
    character (len=*),                   intent(in)    :: field_units
    integer (kind=int_kind),             intent(in)    :: rank
    real(kind=r8),           optional,   intent(in)    :: unit_conv_factor
    real(kind=r8),           optional,   intent(in)    :: field_constant
    character (len=*),       optional,   intent(in)    :: driver_varname
    character (len=*),       optional,   intent(in)    :: named_field
    character (len=*),       optional,   intent(in)    :: filename
    character (len=*),       optional,   intent(in)    :: file_varname
    integer (kind=int_kind), optional,   intent(in)    :: year_first
    integer (kind=int_kind), optional,   intent(in)    :: year_last
    integer (kind=int_kind), optional,   intent(in)    :: year_align
    character(len=*),        optional,   intent(in)    :: tintalgo
    character(len=*),        optional,   intent(in)    :: taxMode
    type(strdata_input_type), optional, intent(inout), pointer   :: strdata_inputlist_ptr(:)
    type(forcing_monthly_every_ts), optional, target, intent(in) :: forcing_calendar_name

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:forcing_field_metadata_set'

    logical(log_kind)       :: has_valid_inputs
    !-----------------------------------------------------------------------

    call document(subname, "marbl_varname", marbl_varname)
    call document(subname, "field_source", field_source)

    ! required variables for all forcing field sources
    this%field_source  = trim(field_source)
    this%marbl_varname = marbl_varname
    this%field_units   = field_units

    ! optional variables
    this%unit_conv_factor = c1
    if (present(unit_conv_factor)) this%unit_conv_factor = unit_conv_factor

    ! optional variables for forcing field type

    ! each forcing type has its own requirements - if we check here, then the
    ! separate type inits can have fewer optional arguments

    has_valid_inputs = .true.

    select case (trim(field_source))

    case('const')
       this%ltime_varying = .false.
       if (.not.present(field_constant)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          call this%field_constant_info%initialize(field_constant)
       endif

    case('zero')
       this%ltime_varying = .false.
       call this%field_constant_info%initialize(c0)

    case('internal')
       this%ltime_varying = .true.
       if (.not.present(driver_varname)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          call this%field_driver_info%initialize(driver_varname)
       endif

    case('named_field')
       this%ltime_varying = .true.
       if (.not.present(named_field)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          call this%field_named_info%initialize(named_field)
       endif

    case('file_time_invariant')
       this%ltime_varying = .false.
       if (.not.present(filename))     has_valid_inputs = .false.
       if (.not.present(file_varname)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          call this%field_file_info%initialize(&
               filename, file_varname, rank)
       endif

    case('shr_stream')
       this%ltime_varying = .true.
       if (.not.present(filename))              has_valid_inputs = .false.
       if (.not.present(file_varname))          has_valid_inputs = .false.
       if (.not.present(year_first))            has_valid_inputs = .false.
       if (.not.present(year_last))             has_valid_inputs = .false.
       if (.not.present(year_align))            has_valid_inputs = .false.
       if (.not.present(strdata_inputlist_ptr)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          call this%field_file_info%initialize(&
               filename, file_varname, rank, &
               year_first=year_first, year_last=year_last, year_align=year_align, &
               tintalgo=tintalgo, taxMode=taxMode, &
               strdata_inputlist_ptr=strdata_inputlist_ptr)
       endif

    case('POP monthly calendar')
       this%ltime_varying = .true.
       if (.not.present(forcing_calendar_name)) has_valid_inputs = .false.
       if (has_valid_inputs) then
          call this%field_monthly_calendar_info%initialize(forcing_calendar_name)
       endif

     case default
       call document(subname, "unknown field_source")
       call exit_POP(sigAbort, 'Stopping in ' // subname)

    end select

    if (.not.has_valid_inputs) then
       call document(subname, "required optional arguments for field_source not provided")
       call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

   end subroutine forcing_field_metadata_set

  !*****************************************************************************

  subroutine forcing_fields_add(this,       &
       field_source,                        &
       marbl_varname,                       &
       field_units,                         &
       id,                                  &
       rank,                                &
       unit_conv_factor,                    &
       dim3_len,                            &
       ldim3_is_depth,                      &
       field_constant,                      &
       driver_varname,                      &
       named_field,                         &
       filename, file_varname,              &
       year_first, year_last, year_align,   &
       tintalgo, taxMode,                   &
       strdata_inputlist_ptr,               &
       forcing_calendar_name)

    class(forcing_fields_type)      , intent(inout) :: this
    character(len=*)                , intent(in)    :: field_source
    character(len=*)                , intent(in)    :: marbl_varname
    character(len=*)                , intent(in)    :: field_units
    integer(kind=int_kind)          , intent(in)    :: id
    integer(kind=int_kind)          , intent(in)    :: rank
    real(kind=r8),          optional, intent(in)    :: unit_conv_factor
    integer(kind=int_kind), optional, intent(in)    :: dim3_len
    logical(kind=log_kind), optional, intent(in)    :: ldim3_is_depth
    real(kind=r8),          optional, intent(in)    :: field_constant
    character(len=*),       optional, intent(in)    :: driver_varname
    character(len=*),       optional, intent(in)    :: named_field
    character(len=*),       optional, intent(in)    :: filename
    character(len=*),       optional, intent(in)    :: file_varname
    integer(kind=int_kind), optional, intent(in)    :: year_first
    integer(kind=int_kind), optional, intent(in)    :: year_last
    integer(kind=int_kind), optional, intent(in)    :: year_align
    character(len=*),       optional, intent(in)    :: tintalgo
    character(len=*),       optional, intent(in)    :: taxMode
    type(strdata_input_type), optional, intent(inout), pointer   :: strdata_inputlist_ptr(:)
    type(forcing_monthly_every_ts), optional, target, intent(in) :: forcing_calendar_name

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'ecosys_forcing_mod:forcing_fields_add'
    !-----------------------------------------------------------------------

    this%rank = rank

    ! confirm that rank == 2 or 3
    if (rank /= 2 .and. rank /= 3) then
      call document(subname, 'marbl_varname', marbl_varname)
      call document(subname, 'rank must be 2 or 3, rank', rank)
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    ! confirm that dim3_len is present if rank == 3
    if (rank == 3 .and. .not. present(dim3_len)) then
      call document(subname, 'marbl_varname', marbl_varname)
      call document(subname, 'rank=3, but dim3_len argument is not present')
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    this%ldim3_is_depth = (rank == 3)
    if (present(ldim3_is_depth)) this%ldim3_is_depth = ldim3_is_depth

    ! confirm that rank == 3 if ldim3_is_depth is .true.
    if (this%ldim3_is_depth .and. rank == 2) then
      call document(subname, 'marbl_varname', marbl_varname)
      call document(subname, 'ldim3_is_depth is true, but rank == 2')
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    ! confirm that dim3_len == km if ldim3_is_depth is .true.
    if (this%ldim3_is_depth) then
      if (dim3_len /= km) then
        call document(subname, 'marbl_varname', marbl_varname)
        call document(subname, 'dim3_len', dim3_len)
        call document(subname, 'ldim3_is_depth is true, but dim3 size /= km')
        call exit_POP(sigAbort, 'Stopping in ' // subname)
      end if
    end if

    ! confirm that source == 'internal' if rank is 3 and ldim3_is_depth is .false.
    if (rank == 3 .and. .not. this%ldim3_is_depth .and. trim(field_source) /= 'internal') then
      call document(subname, 'marbl_varname', marbl_varname)
      call document(subname, 'field_source', field_source)
      call document(subname, 'rank=3 and ldim3_is_depth is false, but field_source /= internal')
      call exit_POP(sigAbort, 'Stopping in ' // subname)
    end if

    call this%metadata%set(                              &
         field_source, marbl_varname,                    &
         field_units, rank,                              &
         unit_conv_factor=unit_conv_factor,              &
         field_constant=field_constant,                  &
         driver_varname=driver_varname,                  &
         named_field=named_field,                        &
         filename=filename, file_varname=file_varname,   &
         year_first=year_first,                          &
         year_last=year_last,                            &
         year_align=year_align,                          &
         tintalgo=tintalgo, taxMode=taxMode,             &
         strdata_inputlist_ptr=strdata_inputlist_ptr,    &
         forcing_calendar_name=forcing_calendar_name)

    if (rank == 2) then
      allocate(this%field_0d(nx_block, ny_block, nblocks_clinic))
      this%field_0d = c0
    end if

    if (rank == 3) then
      allocate(this%field_1d(nx_block, ny_block, dim3_len, nblocks_clinic))
      this%field_1d = c0
    end if

  end subroutine forcing_fields_add

  !*****************************************************************************

  subroutine init_monthly_surface_flux_forcing_metadata(var)

    implicit none

    type(forcing_monthly_every_ts), intent(inout) :: var

    var%interp_type = 'linear'
    var%data_type   = 'monthly-calendar'
    var%interp_freq = 'every-timestep'
    var%filename    = 'not-used-for-monthly'
    var%data_label  = 'not-used-for-monthly'

  end subroutine init_monthly_surface_flux_forcing_metadata

  !*****************************************************************************

  function ecosys_forcing_tracer_ref_val(ind)
    !
    ! !DESCRIPTION:
    !  return reference value for tracer with global tracer index ind
    !  this is used in virtual flux computations

    implicit none

    integer(int_kind) , intent(in) :: ind
    real(r8) :: ecosys_forcing_tracer_ref_val

    ecosys_forcing_tracer_ref_val = surf_avg(ind)

  end function ecosys_forcing_tracer_ref_val

  !*****************************************************************************

end module ecosys_forcing_mod
