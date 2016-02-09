module marbl_interface_types
  ! module for definitions of types that are shared between marbl interior and the driver.

  use marbl_kinds_mod, only : c0, c1, r8, log_kind, int_kind, char_len
  use marbl_interface_constants, only : marbl_str_length

  implicit none

  private

  !*****************************************************************************

  type, public :: marbl_domain_type
     integer(int_kind)     :: km                            ! number of vertical grid cells
     integer(int_kind)     :: kmt                           ! index of ocean floor
     integer(int_kind)     :: num_PAR_subcols               ! number of PAR subcols
     integer(int_kind)     :: num_elements_surface_forcing  ! number of surface forcing columns
     integer(int_kind)     :: num_elements_interior_forcing ! number of interior forcing columns
     real(r8), allocatable :: zt(:)                         ! (km) vert dist from sfc to midpoint of layer
     real(r8), allocatable :: zw(:)                         ! (km) vert dist from sfc to bottom of layer
     real(r8), allocatable :: delta_z(:)                    ! (km) delta z - different values for partial bottom cells
     real(r8), allocatable :: dz(:)                         ! (km) delta z - same values for partial bottom cells
     logical(log_kind)     :: land_mask
   contains
     procedure, public :: construct => marbl_domain_constructor
  end type marbl_domain_type

  !*****************************************************************************

  type, public :: marbl_tracer_metadata_type
     character(char_len) :: short_name
     character(char_len) :: long_name
     character(char_len) :: units
     character(char_len) :: tend_units
     character(char_len) :: flux_units
     logical             :: lfull_depth_tavg
     real(r8)            :: scale_factor
  end type marbl_tracer_metadata_type

  !*****************************************************************************

  type, public :: marbl_tracer_read_type
     character(char_len) :: mod_varname
     character(char_len) :: filename
     character(char_len) :: file_varname
     character(char_len) :: file_fmt
     real(r8)            :: scale_factor
     real(r8)            :: default_val
  end type marbl_tracer_read_type

  !*****************************************************************************
  ! FIXME(mnl,2016-01) move PAR_col_frac and surf_shortwave into this datatype
  !                    and come up with better name
  type, public :: marbl_gcm_state_type
     real(r8), allocatable :: temperature(:)    ! (km)
     real(r8), allocatable :: salinity(:)       ! (km)
     real(r8), allocatable :: PAR_col_frac(:)   ! column fraction occupied by each sub-column
     real(r8), allocatable :: surf_shortwave(:) ! surface shortwave for each sub-column (W/m^2)
   contains
     procedure, public :: construct => marbl_gcm_state_constructor
  end type marbl_gcm_state_type

  !*****************************************************************************

  type, private :: marbl_single_diagnostic_type
     ! marbl_singl_diagnostic : 
     ! a private type, this contains both the metadata
     ! and the actual diagnostic data for a single
     ! diagnostic quantity. Data must be accessed via
     ! the marbl_diagnostics_type data structure.
     character(len=char_len) :: long_name
     character(len=char_len) :: short_name
     character(len=char_len) :: units
     character(len=char_len) :: vertical_grid ! 'none', 'layer_avg', 'layer_iface'
     logical(log_kind)       :: compute_now
     logical(log_kind)       :: ltruncated_vertical_extent
     real(r8), allocatable, dimension(:) :: field_2d
     real(r8), allocatable, dimension(:,:) :: field_3d

   contains
     procedure :: initialize  => marbl_single_diag_init
  end type marbl_single_diagnostic_type

  type, public :: marbl_diagnostics_type
     ! marbl_diagnostics : 
     ! used to pass diagnostic information from marbl back to
     ! the driver.
     integer :: diag_cnt
     integer :: num_elements
     integer :: num_levels
     type(marbl_single_diagnostic_type), dimension(:), allocatable :: diags

   contains
     procedure, public :: construct      => marbl_diagnostics_constructor
     procedure, public :: set_to_zero    => marbl_diagnostics_set_to_zero
     procedure, public :: add_diagnostic => marbl_diagnostics_add
     procedure, public :: deconstruct    => marbl_diagnostics_deconstructor
  end type marbl_diagnostics_type

  !*****************************************************************************

  type, public :: marbl_forcing_input_type
     logical (log_kind), allocatable, dimension(:)   :: land_mask
     real (r8) :: seconds_in_year  ! this is set by the gcm and can change in time if leap year
     real (r8) :: d14c_glo_avg     ! this is computed by the gcm            

     real (r8), allocatable, dimension(:)   :: u10_sqr         
     real (r8), allocatable, dimension(:)   :: ifrac           ! ice fraction
     real (r8), allocatable, dimension(:)   :: atm_press
     real (r8), allocatable, dimension(:)   :: sst             
     real (r8), allocatable, dimension(:)   :: sss             
     real (r8), allocatable, dimension(:)   :: xco2            
     real (r8), allocatable, dimension(:)   :: xco2_alt_co2    
     real (r8), allocatable, dimension(:)   :: ap         
     real (r8), allocatable, dimension(:)   :: dust_flux
     real (r8), allocatable, dimension(:)   :: xkw        
     real (r8), allocatable, dimension(:)   :: iron_flux    
     real (r8), allocatable, dimension(:)   :: ph_prev         
     real (r8), allocatable, dimension(:)   :: ph_prev_alt_co2 
     real (r8), allocatable, dimension(:)   :: d13c
     real (r8), allocatable, dimension(:)   :: d14c

     real (r8), allocatable, dimension(:,:) :: input_forcings
     real (r8), allocatable, dimension(:,:) :: surface_vals
   contains
     procedure, public :: construct => marbl_forcing_input_constructor
  end type marbl_forcing_input_type

  !*****************************************************************************

  type, public :: marbl_forcing_output_type
     real (r8), allocatable, dimension(:)   :: ph_prev         
     real (r8), allocatable, dimension(:)   :: ph_prev_alt_co2 
     real (r8), allocatable, dimension(:)   :: iron_flux    
     real (r8), allocatable, dimension(:)   :: flux_o2     
     real (r8), allocatable, dimension(:)   :: flux_co2     
     real (r8), allocatable, dimension(:)   :: flux_alt_co2 ! tracer flux alternative CO2 (nmol/cm^2/s)
     real (r8), allocatable, dimension(:)   :: co2star
     real (r8), allocatable, dimension(:)   :: dco2star
     real (r8), allocatable, dimension(:)   :: pco2surf
     real (r8), allocatable, dimension(:)   :: dpco2
     real (r8), allocatable, dimension(:)   :: co3
     real (r8), allocatable, dimension(:)   :: co2star_alt
     real (r8), allocatable, dimension(:)   :: dco2star_alt
     real (r8), allocatable, dimension(:)   :: pco2surf_alt
     real (r8), allocatable, dimension(:)   :: dpco2_alt
     real (r8), allocatable, dimension(:)   :: schmidt_co2  ! Schmidt number
     real (r8), allocatable, dimension(:)   :: schmidt_o2   ! Schmidt number
     real (r8), allocatable, dimension(:)   :: pv_o2        ! piston velocity (cm/s)
     real (r8), allocatable, dimension(:)   :: pv_co2       ! piston velocity (cm/s)
     real (r8), allocatable, dimension(:)   :: o2sat        ! used O2 saturation (mmol/m^3)
     real (r8), allocatable, dimension(:,:) :: stf_module
     real (r8), allocatable, dimension(:,:) :: stf_ciso
   contains
     procedure, public :: construct => marbl_forcing_output_constructor
  end type marbl_forcing_output_type

  !*****************************************************************************

  type, public :: forcing_monthly_every_ts
     type(marbl_tracer_read_type)            :: input
     logical(log_kind)                       :: has_data
     real(r8), dimension(:,:,:,:,:), pointer :: DATA
     character(char_len)                     :: interp_type        ! = 'linear'
     character(char_len)                     :: data_type          ! = 'monthly-calendar'
     character(char_len)                     :: interp_freq        ! = 'every-timestep'
     character(char_len)                     :: filename           ! = 'not-used-for-monthly'
     character(char_len)                     :: data_label         ! = 'not-used-for-monthly'
     real(r8), dimension(12)                 :: data_time          ! times where DATA is given
     real(r8), dimension(20)                 :: data_renorm        ! not used for monthly
     real(r8)                                :: data_inc           ! not used for monthly data
     real(r8)                                :: data_next          ! time that will be used for the next
                                                                   ! value of forcing data that is needed
     real(r8)                                :: data_update        ! time when the a new forcing value
                                                                   ! needs to be added to interpolation set
     real(r8)                                :: interp_inc         ! not used for 'every-timestep' interp
     real(r8)                                :: interp_next        ! not used for 'every-timestep' interp
     real(r8)                                :: interp_last        ! not used for 'every-timestep' interp
     integer(int_kind)                       :: data_time_min_loc  ! index of the third dimension of data_time
                                                                   ! containing the minimum forcing time
  end type forcing_monthly_every_ts

  !*****************************************************************************

  !-------------------------------------------------
  type, private :: marbl_forcing_constant_type
     real(KIND=r8) :: field_constant           ! constant value for field_source
   contains
      procedure :: initialize  => marbl_forcing_constant_init
  end type marbl_forcing_constant_type
  !-------------------------------------------------

  !-------------------------------------------------
  type, private :: marbl_forcing_driver_type
     character(char_len) :: marbl_driver_varname
   contains
     procedure :: initialize  => marbl_forcing_driver_init
  end type marbl_forcing_driver_type
  !-------------------------------------------------

  !-------------------------------------------------
  type, private :: marbl_forcing_file_type
     character(char_len)    :: filename
     character(char_len)    :: file_varname
     character(char_len)    :: temporal      ! temporarily to support current I/O routines
     integer(KIND=int_kind) :: year_first
     integer(KIND=int_kind) :: year_last
     integer(KIND=int_kind) :: year_align
     integer(KIND=int_kind) :: date
     integer(KIND=int_kind) :: time
   contains
     procedure :: initialize  => marbl_forcing_file_init
  end type marbl_forcing_file_type
  !-------------------------------------------------

  !-------------------------------------------------
  type, private :: marbl_forcing_monthly_calendar_type
     type (forcing_monthly_every_ts), pointer :: marbl_forcing_calendar_name
   contains
     procedure :: initialize  => marbl_forcing_monthly_calendar_init
  end type marbl_forcing_monthly_calendar_type
  !-------------------------------------------------

  !-------------------------------------------------
  ! single_forcing_field_type (contains the above 4 type definitions)
  type, private :: marbl_single_forcing_field_type
     character(char_len)                        :: marbl_varname
     character(char_len)                        :: field_units          ! units represent what is in field_data,
                                                                        ! not the file (up to driver to do unit conversion)
     character(char_len)                        :: field_source         ! "file", "driver", "POP monthly calendar", 
                                                                        ! "constant", "none"
     character(char_len)                        :: temporal_interp      ! information on interpolation scheme used
                                                                        ! to populate field_data
     real(KIND=r8)                              :: unit_conv_factor     ! unit conversion factor, incorporates scale_factor
     logical (log_kind)                         :: has_data             ! TODO need this?
     real(KIND=r8), dimension(:), allocatable   :: field_data           ! only allocate if field_source != 'none'
     type (marbl_forcing_constant_type)         :: field_constant_info
     type (marbl_forcing_driver_type)           :: field_driver_info
     type (marbl_forcing_file_type)             :: field_file_info
     type (marbl_forcing_monthly_calendar_type) :: field_monthly_calendar_info
   contains
     procedure :: initialize  => marbl_single_forcing_field_init
  end type marbl_single_forcing_field_type
  !-------------------------------------------------

  !-------------------------------------------------
  type, public :: marbl_forcing_fields_type
     integer(KIND=int_kind) :: num_elements
     integer(KIND=int_kind) :: forcing_field_cnt
     type(marbl_single_forcing_field_type), dimension(:), allocatable :: forcing_fields
   contains
     procedure, public :: construct         => marbl_forcing_fields_constructor
     procedure, public :: add_forcing_field => marbl_forcing_fields_add
     procedure, public :: deconstruct       => marbl_forcing_fields_deconstructor
  end type marbl_forcing_fields_type
  !-------------------------------------------------

  !*****************************************************************************

contains

  !*****************************************************************************

  subroutine marbl_domain_constructor(this, &
       num_levels, num_PAR_subcols, num_elements_surface_forcing, num_elements_interior_forcing)

    class(marbl_domain_type), intent(inout) :: this
    integer , intent(in)    :: num_levels
    integer , intent(in)    :: num_PAR_subcols
    integer , intent(in)    :: num_elements_surface_forcing
    integer , intent(in)    :: num_elements_interior_forcing

    this%km = num_levels
    this%num_PAR_subcols = num_PAR_subcols
    this%num_elements_surface_forcing = num_elements_surface_forcing
    this%num_elements_interior_forcing = num_elements_interior_forcing

    allocate(this%dz(num_levels))
    allocate(this%delta_z(num_levels))
    allocate(this%zw(num_levels))
    allocate(this%zt(num_levels))

  end subroutine marbl_domain_constructor
  
  !*****************************************************************************

  subroutine marbl_gcm_state_constructor(this, num_levels, num_PAR_subcols)

    class(marbl_gcm_state_type) , intent(inout) :: this
    integer , intent(in)    :: num_levels
    integer , intent(in)    :: num_PAR_subcols

    allocate(this%temperature(num_levels))
    allocate(this%salinity(num_levels))

    allocate(this%PAR_col_frac(num_PAR_subcols))
    allocate(this%surf_shortwave(num_PAR_subcols))

  end subroutine marbl_gcm_state_constructor

  !*****************************************************************************

  subroutine marbl_single_diag_init(this, lname, sname, units, vgrid,         &
             truncate, num_elements, num_levels)

    class(marbl_single_diagnostic_type) , intent(inout) :: this
    character(len=char_len) , intent(in)    :: lname, sname, units, vgrid
    logical                 , intent(in)    :: truncate
    integer                 , intent(in)    :: num_elements
    integer                 , intent(in)    :: num_levels

    ! Allocate column memory for 3D vars or num_elements memory for 2D vars
    select case (trim(vgrid))
      case ('layer_avg')
        allocate(this%field_3d(num_levels, num_elements))
      case ('layer_iface')
        allocate(this%field_3d(num_levels+1, num_elements))
      case ('none')
        allocate(this%field_2d(num_elements))
      case DEFAULT
        print*, "ERROR: ", trim(vgrid), " is not a valid vertical grid for MARBL"
        ! FIXME: abort
    end select

    this%compute_now = .true.
    this%long_name = trim(lname)
    this%short_name = trim(sname)
    this%units = trim(units)
    this%vertical_grid = trim(vgrid)
    this%ltruncated_vertical_extent = truncate

  end subroutine marbl_single_diag_init

  !*****************************************************************************

  subroutine marbl_forcing_input_constructor(this, &
       num_elements, num_surface_vals, num_input_forcings, ciso_on)
    
    class(marbl_forcing_input_type), intent(inout) :: this
    integer (int_kind), intent(in) :: num_elements
    integer (int_kind), intent(in) :: num_surface_vals
    integer (int_kind), intent(in) :: num_input_forcings
    logical(log_kind) , intent(in) :: ciso_on

    allocate(this%u10_sqr         (num_elements))         
    allocate(this%ifrac           (num_elements))           
    allocate(this%land_mask       (num_elements))
    allocate(this%sst             (num_elements))             
    allocate(this%sss             (num_elements))             
    allocate(this%xco2            (num_elements))            
    allocate(this%atm_press       (num_elements))         
    allocate(this%xco2_alt_co2    (num_elements))    
    allocate(this%xkw             (num_elements))        
    allocate(this%dust_flux       (num_elements))    
    allocate(this%iron_flux       (num_elements))
    allocate(this%ph_prev         (num_elements))
    allocate(this%ph_prev_alt_co2 (num_elements)) 
    allocate(this%input_forcings  (num_elements, num_input_forcings))

    if (ciso_on) then
       allocate(this%d13c(num_elements))
       allocate(this%d14c(num_elements))
    end if

    allocate(this%surface_vals(num_elements, num_surface_vals))
  end subroutine marbl_forcing_input_constructor

  !*****************************************************************************

  subroutine marbl_forcing_output_constructor(this, &
       num_elements, num_surface_vals)

    class(marbl_forcing_output_type), intent(inout) :: this
    integer (int_kind), intent(in) :: num_elements
    integer (int_kind), intent(in) :: num_surface_vals

    allocate(this%ph_prev         (num_elements))
    allocate(this%ph_prev_alt_co2 (num_elements)) 
    allocate(this%iron_flux       (num_elements)) 
    allocate(this%flux_co2        (num_elements))            
    allocate(this%flux_o2         (num_elements))            
    allocate(this%co2star         (num_elements))
    allocate(this%dco2star        (num_elements))
    allocate(this%pco2surf        (num_elements))
    allocate(this%dpco2           (num_elements))
    allocate(this%co3             (num_elements))
    allocate(this%co2star_alt     (num_elements))
    allocate(this%dco2star_alt    (num_elements))
    allocate(this%pco2surf_alt    (num_elements))
    allocate(this%dpco2_alt       (num_elements))
    allocate(this%schmidt_co2     (num_elements)) 
    allocate(this%schmidt_o2      (num_elements))  
    allocate(this%pv_o2           (num_elements))       
    allocate(this%pv_co2          (num_elements))      
    allocate(this%o2sat           (num_elements))       
    allocate(this%flux_alt_co2    (num_elements))
    allocate(this%stf_module      (num_elements, num_surface_vals))

  end subroutine marbl_forcing_output_constructor

  !*****************************************************************************

  subroutine marbl_diagnostics_constructor(this, num_diags, num_elements,     &
             num_levels)

    class(marbl_diagnostics_type), intent(inout) :: this
    integer (int_kind),            intent(in)    :: num_diags
    integer (int_kind),            intent(in)    :: num_elements
    integer (int_kind),            intent(in)    :: num_levels

    allocate(this%diags(num_diags))
    this%diag_cnt = 0
    this%num_elements = num_elements
    this%num_levels = num_levels

  end subroutine marbl_diagnostics_constructor

  !*****************************************************************************

  subroutine marbl_diagnostics_set_to_zero(this)

    class(marbl_diagnostics_type), intent(inout) :: this

    integer (int_kind) :: n

    do n=1,this%diag_cnt
      if (allocated(this%diags(n)%field_2d)) then
        this%diags(n)%field_2d(:) = c0
      elseif (allocated(this%diags(n)%field_3d)) then
        this%diags(n)%field_3d(:, :) = c0
      else
        ! TODO abort abort abort
        write(*,*) "ERROR: neither field_2d nor field_3d are allocated"
        write(*,*) "Diag short name = ", trim(this%diags(n)%short_name)
        write(*,*) "Diag long name = ", trim(this%diags(n)%long_name)
      end if
    end do

  end subroutine marbl_diagnostics_set_to_zero

  !*****************************************************************************

  subroutine marbl_diagnostics_add(this, lname, sname, units, vgrid,          &
             truncate, id)

    class(marbl_diagnostics_type) , intent(inout) :: this
    character(len=char_len)       , intent(in)    :: lname, sname, units, vgrid
    logical (int_kind)            , intent(in)    :: truncate
    integer (int_kind)            , intent(out)   :: id

    this%diag_cnt = this%diag_cnt + 1
    id = this%diag_cnt
    if (id.gt.size(this%diags)) then
      print*, "ERROR: increase max number of diagnostics!"
      ! FIXME: abort
    end if
    call this%diags(id)%initialize(lname, sname, units, vgrid, truncate,      &
         this%num_elements, this%num_levels)

  end subroutine marbl_diagnostics_add

  !*****************************************************************************

  subroutine marbl_diagnostics_deconstructor(this)

    class(marbl_diagnostics_type), intent(inout) :: this

    integer :: n

    do n=1,size(this%diags)
      if (allocated(this%diags(n)%field_2d)) then
        deallocate(this%diags(n)%field_2d)
      end if
      if (allocated(this%diags(n)%field_3d)) then
        deallocate(this%diags(n)%field_3d)
      end if
    end do
    deallocate(this%diags)

  end subroutine marbl_diagnostics_deconstructor

  !*****************************************************************************

  subroutine marbl_forcing_constant_init(this, field_constant)

    class(marbl_forcing_constant_type), intent(inout) :: this
    real(KIND=r8),                      intent(in)    :: field_constant

    this%field_constant = field_constant

  end subroutine marbl_forcing_constant_init

  !*****************************************************************************

  subroutine marbl_forcing_driver_init(this, marbl_driver_varname)

    class(marbl_forcing_driver_type), intent(inout) :: this
    character(char_len),              intent(in)    :: marbl_driver_varname

    this%marbl_driver_varname = marbl_driver_varname

  end subroutine marbl_forcing_driver_init

  !*****************************************************************************

  subroutine marbl_forcing_file_init(this, filename, file_varname, temporal, year_first, &
                                     year_last, year_align, date, time)


    class(marbl_forcing_file_type),   intent(inout) :: this
    character(char_len),              intent(in)    :: filename
    character(char_len),              intent(in)    :: file_varname
    character(char_len),    optional, intent(in)    :: temporal
    integer(KIND=int_kind), optional, intent(in)    :: year_first
    integer(KIND=int_kind), optional, intent(in)    :: year_last
    integer(KIND=int_kind), optional, intent(in)    :: year_align
    integer(KIND=int_kind), optional, intent(in)    :: date
    integer(KIND=int_kind), optional, intent(in)    :: time

    this%filename     = filename
    this%file_varname = file_varname
    if (present(temporal  )) this%temporal   = temporal
    if (present(year_first)) this%year_first = year_first
    if (present(year_last )) this%year_last  = year_last
    if (present(year_align)) this%year_align = year_align
    if (present(date      )) this%date       = date
    if (present(time      )) this%time       = time

  end subroutine marbl_forcing_file_init

  !*****************************************************************************

  subroutine marbl_forcing_monthly_calendar_init(this, marbl_forcing_calendar_name)

    class(marbl_forcing_monthly_calendar_type), intent(inout) :: this
    type (forcing_monthly_every_ts),    target, intent(in)    :: marbl_forcing_calendar_name

    this%marbl_forcing_calendar_name => marbl_forcing_calendar_name

  end subroutine marbl_forcing_monthly_calendar_init

  !*****************************************************************************

  subroutine marbl_single_forcing_field_init(this, num_elements, field_source, marbl_varname, &
                                             field_units, unit_conv_factor, temporal_interp,  &
                                             field_constant, marbl_driver_varname, filename,  &
                                             file_varname, temporal, year_first, year_last,   &
                                             year_align, date, time, marbl_forcing_calendar_name)

    class(marbl_single_forcing_field_type), intent(inout) :: this
    integer (KIND=int_kind),                intent(in)    :: num_elements
    character (char_len),                   intent(in)    :: field_source
    character (char_len),                   intent(in)    :: marbl_varname
    character (char_len),                   intent(in)    :: field_units
    real(KIND=r8),           optional,      intent(in)    :: unit_conv_factor
    character (char_len),    optional,      intent(in)    :: temporal_interp
    real(KIND=r8),           optional,      intent(in)    :: field_constant
    character (char_len),    optional,      intent(in)    :: marbl_driver_varname
    character (char_len),    optional,      intent(in)    :: filename
    character (char_len),    optional,      intent(in)    :: file_varname
    character (char_len),    optional,      intent(in)    :: temporal
    integer (KIND=int_kind), optional,      intent(in)    :: year_first
    integer (KIND=int_kind), optional,      intent(in)    :: year_last
    integer (KIND=int_kind), optional,      intent(in)    :: year_align
    integer (KIND=int_kind), optional,      intent(in)    :: date
    integer (KIND=int_kind), optional,      intent(in)    :: time
    type (forcing_monthly_every_ts), optional, target, intent(in) :: marbl_forcing_calendar_name

    character(len=char_len), dimension(6) :: valid_field_sources
    integer (KIND=int_kind) :: n
    logical (log_kind)      :: has_valid_source
    logical (log_kind)      :: has_valid_inputs

    valid_field_sources(1) = "constant"
    valid_field_sources(2) = "driver"
    valid_field_sources(3) = "file"
    valid_field_sources(4) = "marbl"
    valid_field_sources(5) = "POP monthly calendar"
    valid_field_sources(6) = "none"

    ! set defaults
    this%unit_conv_factor = c1
    this%temporal_interp  = ''
    this%has_data         = .false.

    ! check for valid source
    has_valid_source = .false.
    do n = 1,size(valid_field_sources)
      if (trim(field_source).EQ.trim(valid_field_sources(n))) has_valid_source = .true.
    enddo
    if (.NOT.has_valid_source) then
      write(*,*) "ERROR: ", trim(field_source), "is not a valid field source for MARBL"
      ! FIXME: return error code
    endif

    this%field_source = trim(field_source)
    if (trim(field_source) .NE. "none") then
      this%has_data      = .true.
      allocate(this%field_data(num_elements))
      this%field_data(:) = c0
    endif

    ! required variables for all forcing field sources
    this%marbl_varname = marbl_varname
    this%field_units   = field_units

    ! optional variables for forcing field type
    if (present(unit_conv_factor)) this%unit_conv_factor = unit_conv_factor
    if (present(temporal_interp )) this%temporal_interp  = temporal_interp

    ! each forcing type has its own requirements - if we check here, then the
    ! separate type inits can have fewer optional arguments
    has_valid_inputs = .true.
    if (trim(field_source) .EQ. "constant") then
      if (.NOT.present(field_constant)) has_valid_inputs = .false.
      if (has_valid_inputs) then
        write(*,*) "Adding constant forcing_field_type for ", this%marbl_varname 
!JW        call this%field_constant_info%initialize(field_constant)
        call marbl_forcing_constant_init(this%field_constant_info, field_constant)
      else
        write(*,*) "ERROR: Call to MARBL does not have the correct optional arguments for ", trim(field_source)
        ! FIXME: return error code
      endif
    endif

    if (trim(field_source) .EQ. "driver") then
      if (.NOT.present(marbl_driver_varname)) has_valid_inputs = .false.
      if (has_valid_inputs) then
        write(*,*) "Adding driver forcing_field_type for ", this%marbl_varname 
        call this%field_driver_info%initialize(marbl_driver_varname)
      else
        write(*,*) "ERROR: Call to MARBL does not have the correct optional arguments for ", trim(field_source)
        ! FIXME: return error code
      endif
    endif
    if (trim(field_source) .EQ. "file") then
      if (.NOT.present(filename))     has_valid_inputs = .false.
      if (.NOT.present(file_varname)) has_valid_inputs = .false.
      if (has_valid_inputs) then
        write(*,*) "Adding file forcing_field_type for ", this%marbl_varname 
        call this%field_file_info%initialize(filename, file_varname, &
                                             temporal=temporal, year_first=year_first,   &
                                             year_last=year_last, year_align=year_align, &
                                             date=date, time=time)
      else
        write(*,*) "ERROR: Call to MARBL does not have the correct optional arguments for ", trim(field_source)
        ! FIXME: return error code
      endif
    endif
    if (trim(field_source) .EQ. "POP monthly calendar") then
      if (.NOT.present(marbl_forcing_calendar_name)) has_valid_inputs = .false.
      if (has_valid_inputs) then
        write(*,*) "Adding calendar forcing_field_type for ", this%marbl_varname 
        call this%field_monthly_calendar_info%initialize(marbl_forcing_calendar_name)
      else
        write(*,*) "ERROR: Call to MARBL does not have the correct optional arguments for ", trim(field_source)
        ! FIXME: return error code
      endif
    endif

   end subroutine marbl_single_forcing_field_init

  !*****************************************************************************

  subroutine marbl_forcing_fields_constructor(this, num_elements, num_forcing_fields)

    class(marbl_forcing_fields_type), intent(inout) :: this
    integer (int_kind),               intent(in)    :: num_elements
    integer (int_kind),               intent(in)    :: num_forcing_fields

    allocate(this%forcing_fields(num_forcing_fields))
    this%forcing_field_cnt = 0
    this%num_elements      = num_elements
    !TODO: initialize forcing fields to null?

  end subroutine marbl_forcing_fields_constructor

  !*****************************************************************************

  subroutine marbl_forcing_fields_add(this, field_source, marbl_varname, field_units,    &
                                      unit_conv_factor, temporal_interp, field_constant, &
                                      marbl_driver_varname, filename, file_varname,      &
                                      temporal, year_first, year_last, year_align,       &
                                      date, time, marbl_forcing_calendar_name, id)

    class(marbl_forcing_fields_type) , intent(inout) :: this
    character (char_len)             , intent(in)    :: field_source
    character (char_len)             , intent(in)    :: marbl_varname
    character (char_len)             , intent(in)    :: field_units
    real(KIND=r8),           optional, intent(in)    :: unit_conv_factor
    character (char_len),    optional, intent(in)    :: temporal_interp
    real(KIND=r8),           optional, intent(in)    :: field_constant
    character (char_len),    optional, intent(in)    :: marbl_driver_varname
    character (char_len),    optional, intent(in)    :: filename
    character (char_len),    optional, intent(in)    :: file_varname
    character (char_len),    optional, intent(in)    :: temporal
    integer (KIND=int_kind), optional, intent(in)    :: year_first
    integer (KIND=int_kind), optional, intent(in)    :: year_last
    integer (KIND=int_kind), optional, intent(in)    :: year_align
    integer (KIND=int_kind), optional, intent(in)    :: date
    integer (KIND=int_kind), optional, intent(in)    :: time
    type (forcing_monthly_every_ts), optional, target, intent(in) :: marbl_forcing_calendar_name
    integer (KIND=int_kind)          , intent(out)   :: id

    
    integer (KIND=int_kind) :: num_elem

    this%forcing_field_cnt = this%forcing_field_cnt + 1
    id = this%forcing_field_cnt
    if (id .gt. size(this%forcing_fields)) then
      print*, "ERROR: increase max number of forcing fields!"
      ! FIXME: abort
    end if
    num_elem = this%num_elements

    call marbl_single_forcing_field_init(this%forcing_fields(id), &
                                         num_elem, field_source, marbl_varname, &
                                         field_units, unit_conv_factor=unit_conv_factor, &
                                         temporal_interp=temporal_interp,                &
                                         field_constant=field_constant,                  &
                                         marbl_driver_varname=marbl_driver_varname,      &
                                         filename=filename, file_varname=file_varname,   &
                                         temporal=temporal, year_first=year_first,       &
                                         year_last=year_last, year_align=year_align,     &
                                         date=date, time=time,                           &
                                         marbl_forcing_calendar_name=marbl_forcing_calendar_name)
!JW    call this%forcing_fields(id)%initialize(num_elem, field_source, marbl_varname, &
!JW                                            field_units, unit_conv_factor=unit_conv_factor, &
!JW                                            temporal_interp=temporal_interp,                &
!JW                                            field_constant=field_constant,                  &
!JW                                            marbl_driver_varname=marbl_driver_varname,      &
!JW                                            filename=filename, file_varname=file_varname,   &
!JW                                            temporal=temporal, year_first=year_first,       &
!JW                                            year_last=year_last, year_align=year_align,     &
!JW                                            date=date, time=time,                           &
!JW                                            marbl_forcing_calendar_name=marbl_forcing_calendar_name)

  end subroutine marbl_forcing_fields_add

  !*****************************************************************************

  subroutine marbl_forcing_fields_deconstructor(this)

    class(marbl_forcing_fields_type), intent(inout) :: this

    integer (KIND=int_kind) :: n

    do n = 1,size(this%forcing_fields)
      if (allocated(this%forcing_fields(n)%field_data)) then
         deallocate(this%forcing_fields(n)%field_data)
      end if
    end do
    deallocate(this%forcing_fields)

  end subroutine marbl_forcing_fields_deconstructor

  !*****************************************************************************

end module marbl_interface_types
