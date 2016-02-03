module marbl_interface_types
  ! module for definitions of types that are shared between marbl interior and the driver.

  use marbl_kinds_mod, only : c0, r8, log_kind, int_kind, char_len
  use marbl_interface_constants, only : marbl_str_length
  use ecosys_constants, only : autotroph_cnt, zooplankton_cnt, ecosys_tracer_cnt

  implicit none

  private

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
  type, public :: marbl_column_domain_type
     logical(log_kind)     :: land_mask
     integer(int_kind)     :: PAR_nsubcols      ! number of sub-column values for PAR
     integer(int_kind)     :: km                ! number of vertical grid cells
     integer(int_kind)     :: kmt               ! index of ocean floor
     real(r8), allocatable :: zt(:)             ! (km) vert dist from sfc to midpoint of layer
     real(r8), allocatable :: zw(:)             ! (km) vert dist from sfc to bottom of layer
     real(r8), allocatable :: delta_z(:)        ! (km) delta z - different values for partial bottom cells
     real(r8), allocatable :: dz(:)             ! (km) delta z - same values for partial bottom cells
     real(r8), allocatable :: PAR_col_frac(:)   ! column fraction occupied by each sub-column
     real(r8), allocatable :: surf_shortwave(:) ! surface shortwave for each sub-column (W/m^2)
  end type marbl_column_domain_type

  !*****************************************************************************
  ! FIXME(mnl,2016-01) move PAR_col_frac and surf_shortwave into this datatype
  !                    and come up with better name
  type, public :: marbl_gcm_state_type
     real(r8), allocatable :: temperature(:) ! (km)
     real(r8), allocatable :: salinity(:)    ! (km)
  end type marbl_gcm_state_type

  !*****************************************************************************
  type :: marbl_single_diagnostic_type
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

  !*****************************************************************************
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
  type, public :: carbonate_type
     ! FIXME(bja, 2015-07) remove alt_co2 variables, and just reuse
     ! the type as type(column_carbonate_type) :: carbonate, carbonate_alt
     real (r8) :: CO3 ! carbonate ion
     real (r8) :: HCO3 ! bicarbonate ion
     real (r8) :: H2CO3 ! carbonic acid
     real (r8) :: pH
     real (r8) :: CO3_sat_calcite
     real (r8) :: CO3_sat_aragonite
     real (r8) :: CO3_ALT_CO2 ! carbonate ion, alternative CO2
     real (r8) :: HCO3_ALT_CO2 ! bicarbonate ion, alternative CO2
     real (r8) :: H2CO3_ALT_CO2  ! carbonic acid, alternative CO2
     real (r8) :: pH_ALT_CO2
  end type carbonate_type

  !*****************************************************************************
  type, public :: autotroph_secondary_species_type
     real (r8) :: thetaC          ! local Chl/C ratio (mg Chl/mmol C)
     real (r8) :: QCaCO3          ! CaCO3/C ratio (mmol CaCO3/mmol C)
     real (r8) :: Qfe             ! init fe/C ratio (mmolFe/mmolC)
     real (r8) :: gQfe            ! fe/C for growth
     real (r8) :: Qsi             ! initial Si/C ratio (mmol Si/mmol C)
     real (r8) :: gQsi            ! diatom Si/C ratio for growth (new biomass)
     real (r8) :: VNO3            ! NH4 uptake rate (non-dim)
     real (r8) :: VNH4            ! NO3 uptake rate (non-dim)
     real (r8) :: VNtot           ! total N uptake rate (non-dim)
     real (r8) :: NO3_V           ! nitrate uptake (mmol NO3/m^3/sec)
     real (r8) :: NH4_V           ! ammonium uptake (mmol NH4/m^3/sec)
     real (r8) :: PO4_V           ! PO4 uptake (mmol PO4/m^3/sec)
     real (r8) :: DOP_V           ! DOP uptake (mmol DOP/m^3/sec)
     real (r8) :: VPO4            ! C-specific PO4 uptake (non-dim)
     real (r8) :: VDOP            ! C-specific DOP uptake rate (non-dim)
     real (r8) :: VPtot           ! total P uptake rate (non-dim)
     real (r8) :: f_nut           ! nut limitation factor, modifies C fixation (non-dim)
     real (r8) :: VFe             ! C-specific Fe uptake (non-dim)
     real (r8) :: VSiO3           ! C-specific SiO3 uptake (non-dim)
     real (r8) :: light_lim       ! light limitation factor
     real (r8) :: PCphoto         ! C-specific rate of photosynth. (1/sec)
     real (r8) :: photoC          ! C-fixation (mmol C/m^3/sec)
     real (r8) :: photoFe         ! iron uptake
     real (r8) :: photoSi         ! silicon uptake (mmol Si/m^3/sec)
     real (r8) :: photoacc        ! Chl synth. term in photoadapt. (GD98) (mg Chl/m^3/sec)
     real (r8) :: auto_loss       ! autotroph non-grazing mort (mmol C/m^3/sec)
     real (r8) :: auto_loss_poc   ! auto_loss routed to poc (mmol C/m^3/sec)
     real (r8) :: auto_loss_doc   ! auto_loss routed to doc (mmol C/m^3/sec)
     real (r8) :: auto_loss_dic   ! auto_loss routed to dic (mmol C/m^3/sec)
     real (r8) :: auto_agg        ! autotroph aggregation (mmol C/m^3/sec)
     real (r8) :: auto_graze      ! autotroph grazing rate (mmol C/m^3/sec)
     real (r8) :: auto_graze_zoo  ! auto_graze routed to zoo (mmol C/m^3/sec)
     real (r8) :: auto_graze_poc  ! auto_graze routed to poc (mmol C/m^3/sec)
     real (r8) :: auto_graze_doc  ! auto_graze routed to doc (mmol C/m^3/sec)
     real (r8) :: auto_graze_dic  ! auto_graze routed to dic (mmol C/m^3/sec)
     real (r8) :: Pprime          ! used to limit autotroph mort at low biomass (mmol C/m^3)
     real (r8) :: CaCO3_PROD      ! prod. of CaCO3 by small phyto (mmol CaCO3/m^3/sec)
     real (r8) :: Nfix            ! total Nitrogen fixation (mmol N/m^3/sec)
     real (r8) :: Nexcrete        ! fixed N excretion
     real (r8) :: remaining_P_dop ! remaining_P from mort routed to DOP pool
     real (r8) :: remaining_P_dip ! remaining_P from mort routed to remin
  end type autotroph_secondary_species_type

  !*****************************************************************************
  type, public :: zooplankton_secondary_species_type
     real (r8):: f_zoo_detr       ! frac of zoo losses into large detrital pool (non-dim)
     real (r8):: x_graze_zoo      ! {auto, zoo}_graze routed to zoo (mmol C/m^3/sec)
     real (r8):: zoo_graze        ! zooplankton losses due to grazing (mmol C/m^3/sec)
     real (r8):: zoo_graze_zoo    ! grazing of zooplankton routed to zoo (mmol C/m^3/sec)
     real (r8):: zoo_graze_poc    ! grazing of zooplankton routed to poc (mmol C/m^3/sec)
     real (r8):: zoo_graze_doc    ! grazing of zooplankton routed to doc (mmol C/m^3/sec)
     real (r8):: zoo_graze_dic    ! grazing of zooplankton routed to dic (mmol C/m^3/sec)
     real (r8):: zoo_loss         ! mortality & higher trophic grazing on zooplankton (mmol C/m^3/sec)
     real (r8):: zoo_loss_poc     ! zoo_loss routed to poc (mmol C/m^3/sec)
     real (r8):: zoo_loss_doc     ! zoo_loss routed to doc (mmol C/m^3/sec)
     real (r8):: zoo_loss_dic     ! zoo_loss routed to dic (mmol C/m^3/sec)
     real (r8):: Zprime           ! used to limit zoo mort at low biomass (mmol C/m^3)
  end type zooplankton_secondary_species_type

  !*****************************************************************************
  type, public :: photosynthetically_available_radiation_type
     real(r8), allocatable :: col_frac(:)    ! column fraction occupied by each sub-column, dimension is (PAR_nsubcols)
     real(r8), allocatable :: interface(:,:) ! PAR at layer interfaces, dimensions are (0:km,PAR_nsubcols)
     real(r8), allocatable :: avg(:,:)       ! PAR averaged over layer, dimensions are (km,PAR_nsubcols)
     real(r8), allocatable :: KPARdz(:)      ! PAR adsorption coefficient times dz (cm), dimension is (km)
  end type photosynthetically_available_radiation_type

  !*****************************************************************************
  type, public :: dissolved_organic_matter_type
     real (r8) :: DOC_prod         ! production of DOC (mmol C/m^3/sec)
     real (r8) :: DOC_remin        ! remineralization of DOC (mmol C/m^3/sec)
     real (r8) :: DOCr_remin       ! remineralization of DOCr
     real (r8) :: DON_prod         ! production of DON
     real (r8) :: DON_remin        ! remineralization of DON
     real (r8) :: DONr_remin       ! remineralization of DONr
     real (r8) :: DOP_prod         ! production of DOP
     real (r8) :: DOP_remin        ! remineralization of DOP
     real (r8) :: DOPr_remin       ! remineralization of DOPr
  end type dissolved_organic_matter_type

  !*****************************************************************************

contains

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

end module marbl_interface_types
