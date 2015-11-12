module marbl_interface_types
  ! module for definitions of types that are shared between marbl interior and the driver.

  use marbl_kinds_mod, only : c0, r8, log_kind, int_kind, char_len
  use marbl_interface_constants, only : marbl_str_length
  use marbl_share_mod, only : autotroph_cnt, zooplankton_cnt

  use domain_size, only : km

  implicit none

  private

  type, public :: marbl_status_type
     integer :: status
     character(marbl_str_length) :: message
  end type marbl_status_type

  type, public :: marbl_column_domain_type
     logical(log_kind) :: land_mask
     integer(int_kind) :: km ! number of vertical grid cells
     integer(int_kind) :: kmt ! index of ocean floor
     real(r8), allocatable :: dzt(:) ! (km) delta z for partial bottom cells
     real(r8), allocatable :: dz(:) ! (km) delta z
     real(r8), allocatable :: temperature(:) ! (km)
     real(r8), allocatable :: salinity(:) ! (km)
  end type marbl_column_domain_type

  type, public :: marbl_saved_state_type
     ! this struct is necessary because there is some global state
     ! that needs to be preserved for marbl
     real (r8), dimension(:, :, :), allocatable :: dust_FLUX_IN  ! dust flux not stored in STF since dust is not prognostic
     real (r8), dimension(:, :, :), allocatable :: PAR_out  ! photosynthetically available radiation (W/m^2)
     real (r8), dimension(:, :, :, :), allocatable :: PH_PREV_3D         ! computed pH_3D from previous time step
     real (r8), dimension(:, :, :, :), allocatable :: PH_PREV_ALT_CO2_3D ! computed pH_3D from previous time step, alternative CO2
     logical (log_kind), dimension(:, :, :), allocatable :: LAND_MASK

     
  end type marbl_saved_state_type
  
  ! marbl_diagnostics : used to pass diagnostic information for tavg
  ! etc from marbl back to the driver. The driver is responsible for
  ! sizing these arrays correctly according to the size of the domain
  ! and tracer count returned by marbl.
  type :: marbl_diagnostic_data_and_metadata_type
    character(len=char_len) :: long_name
    character(len=char_len) :: short_name
    character(len=char_len) :: units
    logical(log_kind) :: compute_now
    character(len=char_len) :: vertical_grid ! 'none', 'layer_avg', 'layer_iface'
    logical(log_kind) :: ltruncated_vertical_extent
    real(r8) :: field_2d
    real(r8), allocatable, dimension(:) :: field_3d

  contains
    procedure, public :: initialize  => marbl_diagnostic_metadata_init
  end type marbl_diagnostic_data_and_metadata_type

  type, public :: marbl_diagnostics_type
     type(marbl_diagnostic_data_and_metadata_type), dimension(:), allocatable :: diags

  contains
    procedure, public :: construct      => marbl_diagnostics_constructor
    procedure, public :: initialize     => marbl_diagnostics_init
    procedure, public :: add_diagnostic => marbl_diagnostics_add
    procedure, public :: set_to_zero    => marbl_diagnostics_set_to_zero
    procedure, public :: deconstruct    => marbl_diagnostics_deconstructor
  end type marbl_diagnostics_type

  type, public :: carbonate_type
     ! FIXME(bja, 2015-07) remove alt_co2 variables, and just reuse
     ! the type as type(column_carbonate_type) :: carbonate,
     ! carbonate_alt
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
     real (r8) :: VNC             ! C-specific N uptake rate (mmol N/mmol C/sec)
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

  type, public :: photosynthetically_available_radiation_type
     real(r8) :: in     ! photosynthetically available radiation (W/m^2)
     real(r8) :: KPARdz ! PAR adsorption coefficient (non-dim)
     real(r8) :: avg    ! average PAR over mixed layer depth (W/m^2)
  end type photosynthetically_available_radiation_type

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

  type, public :: dissolved_organic_matter_type
     real (r8) :: DOC_prod         ! production of DOC (mmol C/m^3/sec)
     real (r8) :: DOC_remin        ! remineralization of DOC (mmol C/m^3/sec)
     real (r8) :: DON_prod         ! production of dissolved organic N
     real (r8) :: DON_remin        ! portion of DON remineralized
     real (r8) :: DOFe_prod        ! produciton of dissolved organic Fe
     real (r8) :: DOFe_remin       ! portion of DOFe remineralized
     real (r8) :: DOP_prod         ! production of dissolved organic P
     real (r8) :: DOP_remin        ! portion of DOP remineralized
     real (r8) :: DONr_remin       ! portion of refractory DON remineralized
     real (r8) :: DOPr_remin       ! portion of refractory DOP remineralized
  end type dissolved_organic_matter_type

 integer, parameter :: max_diags = 66 + autotroph_cnt*26 + zooplankton_cnt*8
 integer, public :: diag_cnt

contains

  subroutine marbl_diagnostics_add(this, lname, sname, units, vgrid,          &
                                   truncate, id)

    class(marbl_diagnostics_type), intent(inout) :: this
    character(len=char_len), intent(in) :: lname, sname, units, vgrid
    logical, intent(in) :: truncate
    integer, intent(out) :: id


    diag_cnt = diag_cnt + 1
    id = diag_cnt
    if (id.gt.size(this%diags)) then
      print*, "ERROR: increase max number of diagnostics!"
      ! FIXME: abort
    end if
    call this%diags(id)%initialize(lname, sname, units, vgrid, truncate)

  end subroutine marbl_diagnostics_add

  subroutine marbl_diagnostic_metadata_init(this, lname, sname, units, vgrid, &
                                            truncate)

    class(marbl_diagnostic_data_and_metadata_type), intent(inout) :: this
    character(len=char_len), intent(in) :: lname, sname, units, vgrid
    logical, intent(in) :: truncate

    character(len=char_len), dimension(3) :: valid_vertical_grids
    integer :: n
    logical :: valid_vgrid

    valid_vertical_grids(1) = 'none'
    valid_vertical_grids(2) = 'layer_avg'
    valid_vertical_grids(3) = 'layer_iface'
    valid_vgrid = .false.
    do n=1,size(valid_vertical_grids)
      if (trim(vgrid).eq.trim(valid_vertical_grids(n))) valid_vgrid = .true.
    end do
    if (.not.valid_vgrid) then
      print*, "ERROR: ", trim(vgrid), " is not a valid vertical grid for MARBL"
      ! FIXME: abort
    end if

    this%compute_now = .true.
    this%long_name = trim(lname)
    this%short_name = trim(sname)
    this%units = trim(units)
    this%vertical_grid = trim(vgrid)
    this%ltruncated_vertical_extent = truncate

    ! Allocate column memory for 3D variables
    if (trim(vgrid).eq.'layer_avg') then
      allocate(this%field_3d(km))
    end if

    if (trim(vgrid).eq.'layer_iface') then
      allocate(this%field_3d(km+1))
    end if

  end subroutine marbl_diagnostic_metadata_init

  subroutine marbl_diagnostics_constructor(this, ecosys_tracer_cnt)

    class(marbl_diagnostics_type), intent(inout) :: this
    integer, intent(in) :: ecosys_tracer_cnt

    allocate(this%diags(max_diags))
    diag_cnt = 0

    call this%initialize()

  end subroutine marbl_diagnostics_constructor

  subroutine marbl_diagnostics_init(this)

    class(marbl_diagnostics_type), intent(inout) :: this

    call this%set_to_zero()

  end subroutine marbl_diagnostics_init

  subroutine marbl_diagnostics_set_to_zero(this)

    class(marbl_diagnostics_type), intent(inout) :: this

    integer :: m,n

    do n=1,size(this%diags) ! ecosys_diag_cnt_2d
      this%diags(n)%field_2d = c0
      if (allocated(this%diags(n)%field_3d)) then
        this%diags(n)%field_3d(:) = c0
      end if
    end do

  end subroutine marbl_diagnostics_set_to_zero

  subroutine marbl_diagnostics_deconstructor(this)

    class(marbl_diagnostics_type), intent(inout) :: this

    integer :: n

    do n=1,size(this%diags) ! ecosys_diag_cnt_2d
      if (allocated(this%diags(n)%field_3d)) then
        deallocate(this%diags(n)%field_3d)
      end if
    end do
    deallocate(this%diags)

  end subroutine marbl_diagnostics_deconstructor

end module marbl_interface_types
