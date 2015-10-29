module marbl_interface_types
  ! module for definitions of types that are shared between marbl interior and the driver.

  use marbl_kinds_mod, only : c0, r8, log_kind, int_kind
  use marbl_interface_constants, only : marbl_str_length

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
  type, public :: marbl_diagnostics_type
     real(r8), allocatable :: diags_2d(:,:)        ! (km, ecosys_diag_cnt_2d)
     real(r8), allocatable :: diags_3d(:,:)        ! (km, ecosys_diag_cnt_3d)
     real(r8), allocatable :: auto_diags(:, :, :) ! (km, auto_diag_cnt, autotroph_cnt)
     real(r8), allocatable :: zoo_diags(:, :, :)  ! (km, zoo_diag_cnt, zooplankton_cnt)
     real(r8), allocatable :: part_diags_2d(:,:)      ! (km, part_diag_cnt_2d)
     real(r8), allocatable :: part_diags_3d(:,:)      ! (km, part_diag_cnt_3d)
     real(r8), allocatable :: restore_diags(:,:) ! (km, ecosys_tracer_cnt)

  contains
    procedure, public :: construct   => marbl_diagnostics_constructor
    procedure, public :: set_to_zero => marbl_diagnostics_set_to_zero
    procedure, public :: deconstruct => marbl_diagnostics_deconstructor
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

contains

  subroutine marbl_diagnostics_constructor(this, km, ecosys_diag_cnt_2d,      &
           ecosys_diag_cnt_3d, auto_diag_cnt, zoo_diag_cnt, part_diag_cnt_2d, &
           part_diag_cnt_3d, ecosys_tracer_cnt, autotroph_cnt, zooplankton_cnt)

    class(marbl_diagnostics_type), intent(inout) :: this
    integer, intent(in) :: km
    integer, intent(in) :: ecosys_diag_cnt_2d, ecosys_diag_cnt_3d,            &
                           auto_diag_cnt, zoo_diag_cnt, part_diag_cnt_2d,     &
                           part_diag_cnt_3d
    integer, intent(in) :: ecosys_tracer_cnt, autotroph_cnt, zooplankton_cnt

    allocate(this%diags_2d(km, ecosys_diag_cnt_2d))
    allocate(this%diags_3d(km, ecosys_diag_cnt_3d))
    allocate(this%auto_diags(km, auto_diag_cnt, autotroph_cnt))
    allocate(this%zoo_diags(km, zoo_diag_cnt, zooplankton_cnt))
    allocate(this%part_diags_2d(km, part_diag_cnt_2d))
    allocate(this%part_diags_3d(km, part_diag_cnt_3d))
    allocate(this%restore_diags(km, ecosys_tracer_cnt))

    call this%set_to_zero()

  end subroutine marbl_diagnostics_constructor

  subroutine marbl_diagnostics_set_to_zero(this)

    class(marbl_diagnostics_type), intent(inout) :: this

    this%diags_2d = c0
    this%diags_3d = c0
    this%auto_diags = c0
    this%zoo_diags = c0
    this%part_diags_2d = c0
    this%part_diags_3d = c0
    this%restore_diags = c0

  end subroutine marbl_diagnostics_set_to_zero

  subroutine marbl_diagnostics_deconstructor(this)

    class(marbl_diagnostics_type), intent(inout) :: this

    deallocate(this%diags_2d)
    deallocate(this%diags_3d)
    deallocate(this%auto_diags)
    deallocate(this%zoo_diags)
    deallocate(this%part_diags_2d)
    deallocate(this%part_diags_3d)
    deallocate(this%restore_diags)

  end subroutine marbl_diagnostics_deconstructor

end module marbl_interface_types
