module marbl_interface_types
  ! module for definitions of types that are shared between marbl interior and the driver.

  use marbl_kinds_mod, only : r8, log_kind, int_kind
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
  type, public :: ecosys_diagnostics_type
     ! NOTE(bja, 2015-07) this is the old slab ordering, will either
     ! go away or move to ecosys_driver!
     real(r8), allocatable :: DIAGS(:, :, :)         ! (nx_block, ny_block, ecosys_diag_cnt)
     real(r8), allocatable :: AUTO_DIAGS(:, :, :, :) ! (nx_block, ny_block, autotroph_cnt, auto_diag_cnt)
     real(r8), allocatable :: ZOO_DIAGS(:, :, :, :)  ! (nx_block, ny_block, zooplankton_cnt, zoo_diag_cnt)
     real(r8), allocatable :: PART_DIAGS(:, :, :)    ! (nx_block, ny_block, part_diag_cnt)
     real(r8), allocatable :: restore_diags(:, :, :) ! (nx_block, ny_block, ecosys_diag_cnt)
  end type ecosys_diagnostics_type

  type, public :: marbl_diagnostics_type
     real(r8), allocatable :: diags(:)           ! (ecosys_diag_cnt)
     real(r8), allocatable :: auto_diags(:, :) ! (auto_diag_cnt, autotroph_cnt)
     real(r8), allocatable :: zoo_diags(:, :)  ! (zoo_diag_cnt, zooplankton_cnt)
     real(r8), allocatable :: part_diags(:)      ! (part_diag_cnt)
     real(r8), allocatable :: restore_diags(:) ! (ecosys_diag_cnt, ny_block, nx_block)
  end type marbl_diagnostics_type
  
end module marbl_interface_types
