module marbl_interface
  !
  ! marbl interface
  !
  ! This module contains the public interface to marbl. This is the
  ! only public API and the only interface guaranteed to be documented
  ! and semi-stable. All other marbl rountines are private and may
  ! change at any time.
  !
  ! The following terminology is used in this module:
  !
  ! driver - refers to the code that is calling marbl routines. This
  !          can be anything from a full GCM to a simplified single
  !          column test system
  !
  ! Notes :
  !
  ! (1) ALL marbl routines will use marbl_share_mod - but not ANY GCM routines
  ! as an example - so will pass marbl_private_data%marbl_interio_shar to 
  ! the underlying marbl routines RATHER than marble_private_data directly
  ! this will avoid circular dependencies 
  !
  ! (2) In GCM driver, (e.g.see ecosys_driver), will ONLY have the following 
  ! use statements to marbl at the module level to create data structures
  !      use marbl_interface, only: marbl_private_data_type
  !      .....
  !      type(marbl_private_data_type) :: marbl_private_data
  !      ....
  !
  ! (3) In the subroutine for init (e.g. ecosys_driver_init) will call marbl_init
  ! as follows  
  !      use marbl_interface, only: marbl_init
  !      .....
  !      call marbl_init(marbl_private_data, ...)
  !
  ! (4) In the run subroutines (e.g. ecosys_driver_set_interior) will call
  ! marbl_set_interior as follows
  !      use marbl_interface, only: marbl_set_interior
  !      .....
  !      call marbl_set_interior(marbl_private_data, ....)

  use marbl_kinds_mod           , only : int_kind
  use marbl_interface_constants , only : marbl_status_ok
  
  use exit_mod, only : exit_POP
  use exit_mod, only : sigAbort

  implicit none

  private

  ! NOTE: marbl_private_data is storage that the driver maintains for
  ! marbl to avoid global variables. It is NOT part of the public API
  ! to the marbl library; the driver should consider it to be an
  ! opaque blob of memory. Do not access any data in this
  ! structure. The contents may change at any time between api
  ! revisions, and the data may change at any time during a simulation!

  !-----------------------------------------------------------------------------
  !
  ! The following data structures are part of the public API that the
  ! driver will read/write to interact with marbl.
  !
  !-----------------------------------------------------------------------------
  
  ! marbl_sizes : should be set during marbl_init and returned to the
  ! driver as an intent(out) so it can check that its memory
  ! allocation agrees with marbl's.
  type, public :: marbl_sizes_type
     integer :: ecosys_tracer_cnt
     integer :: autotroph_cnt
     integer :: zooplankton_cnt
     integer :: ecosys_ciso_tracer_cnt
  end type marbl_sizes_type

  type, public :: marbl_interface_class
     ! FIXME(bja, 2015-01) needs to private when all data is isolated!
     type(marbl_sizes_type), private :: marbl_sizes

   contains
     procedure, public :: init => marbl_init
     procedure, public :: shutdown => marbl_shutdown
     procedure, public :: set_interior => marbl_set_interior
     procedure, public :: set_surface_flux => marbl_set_surface_flux
  end type marbl_interface_class
  
  private :: &
       marbl_init, &
       marbl_shutdown, &
       marbl_set_interior, &
       marbl_set_surface_flux

contains

  !-----------------------------------------------------------------------------
  
  subroutine marbl_init(this, marbl_sizes, nl_buffer, marbl_tracer_metadata, marbl_status )

    use marbl_namelist_mod        , only: marbl_nl_cnt
    use marbl_namelist_mod        , only: marbl_nl_buffer_size
    use marbl_interface_types     , only: marbl_status_type
    use marbl_interface_types     , only: marbl_tracer_metadata_type
    use marbl_share_mod           , only: autotroph_cnt, zooplankton_cnt
    use marbl_share_mod           , only: ecosys_tracer_cnt
    use marbl_share_mod           , only: ecosys_ciso_tracer_cnt
    use ecosys_mod                , only: marbl_ecosys_init_nml
    use ecosys_mod                , only: marbl_ecosys_init_tracer_metadata
    
    implicit none

    class(marbl_interface_class)    , intent(inout) :: this
    character(marbl_nl_buffer_size), dimension(marbl_nl_cnt), intent(in) :: nl_buffer
    type(marbl_tracer_metadata_type), intent(inout) :: marbl_tracer_metadata(:)
    type(marbl_sizes_type)          , intent(out)   :: marbl_sizes
    type(marbl_status_type)         , intent(out)   :: marbl_status

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: bid ! local block address for this block

    !-----------------------------------------------------------------------

    ! now we know how many tracers marbl has, we can verify that pop
    ! has the correctly sized data.
    marbl_sizes%ecosys_tracer_cnt      = ecosys_tracer_cnt
    marbl_sizes%autotroph_cnt          = autotroph_cnt
    marbl_sizes%zooplankton_cnt        = zooplankton_cnt
    marbl_sizes%ecosys_ciso_tracer_cnt = ecosys_ciso_tracer_cnt
    
    ! initialize marbl namelists
    call marbl_ecosys_init_nml(nl_buffer, marbl_status)

    if (marbl_status%status /= marbl_status_ok) then
       call exit_POP(sigAbort, &
            'ERROR in marbl_init: returned status: "'//marbl_status%message//'"')
    end if

    ! initialize marbl tracer metadata 
    call marbl_ecosys_init_tracer_metadata(marbl_tracer_metadata)

  end subroutine marbl_init

  !-----------------------------------------------------------------------------
  
  subroutine marbl_shutdown(this)

    implicit none

    class(marbl_interface_class), intent(inout) :: this

    ! free dynamically allocated memory, etc
    
  end subroutine marbl_shutdown

  !-----------------------------------------------------------------------------
  
  subroutine marbl_set_interior(this)

    use ecosys_mod     , only: marbl_ecosys_set_interior
    
    implicit none

    class(marbl_interface_class), intent(inout) :: this

    ! unpack marbl_privat_data and pass the contents through individually.
!!$    call ecosys_set_interior(marbl_private_data%ecosys_interior_share)
    
  end subroutine marbl_set_interior
  
  !-----------------------------------------------------------------------------
  
  subroutine marbl_set_surface_flux(this)

    use ecosys_mod     , only: marbl_ecosys_set_sflux
    
    implicit none

    class(marbl_interface_class), intent(inout) :: this

    ! unpack marbl_privat_data and pass the contents through individually.
!!$    call ecosys_set_sflux(marbl_private_data%ecosys_surface_share)
    
  end subroutine marbl_set_surface_flux
  
end module marbl_interface
