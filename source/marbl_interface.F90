! -*- mode: f90; indent-tabs-mode: nil; f90-do-indent:3; f90-if-indent:3; f90-type-indent:3; f90-program-indent:2; f90-associate-indent:0; f90-continuation-indent:5  -*-
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

  use marbl_kinds_mod           , only : r8, log_kind, int_kind, log_kind
  use marbl_interface_constants , only : marbl_status_ok

  use marbl_sizes               , only : ecosys_tracer_cnt
  use marbl_sizes               , only : ecosys_ciso_tracer_cnt
  use marbl_sizes               , only : ecosys_used_tracer_cnt
  use marbl_sizes               , only : ecosys_ind_beg, ecosys_ind_end
  use marbl_sizes               , only : ecosys_ciso_ind_beg, ecosys_ciso_ind_end
  use marbl_sizes               , only : autotroph_cnt
  use marbl_sizes               , only : zooplankton_cnt

  use marbl_interface_types     , only : marbl_domain_type
  use marbl_interface_types     , only : marbl_gcm_state_type
  use marbl_interface_types     , only : marbl_status_type
  use marbl_interface_types     , only : marbl_tracer_metadata_type
  use marbl_interface_types     , only : marbl_forcing_fields_type
  use marbl_interface_types     , only : marbl_forcing_input_type
  use marbl_interface_types     , only : marbl_forcing_output_type
  use marbl_interface_types     , only : marbl_diagnostics_type

  use marbl_interface_types     , only : marbl_domain_type
  use marbl_interface_types     , only : marbl_input_type
  use marbl_interface_types     , only : marbl_output_type

  use exit_mod                  , only : exit_POP  !FIXME
  use exit_mod                  , only : sigAbort  !FIXME

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
  
  type, public :: marbl_interface_class
   contains
     procedure, public :: init             
     procedure, public :: set_interior_forcing     
     procedure, public :: set_surface_forcing
     procedure, public :: shutdown         
  end type marbl_interface_class
  
  private :: init
  private :: set_interior_forcing
  private :: set_surface_forcing
  private :: shutdown

  !***********************************************************************

contains

  !***********************************************************************
  
  subroutine init(this,       &
       ciso_on,               &
       nl_buffer,             &
       marbl_domain,          &
       marbl_tracer_metadata, &
       marbl_input,           &
       marbl_output,          &
       marbl_status)

    use marbl_parms            , only : total_input_forcing_cnt
    use marbl_namelist_mod     , only : marbl_nl_cnt
    use marbl_namelist_mod     , only : marbl_nl_buffer_size
    use marbl_ciso_mod         , only : marbl_ciso_init_nml
    use marbl_ciso_mod         , only : marbl_ciso_init_tracer_metadata
    use ecosys_mod             , only : marbl_init_nml
    use ecosys_mod             , only : marbl_sflux_forcing_fields_init
    use ecosys_mod             , only : marbl_init_tracer_metadata
    use ecosys_diagnostics_mod , only : marbl_diagnostics_init  
    use marbl_share_mod        , only : marbl_seconds_in_year
    
    implicit none

    class     (marbl_interface_class)      , intent(inout) :: this
    logical   (log_kind)                   , intent(in)    :: ciso_on
    character (marbl_nl_buffer_size)       , intent(in)    :: nl_buffer(marbl_nl_cnt)
    real      (r8)                         , intent(in)    :: seconds_in_year
    type      (marbl_domain_type)          , intent(in)    :: marbl_domain 
    type      (marbl_gcm_state_type)       , intent(inout) :: marbl_gcm_state
    type      (marbl_tracer_metadata_type) , intent(inout) :: marbl_tracer_metadata(:)
    type      (marbl_diagnostics_type)     , intent(inout) :: marbl_interior_diags
    type      (marbl_diagnostics_type)     , intent(inout) :: marbl_restore_diags
    type      (marbl_diagnostics_type)     , intent(inout) :: marbl_forcing_diags
    type      (marbl_forcing_fields_type)  , intent(inout) :: marbl_forcing_fields
    type      (marbl_forcing_input_type)   , intent(inout) :: marbl_forcing_input
    type      (marbl_forcing_output_type)  , intent(inout) :: marbl_forcing_output
    type      (marbl_status_type)          , intent(out)   :: marbl_status

    !-----------------------------------------------------------------------

    marbl_seconds_in_year = seconds_in_year

    marbl_status%status = marbl_status_ok
    marbl_status%message = ''

    associate(                                                                &
         num_elements_interior => marbl_domain%num_elements_interior_forcing, &
         num_elements_forcing  => marbl_domain%num_elements_surface_forcing,  &
         num_PAR_subcols       => marbl_domain%num_PAR_subcols,               &
         num_levels            => marbl_domain%km                             &
         )

    !--------------------------------------------------------------------
    ! initialize marbl namelists
    !--------------------------------------------------------------------

    call marbl_init_nml(nl_buffer, marbl_status)

    if (marbl_status%status /= marbl_status_ok) then
       call exit_POP(sigAbort, &
            'ERROR in marbl_init_nml: returned status: "'//marbl_status%message//'"')
    end if

    if (ciso_on) then
       call marbl_ciso_init_nml(nl_buffer, marbl_status)

       if (marbl_status%status /= marbl_status_ok) then
          call exit_POP(sigAbort, &
               'ERROR in marbl_ciso_init_nml: returned status: "'//marbl_status%message//'"')
       end if
    end if

    !--------------------------------------------------------------------
    ! initialize marbl gcm_state 
    !--------------------------------------------------------------------

    call marbl_gcm_state%construct(num_levels, num_PAR_subcols)

    !--------------------------------------------------------------------
    ! initialize marbl tracer metadata 
    !--------------------------------------------------------------------

    call marbl_init_tracer_metadata(marbl_tracer_metadata)

    if (ciso_on) then
       call marbl_ciso_init_tracer_metadata(&
            marbl_tracer_metadata(ecosys_ciso_ind_beg:ecosys_ciso_ind_end))
    end if

    !--------------------------------------------------------------------
    ! Initialize marbl diagnostics
    !--------------------------------------------------------------------

    call marbl_diagnostics_init(                                        &
         ciso_on=ciso_on,                                               &
         marbl_domain=marbl_domain, &
         marbl_tracer_metadata=marbl_tracer_metadata(ecosys_ind_beg:ecosys_ind_end), &
         marbl_interior_diags=marbl_interior_diags,                     &
         marbl_restore_diags=marbl_restore_diags,                       &
         marbl_forcing_diags=marbl_forcing_diags)

    !--------------------------------------------------------------------
    ! initialize marbl surface forcing input and output
    !--------------------------------------------------------------------

    call marbl_forcing_input%construct(              &
         num_elements_forcing,                       &
         num_surface_vals=ecosys_used_tracer_cnt,    &
         num_input_forcings=total_input_forcing_cnt, &  
         ciso_on=ciso_on)

    call marbl_sflux_forcing_fields_init(            &
         num_elements_forcing,                       &
         marbl_forcing_fields)

    call marbl_forcing_output%construct(             &
         num_elements_forcing,                       &
         num_surface_vals=ecosys_used_tracer_cnt)
    
    end associate

  end subroutine init

  !***********************************************************************
  
  subroutine set_interior_forcing(this, &
       ciso_on,                       &
       marbl_domain,                  &
       marbl_gcm_state,               &
       column_restore,                &
       column_dust_flux_in,           &
       column_fesedflux,              &
       column_tracers,                &
       marbl_interior_diags,          &
       marbl_restore_diags,           &
       column_ph_prev_3d,             &
       column_ph_prev_alt_co2_3d,     &
       column_dtracers)

    use ecosys_mod          , only : marbl_set_interior_forcing
    use marbl_ciso_mod      , only : marbl_ciso_set_interior_forcing
    use marbl_internal_types, only : marbl_interior_share_type
    use marbl_internal_types, only : marbl_autotroph_share_type
    use marbl_internal_types, only : marbl_zooplankton_share_type
    use marbl_internal_types, only : marbl_particulate_share_type
    
    implicit none

    class   (marbl_interface_class)  , intent(inout) :: this
    logical (log_kind)               , intent(in)    :: ciso_on   
    type    (marbl_domain_type)      , intent(in)    :: marbl_domain                                
    type    (marbl_gcm_state_type)   , intent(in)    :: marbl_gcm_state
    real    (r8)                     , intent(in)    :: column_restore(ecosys_tracer_cnt, marbl_domain%km) 
    real    (r8)                     , intent(in)    :: column_dust_flux_in
    real    (r8)                     , intent(in)    :: column_fesedflux(marbl_domain%km)
    real    (r8)                     , intent(in)    :: column_tracers(ecosys_used_tracer_cnt, marbl_domain%km) 
    real    (r8)                     , intent(inout) :: column_ph_prev_3d(marbl_domain%km)         
    real    (r8)                     , intent(inout) :: column_ph_prev_alt_co2_3d(marbl_domain%km) 
    type    (marbl_diagnostics_type) , intent(inout) :: marbl_interior_diags
    type    (marbl_diagnostics_type) , intent(inout) :: marbl_restore_diags
    real    (r8)                     , intent(out)   :: column_dtracers(ecosys_used_tracer_cnt, marbl_domain%km) 

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    type(marbl_interior_share_type)     :: marbl_interior_share(marbl_domain%km)
    type(marbl_zooplankton_share_type)  :: marbl_zooplankton_share(zooplankton_cnt, marbl_domain%km)
    type(marbl_autotroph_share_type)    :: marbl_autotroph_share(autotroph_cnt, marbl_domain%km)
    type(marbl_particulate_share_type)  :: marbl_particulate_share
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! first allocate memory for local derived data types containing allocatable arrays
    !-----------------------------------------------------------------------

    call marbl_particulate_share%construct(num_levels=marbl_domain%km)

    call marbl_set_interior_forcing( &
         ciso_on,                    &
         marbl_domain,               &
         marbl_gcm_state,            &
         column_restore,             &
         column_dust_flux_in,        &
         column_fesedflux,           &
         column_tracers,             &
         marbl_interior_diags,       &
         marbl_restore_diags,        &
         column_ph_prev_3d,          &
         column_ph_prev_alt_co2_3d,  &
         column_dtracers,            &
         marbl_interior_share,       &
         marbl_zooplankton_share,    &
         marbl_autotroph_share,      &
         marbl_particulate_share)
    
    !  compute time derivatives for ecosystem carbon isotope state variables
    if (ciso_on) then
       call marbl_ciso_set_interior_forcing(                           &
            marbl_domain,                                              &
            marbl_gcm_state,                                           &
            marbl_interior_share,                                      &
            marbl_zooplankton_share,                                   &
            marbl_autotroph_share,                                     &
            marbl_particulate_share,                                   &
            column_tracers(ecosys_ciso_ind_beg:ecosys_ciso_ind_end,:), &
            marbl_interior_diags,                                      &
            column_dtracers(ecosys_ciso_ind_beg:ecosys_ciso_ind_end,:))
    end if

    !-----------------------------------------------------------------------
    ! deallocate memory for local derived data types containing allocatable arrays
    !-----------------------------------------------------------------------

    call marbl_particulate_share%destruct()

  end subroutine set_interior_forcing
  
  !***********************************************************************
  
  subroutine set_surface_forcing(this, &
       num_elements,                   &
       ciso_on,                        &
       marbl_forcing_input,            &
       marbl_forcing_output,           &
       marbl_forcing_diags)

    use marbl_internal_types  , only : marbl_forcing_share_type
    use marbl_ciso_mod        , only : marbl_ciso_set_surface_forcing
    use ecosys_mod            , only : marbl_set_surface_forcing
    
    implicit none

    class   (marbl_interface_class)     , intent(inout) :: this
    integer (int_kind)                  , intent(in)    :: num_elements
    logical (log_kind)                  , intent(in)    :: ciso_on
    type    (marbl_forcing_input_type)  , intent(in)    :: marbl_forcing_input
    type    (marbl_forcing_output_type) , intent(inout) :: marbl_forcing_output
    type    (marbl_diagnostics_type),     intent(inout) :: marbl_forcing_diags

    !-----------------------------------------------------------------------
    !  local variables
    !-----------------------------------------------------------------------
    type(marbl_forcing_share_type)  :: marbl_forcing_share
    !-----------------------------------------------------------------------

    call marbl_forcing_share%construct(num_elements)

    call marbl_set_surface_forcing( &
         ciso_on,                   &
         num_elements,              &
         marbl_forcing_input,       &
         marbl_forcing_output,      &
         marbl_forcing_share,       &
         marbl_forcing_diags)

    if (ciso_on) then
       call marbl_ciso_set_surface_forcing(                                              &
            num_elements,                                                                &
            ecosys_ciso_tracer_cnt,                                                      &
            marbl_forcing_input%land_mask,                                               &
            marbl_forcing_input%sst,                                                     &
            marbl_forcing_input%d13c,                                                    &
            marbl_forcing_input%d14c,                                                    &
            marbl_forcing_input%d14c_glo_avg,                                            &
            marbl_forcing_input%surface_vals(:,ecosys_ciso_ind_beg:ecosys_ciso_ind_end), &
            marbl_forcing_share,                                                         &
            marbl_forcing_output%stf_module(:,ecosys_ciso_ind_beg:ecosys_ciso_ind_end),  &
            marbl_forcing_diags)
    end if

    call marbl_forcing_share%destruct(num_elements)

  end subroutine set_surface_forcing
  
  !***********************************************************************
  
  subroutine shutdown(this)

    implicit none

    class(marbl_interface_class), intent(inout) :: this

    ! free dynamically allocated memory, etc
    
  end subroutine shutdown

end module marbl_interface
