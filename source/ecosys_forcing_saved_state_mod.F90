module ecosys_forcing_saved_state_mod

  !-----------------------------------------------------------------------
  ! forcing related saved state variables
  !-----------------------------------------------------------------------

  use kinds_mod, only : int_kind, log_kind, r8

  implicit none
  private

  public :: ecosys_forcing_saved_state_init
  public :: ecosys_forcing_saved_state_get_var_val
  public :: ecosys_forcing_saved_state_update

  !-----------------------------------------------------------------------
  ! module public variables
  !-----------------------------------------------------------------------

  logical (log_kind),            public :: lbox_atm_co2            ! is box_atm_trace_gas_co2 being simulated
                                                                   ! set in ecosys_forcing_mod.F90

  real (r8),                     public :: box_atm_co2_init_val    ! initial value for box_atm_trace_gas_co2
                                                                   ! set in ecosys_forcing_mod.F90

  integer (int_kind), parameter, public :: box_atm_co2_forcing_saved_state_id = 1

  !-----------------------------------------------------------------------
  ! module private variables
  !-----------------------------------------------------------------------

  integer (int_kind) :: box_atm_co2_atm_trace_gas_index ! index returned from box_atm_trace_gas_mod

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------

  subroutine ecosys_forcing_saved_state_init(ecosys_restart_filename)

    use box_atm_trace_gas_mod, only : box_atm_trace_gas_define_var
    use box_atm_trace_gas_mod, only : box_atm_trace_gas_var_exists_in_file
    use box_atm_trace_gas_mod, only : box_atm_trace_gas_init_var

    character (len=*), intent(in) :: ecosys_restart_filename

    !-----------------------------------------------------------------------

    if (lbox_atm_co2) then
      call box_atm_trace_gas_define_var('ecosys_forcing_atm_co2', 1.0e-9_r8, 1.0e6_r8, box_atm_co2_atm_trace_gas_index)
      if (len_trim(ecosys_restart_filename) > 0) then
        if (box_atm_trace_gas_var_exists_in_file(box_atm_co2_atm_trace_gas_index, ecosys_restart_filename)) then
          call box_atm_trace_gas_init_var(box_atm_co2_atm_trace_gas_index, ecosys_restart_filename)
        else
          call box_atm_trace_gas_init_var(box_atm_co2_atm_trace_gas_index, box_atm_co2_init_val)
        endif
      else
        call box_atm_trace_gas_init_var(box_atm_co2_atm_trace_gas_index, box_atm_co2_init_val)
      end if
    end if

    !-----------------------------------------------------------------------

  end subroutine ecosys_forcing_saved_state_init

  !-----------------------------------------------------------------------

  function ecosys_forcing_saved_state_get_var_val(varid) result(FIELD)

    use blocks, only : nx_block, ny_block
    use box_atm_trace_gas_mod, only : box_atm_trace_gas_get_var_val

    integer (int_kind), intent(in) :: varid
    real (r8)                      :: FIELD(nx_block, ny_block)

    if (varid == box_atm_co2_forcing_saved_state_id) then
      FIELD(:,:) = box_atm_trace_gas_get_var_val(box_atm_co2_atm_trace_gas_index)
    end if

  end function ecosys_forcing_saved_state_get_var_val

  !-----------------------------------------------------------------------

  subroutine ecosys_forcing_saved_state_update(flux_co2)

    use box_atm_trace_gas_mod, only : box_atm_trace_gas_update_var

    real (r8), intent(in) :: flux_co2(:,:,:)

    if (lbox_atm_co2) then
      call box_atm_trace_gas_update_var(box_atm_co2_atm_trace_gas_index, flux_co2)
    end if

  end subroutine ecosys_forcing_saved_state_update

  !-----------------------------------------------------------------------

end module ecosys_forcing_saved_state_mod
