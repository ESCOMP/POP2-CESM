module box_atm_trace_gas_mod

  !-----------------------------------------------------------------------
  ! Provide an interface to simulate single-box atm trace gases.
  !-----------------------------------------------------------------------

  use kinds_mod, only: char_len, int_kind, log_kind, r8
  use io_tools, only: document
  use exit_mod, only: exit_POP, sigAbort

  implicit none
  private

  public :: box_atm_trace_gas_define_var
  public :: box_atm_trace_gas_get_var_index
  public :: box_atm_trace_gas_var_exists_in_file
  public :: box_atm_trace_gas_init_var
  public :: box_atm_trace_gas_update_var
  public :: box_atm_trace_gas_get_var_val
  public :: box_atm_trace_gas_write_restart

  !-----------------------------------------------------------------------
  ! module private types and variables
  !-----------------------------------------------------------------------

  type :: box_atm_trace_gas_type
    character (char_len) :: name             ! user supplied trace gas name
    real (r8)            :: flux_conv_factor ! conversion factor to convert flux to [mol/cm2/s]
                                             !    (e.g. 1e-9 for nmol/cm2/s)
    real (r8)            :: vmr_conv_factor  ! conversion factor to convert internal vmr units [mol/mol]
                                             !    to user vmr units (e.g. 1e6 for ppmv)
    character (char_len) :: file_varname     ! name of variable in restart file
    logical (log_kind)   :: linit            ! has val been initialized?
    real (r8)            :: val
  end type box_atm_trace_gas_type

  integer (int_kind), parameter :: box_atm_trace_gas_cnt_max = 100
  integer (int_kind) :: box_atm_trace_gas_cnt = 0

  type (box_atm_trace_gas_type), target :: &
       box_atm_trace_gas_array(box_atm_trace_gas_cnt_max) ! array of box_atm_trace_gas variables

  !-----------------------------------------------------------------------
  ! generic interface definitions
  !-----------------------------------------------------------------------

  interface box_atm_trace_gas_init_var
    module procedure box_atm_trace_gas_init_var_filename
    module procedure box_atm_trace_gas_init_var_val
  end interface

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------

  subroutine box_atm_trace_gas_define_var(name, flux_conv_factor, vmr_conv_factor, index)

    !-----------------------------------------------------------------------
    ! Define a box_atm_trace_gas var. It is a fatal error to attempt to define a
    ! previously defined name.
    !
    ! This should only be called once per task.
    !-----------------------------------------------------------------------

    use constants, only: c0

    character (len=*), intent(in)   :: name             ! see box_atm_trace_gas_type for description
    real (r8), intent(in)           :: flux_conv_factor ! see box_atm_trace_gas_type for description
    real (r8), intent(in)           :: vmr_conv_factor  ! see box_atm_trace_gas_type for description
    integer (int_kind), intent(out) :: index            ! returned index

    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------

    character (len=*), parameter :: subname = 'box_atm_trace_gas_mod:box_atm_trace_gas_define_var'

    !-----------------------------------------------------------------------

    call document(subname, 'name', name)

    !-----------------------------------------------------------------------
    ! error checking
    !-----------------------------------------------------------------------

    ! avoid spaces in name, these cause problems for creating restart file
    ! variable names based on name
    if (scan(trim(name), ' ') > 0) then
      call exit_POP(sigAbort, 'spaces are not allowed in box_atm_trace_gas names')
    end if

    ! check to see if name is already registered
    call box_atm_trace_gas_get_var_index(name, index, exit_on_err=.false.)
    if (index > 0) then
      call exit_POP(sigAbort, 'box_atm_trace_gas name already defined')
    end if

    box_atm_trace_gas_cnt = box_atm_trace_gas_cnt + 1
    if (box_atm_trace_gas_cnt > box_atm_trace_gas_cnt_max) then
      call document(subname, 'box_atm_trace_gas_cnt_max', box_atm_trace_gas_cnt_max)
      call exit_POP(sigAbort, 'too many box_atm_trace_gas variables defined')
    end if

    !-----------------------------------------------------------------------
    ! setup new box_atm_trace_gas var
    !-----------------------------------------------------------------------

    box_atm_trace_gas_array(box_atm_trace_gas_cnt)%name = trim(name)
    box_atm_trace_gas_array(box_atm_trace_gas_cnt)%flux_conv_factor = flux_conv_factor
    box_atm_trace_gas_array(box_atm_trace_gas_cnt)%vmr_conv_factor = vmr_conv_factor
    box_atm_trace_gas_array(box_atm_trace_gas_cnt)%file_varname = 'box_atm_trace_gas_' /&
         &/ trim(name)
    box_atm_trace_gas_array(box_atm_trace_gas_cnt)%linit = .false.
    box_atm_trace_gas_array(box_atm_trace_gas_cnt)%val = c0

    index = box_atm_trace_gas_cnt

    !-----------------------------------------------------------------------

  end subroutine box_atm_trace_gas_define_var

  !-----------------------------------------------------------------------

  subroutine box_atm_trace_gas_get_var_index(name, index, exit_on_err)

    !-----------------------------------------------------------------------
    ! Search for a box_atm_trace_gas variable by name and return its index.
    ! If the name is not found then exit_POP is called, unless
    !   the optional argument exit_on_err is set to .false., in which
    !   case index is set to 0.
    !-----------------------------------------------------------------------

    character (len=*), intent(in)            :: name        ! name of variable to be looked up
    integer (int_kind), intent(out)          :: index
    logical (log_kind), optional, intent(in) :: exit_on_err ! Is exit_POP called if name not found?

    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------

    character (len=*), parameter :: subname = 'box_atm_trace_gas_mod:box_atm_trace_gas_get_var_index'
    integer (int_kind) :: i               ! loop index
    logical (log_kind) :: loc_exit_on_err ! local copy of exit_on_err

    !-----------------------------------------------------------------------

    do i = 1, box_atm_trace_gas_cnt
      if (trim(box_atm_trace_gas_array(i)%name) == trim(name)) then
        index = i
        return
      end if
    end do

    if (.not. present(exit_on_err)) then
      loc_exit_on_err = .true.
    else
      loc_exit_on_err = exit_on_err
    end if

    ! name not found
    if (loc_exit_on_err) then
      call document(subname, 'name', trim(name))
      call exit_POP(sigAbort, 'name not found')
    end if

    index = 0

    !-----------------------------------------------------------------------

  end subroutine box_atm_trace_gas_get_var_index

  !-----------------------------------------------------------------------

  function box_atm_trace_gas_var_exists_in_file(index, filename) result(exists)

    !-----------------------------------------------------------------------
    ! Determine if a box_atm_trace_gas variable exists in a file.
    !-----------------------------------------------------------------------

    use passive_tracer_tools, only: field_exists_in_file

    integer (int_kind), intent(in) :: index
    character (len=*), intent(in)  :: filename
    logical (log_kind)             :: exists

    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------

    character (len=*), parameter :: subname = 'box_atm_trace_gas_mod:box_atm_trace_gas_var_exists_in_file'

    !-----------------------------------------------------------------------
    ! error checking
    !-----------------------------------------------------------------------

    if (index < 1 .or. index > box_atm_trace_gas_cnt) then
      call document(subname, 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
    end if

    exists = field_exists_in_file('nc', filename, box_atm_trace_gas_array(index)%file_varname)

    !-----------------------------------------------------------------------

  end function box_atm_trace_gas_var_exists_in_file

  !-----------------------------------------------------------------------

  subroutine box_atm_trace_gas_init_var_filename(index, filename)

    !-----------------------------------------------------------------------
    ! Initialize a box_atm_trace_gas variable from a file.
    !-----------------------------------------------------------------------

    use io_types, only: datafile, io_field_desc
    use io_types, only: construct_file, construct_io_field
    use io_types, only: destroy_file, destroy_io_field
    use io, only: data_set

    integer (int_kind), intent(in) :: index
    character (len=*), intent(in)  :: filename

    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------

    character (len=*), parameter :: subname = 'box_atm_trace_gas_mod:box_atm_trace_gas_init_var_filename'
    type (datafile)      :: file      ! io file descriptor
    type (io_field_desc) :: vals_desc ! io field descriptor

    !-----------------------------------------------------------------------
    ! error checking
    !-----------------------------------------------------------------------

    if (index < 1 .or. index > box_atm_trace_gas_cnt) then
      call document(subname, 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
    end if

    !-----------------------------------------------------------------------
    ! initialize values from file
    !-----------------------------------------------------------------------

    call document(subname, 'name', box_atm_trace_gas_array(index)%name)
    call document(subname, 'file_varname', box_atm_trace_gas_array(index)%file_varname)

    file = construct_file('nc', full_name=trim(filename))
    call data_set(file, 'open_read')
    vals_desc = construct_io_field(box_atm_trace_gas_array(index)%file_varname, &
         d0d_array=box_atm_trace_gas_array(index)%val)
    call data_set(file, 'define', vals_desc)
    call data_set(file, 'read', vals_desc)
    call destroy_io_field(vals_desc)
    call data_set(file, 'close')
    call destroy_file(file)

    box_atm_trace_gas_array(index)%linit = .true.

    !-----------------------------------------------------------------------

  end subroutine box_atm_trace_gas_init_var_filename

  !-----------------------------------------------------------------------

  subroutine box_atm_trace_gas_init_var_val(index, val)

    !-----------------------------------------------------------------------
    ! Initialize a box_atm_trace_gas variable from specified value.
    !-----------------------------------------------------------------------

    integer (int_kind), intent(in) :: index
    real (r8), intent(in)          :: val

    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------

    character (len=*), parameter :: subname = 'box_atm_trace_gas_mod:box_atm_trace_gas_init_var_val'

    !-----------------------------------------------------------------------
    ! error checking
    !-----------------------------------------------------------------------

    if (index < 1 .or. index > box_atm_trace_gas_cnt) then
      call document('box_atm_trace_gas_init_var_val', 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
    end if

    call document(subname, 'name', box_atm_trace_gas_array(index)%name)

    !-----------------------------------------------------------------------
    ! initialize value, converting from user vmr units to internal vmr units
    !-----------------------------------------------------------------------

    box_atm_trace_gas_array(index)%val = val / box_atm_trace_gas_array(index)%vmr_conv_factor

    box_atm_trace_gas_array(index)%linit = .true.

    !-----------------------------------------------------------------------

  end subroutine box_atm_trace_gas_init_var_val

  !-----------------------------------------------------------------------

  subroutine box_atm_trace_gas_update_var(index, flux_vals)

    !-----------------------------------------------------------------------
    ! Update a box_atm_trace_gas variable using specified fluxes.
    !-----------------------------------------------------------------------

    use grid, only: TAREA, RCALCT
    use domain, only: distrb_clinic
    use constants, only: c1, p5, field_loc_center
    use global_reductions, only : global_sum_prod
    use time_management, only: lpre_time_manager, avg_ts_next, back_to_back_next, avg_ts, back_to_back, dtt

    integer (int_kind), intent(in) :: index
    real (r8), intent(in)          :: flux_vals(:,:,:)

    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------

    character (len=*), parameter :: subname = 'box_atm_trace_gas_mod:box_atm_trace_gas_update_var'

    real (r8) :: gsum   ! global integral of flux_vals [mol/s]
    real (r8) :: dt_loc ! timestep [s]

    real (r8), parameter :: atm_mol_w    = 28.97_r8     ! molecular weight of dry air [g/mol]
    real (r8), parameter :: atm_dry_mass = 5.1352e21_r8 ! dry mass of atmsphere [g]
    real (r8), parameter :: atm_moles    = atm_dry_mass / atm_mol_w ! number of dry moles in atmosphere
    real (r8), parameter :: conv_factor  = c1 / atm_moles

    !-----------------------------------------------------------------------
    ! error checking
    !-----------------------------------------------------------------------

    if (index < 1 .or. index > box_atm_trace_gas_cnt) then
      call document(subname, 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
    end if

    if (.not. box_atm_trace_gas_array(index)%linit) then
      call document(subname, 'index', index)
      call document(subname, 'name', box_atm_trace_gas_array(index)%name)
      call exit_POP(sigAbort, 'updating uninitialized box_atm_trace_gas variable')
    end if

    !-----------------------------------------------------------------------
    ! compute global integral of flux_vals, convert to mol/s
    !-----------------------------------------------------------------------

    gsum = global_sum_prod(flux_vals, TAREA, distrb_clinic, field_loc_center, RCALCT)
    gsum = box_atm_trace_gas_array(index)%flux_conv_factor * gsum

    !-----------------------------------------------------------------------
    ! convert integral into change in vmr and increment value
    !-----------------------------------------------------------------------

    if (lpre_time_manager) then
      if (avg_ts_next .or. back_to_back_next) then
        dt_loc = -p5*dtt
      else
        dt_loc = -dtt
      end if
    else
      if (avg_ts .or. back_to_back) then
        dt_loc = -p5*dtt
      else
        dt_loc = -dtt
      end if
    end if

    box_atm_trace_gas_array(index)%val = box_atm_trace_gas_array(index)%val &
         + dt_loc * conv_factor * gsum

    !-----------------------------------------------------------------------

  end subroutine box_atm_trace_gas_update_var

  !-----------------------------------------------------------------------

  function box_atm_trace_gas_get_var_val(index) result(val)

    !-----------------------------------------------------------------------
    ! Get value of a box_atm_trace_gas variable
    !-----------------------------------------------------------------------

    integer (int_kind), intent(in) :: index
    real (r8)                      :: val

    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------

    character (len=*), parameter :: subname = 'box_atm_trace_gas_mod:box_atm_trace_gas_get_var_val'

    !-----------------------------------------------------------------------
    ! error checking
    !-----------------------------------------------------------------------

    if (index < 1 .or. index > box_atm_trace_gas_cnt) then
      call document(subname, 'index', index)
      call exit_POP(sigAbort, 'index out of bounds')
    end if

    if (.not. box_atm_trace_gas_array(index)%linit) then
      call document(subname, 'index', index)
      call document(subname, 'name', box_atm_trace_gas_array(index)%name)
      call exit_POP(sigAbort, 'box_atm_trace_gas variable is not initialized')
    end if

    !-----------------------------------------------------------------------
    ! get value, converting from internal vmr units to user vmr units
    !-----------------------------------------------------------------------

    val = box_atm_trace_gas_array(index)%val * box_atm_trace_gas_array(index)%vmr_conv_factor

    !-----------------------------------------------------------------------

  end function box_atm_trace_gas_get_var_val

  !-----------------------------------------------------------------------

  subroutine box_atm_trace_gas_write_restart(restart_file, action)

    !-----------------------------------------------------------------------
    ! Write necessary data to the restart file
    !-----------------------------------------------------------------------

    use io_types, only: datafile, io_field_desc
    use io_types, only: construct_io_field, destroy_io_field
    use io, only: data_set

    type (datafile), intent(inout) :: restart_file
    character (len=*), intent(in)  :: action

    !-----------------------------------------------------------------------
    ! local variables
    !-----------------------------------------------------------------------

    integer (int_kind) :: index
    type (io_field_desc), allocatable, save :: vals_desc(:) ! io field descriptor

    !-----------------------------------------------------------------------

    if (trim(action) == 'define') then
       allocate(vals_desc(box_atm_trace_gas_cnt))

      do index = 1, box_atm_trace_gas_cnt
        vals_desc(index) = construct_io_field(box_atm_trace_gas_array(index)%file_varname, &
             d0d_array=box_atm_trace_gas_array(index)%val)
        call data_set(restart_file, 'define', vals_desc(index))
      end do
    end if

    if (trim(action) == 'write') then
      do index = 1, box_atm_trace_gas_cnt
        call data_set(restart_file, 'write', vals_desc(index))
        call destroy_io_field(vals_desc(index))
      end do

      deallocate(vals_desc)
    end if

    !-----------------------------------------------------------------------

  end subroutine box_atm_trace_gas_write_restart

  !-----------------------------------------------------------------------

end module box_atm_trace_gas_mod
