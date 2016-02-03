module marbl_logging

  use marbl_kinds_mod, only : char_len
  implicit none
  private
  save

  ! MNL TO-DO: do we want a wrapper class that contains a linked list of
  !            marbl_log_type? Then labort_marbl can be in the wrapped class
  !            (and so can the status level definitions); would make this more
  !            in line with the diagnostic and forcing field types, as well
  type, public :: marbl_status_log_entry_type
    integer :: StatusLevel
    character(len=char_len) :: LogMessage,   &! Message text
                               ModName,      &! Module writing the log
                               SubName        ! Subroutine writing the log
    type(marbl_status_log_entry_type), pointer :: next
  end type marbl_status_log_entry_type

  type :: marbl_log_status_level_type
    integer :: VerboseCode
    integer :: NameListCode
    integer :: DiagnosticCode
    integer :: WarningCode
    integer :: ErrorCode
  contains
    procedure :: construct => marbl_status_level_constructor
  end type marbl_log_status_level_type

  type, public :: marbl_log_type
    logical :: labort_marbl ! True => driver should abort GCM
    type(marbl_log_status_level_type) :: marbl_status_levels
    type(marbl_status_log_entry_type), pointer :: FullLog
    type(marbl_status_log_entry_type), pointer :: LastEntry
  contains
    procedure :: construct => marbl_log_constructor
  end type marbl_log_type


  public :: marbl_log_namelist

contains

  subroutine marbl_status_level_constructor(this, VC, NLC, DC, WC, EC)

    class(marbl_log_status_level_type), intent(inout) :: this
    integer, intent(in) :: VC, NLC, DC, WC, EC

    this%VerboseCode = VC
    this%NamelistCode = NLC
    this%DiagnosticCode = DC
    this%WarningCode = WC
    this%ErrorCode = EC

  end subroutine marbl_status_level_constructor

  subroutine marbl_log_constructor(this)

    class(marbl_log_type), intent(inout) :: this

    this%labort_marbl = .false.
    nullify(this%FullLog)
    nullify(this%LastEntry)
    call this%marbl_status_levels%construct(1,2,3,4,5)

  end subroutine marbl_log_constructor

  subroutine marbl_log_namelist()

  end subroutine marbl_log_namelist

end module marbl_logging
