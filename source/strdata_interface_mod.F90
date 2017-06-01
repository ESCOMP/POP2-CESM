module strdata_interface_mod

!BOP
! !MODULE: strdata_interface_mod
!
! !DESCRIPTION:
!
!  Provide an interface to shr_strdata to use to read input data sets such as
!  initial data fields for ecosystem tracer fields.
!
!  Initial work done by mlevy (July 2015)

  use shr_strdata_mod, only : shr_strdata_type
  use shr_strdata_mod, only : shr_strdata_create
  use shr_strdata_mod, only : shr_strdata_print
  use shr_strdata_mod, only : shr_strdata_advance
  use shr_pio_mod,     only : shr_pio_getiotype
  use shr_pio_mod,     only : shr_pio_getiosys

  use kinds_mod,        only : char_len
  use kinds_mod,        only : int_kind
  use kinds_mod,        only : log_kind
  use domain_size,      only : nx_global
  use domain_size,      only : ny_global
  use domain_size,      only : km
  use communicate,      only : my_task
  use communicate,      only : master_task
  use POP_IOUnitsMod,   only : inst_name
  use POP_CommMod,      only : POP_communicator
  use POP_MCT_vars_mod, only : POP_MCT_OCNID
  use POP_MCT_vars_mod, only : POP_MCT_gsMap_o
  use POP_MCT_vars_mod, only : POP_MCT_gsMap3d_o
  use POP_MCT_vars_mod, only : POP_MCT_dom_o
  use POP_MCT_vars_mod, only : POP_MCT_dom3d_o

  implicit none
  private

  public :: strdata_input_type
  public :: POP_strdata_type_set
  public :: POP_strdata_type_match
  public :: POP_strdata_type_cp
  public :: POP_strdata_type_append_field
  public :: POP_strdata_type_field_count
  public :: POP_strdata_create
  public :: POP_strdata_advance

  interface POP_strdata_advance
    module procedure POP_strdata_advance_scalar
    module procedure POP_strdata_advance_array
  end interface

  type strdata_input_type
    ! Contains arguments for shr_strdata_* that are unique to variables being read
    ! all type components, except for sdat, set by POP_strdata_type_set
    ! sdat set by POP_strdata_create
    character(len=char_len) :: file_name
    character(len=char_len) :: field_list
    character(len=char_len) :: timer_label
    integer(int_kind)       :: year_first
    integer(int_kind)       :: year_last
    integer(int_kind)       :: year_align
    logical(log_kind)       :: depth_flag
    character(len=char_len) :: tintalgo
    character(len=char_len) :: taxMode
    type(shr_strdata_type)  :: sdat
  end type strdata_input_type

contains

  !***********************************************************************

  subroutine POP_strdata_type_set(strdata_input_var, &
       file_name, field, timer_label, year_first, year_last, year_align, &
       depth_flag, tintalgo, taxMode)

    type (strdata_input_type), intent(out)   :: strdata_input_var
    character (len=*),  intent(in)           :: file_name
    character (len=*),  intent(in)           :: field
    character (len=*),  intent(in)           :: timer_label
    integer (int_kind), intent(in)           :: year_first
    integer (int_kind), intent(in)           :: year_last
    integer (int_kind), intent(in)           :: year_align
    logical (log_kind), optional, intent(in) :: depth_flag
    character (len=*),  optional, intent(in) :: tintalgo
    character (len=*),  optional, intent(in) :: taxMode

    strdata_input_var%file_name    = trim(file_name)
    strdata_input_var%field_list   = trim(field)
    strdata_input_var%timer_label  = trim(timer_label)
    strdata_input_var%year_first   = year_first
    strdata_input_var%year_last    = year_last
    strdata_input_var%year_align   = year_align

    strdata_input_var%depth_flag   = .false.
    if (present(depth_flag)) then
      strdata_input_var%depth_flag = depth_flag
    endif

    strdata_input_var%tintalgo     = 'linear'
    if (present(tintalgo)) then
      strdata_input_var%tintalgo   = tintalgo
    endif

    strdata_input_var%taxMode      = 'cycle'
    if (present(taxMode)) then
      strdata_input_var%taxMode    = taxMode
    endif

  end subroutine POP_strdata_type_set

  !***********************************************************************

  function POP_strdata_type_match(strdata_input_var1, strdata_input_var2)

    type (strdata_input_type), intent(in) :: strdata_input_var1
    type (strdata_input_type), intent(in) :: strdata_input_var2

    logical (log_kind) :: POP_strdata_type_match

    POP_strdata_type_match = &
      strdata_input_var1%file_name     == strdata_input_var2%file_name  .and. &
      strdata_input_var1%year_first    == strdata_input_var2%year_first .and. &
      strdata_input_var1%year_last     == strdata_input_var2%year_last  .and. &
      strdata_input_var1%year_align    == strdata_input_var2%year_align .and. &
      strdata_input_var1%depth_flag .eqv. strdata_input_var2%depth_flag .and. &
      strdata_input_var1%tintalgo      == strdata_input_var2%tintalgo   .and. &
      strdata_input_var1%taxMode       == strdata_input_var2%taxMode

  end function POP_strdata_type_match

  !***********************************************************************

  subroutine POP_strdata_type_cp(strdata_input_var_in, strdata_input_var_out)

    type (strdata_input_type), intent(in)  :: strdata_input_var_in
    type (strdata_input_type), intent(out) :: strdata_input_var_out

    strdata_input_var_out%file_name   = strdata_input_var_in%file_name
    strdata_input_var_out%field_list  = strdata_input_var_in%field_list
    strdata_input_var_out%timer_label = strdata_input_var_in%timer_label
    strdata_input_var_out%year_first  = strdata_input_var_in%year_first
    strdata_input_var_out%year_last   = strdata_input_var_in%year_last
    strdata_input_var_out%year_align  = strdata_input_var_in%year_align
    strdata_input_var_out%depth_flag  = strdata_input_var_in%depth_flag
    strdata_input_var_out%tintalgo    = strdata_input_var_in%tintalgo
    strdata_input_var_out%taxMode     = strdata_input_var_in%taxMode

  end subroutine POP_strdata_type_cp

  !***********************************************************************

  subroutine POP_strdata_type_append_field(field, strdata_input_var)

    type (strdata_input_type), intent(inout) :: strdata_input_var
    character (len=*),         intent(in)    :: field

    strdata_input_var%field_list = trim(strdata_input_var%field_list) // ':' // trim(field)

  end subroutine POP_strdata_type_append_field

  !***********************************************************************

  function POP_strdata_type_field_count(strdata_input_var)

    use shr_string_mod, only : shr_string_countChar

    type (strdata_input_type), intent(in) :: strdata_input_var

    integer (int_kind) :: POP_strdata_type_field_count

    POP_strdata_type_field_count = shr_string_countChar(strdata_input_var%field_list, ':') + 1

  end function POP_strdata_type_field_count

  !***********************************************************************

  subroutine POP_strdata_create(inputlist)

    type(strdata_input_type), intent(inout) :: inputlist

    if (inputlist%depth_flag) then
       !--- include nzg and domZvarName in call
       call shr_strdata_create(inputlist%sdat, name='not_used',               &
                            mpicom=POP_communicator,                          &
                            compid=POP_MCT_OCNID,                             &
                            gsmap=POP_MCT_gsMap3d_o, ggrid=POP_MCT_dom3d_o,   &
                            nxg=nx_global, nyg=ny_global, nzg=km,             &
                            yearFirst=inputlist%year_first,                   &
                            yearLast=inputlist%year_last,                     &
                            yearAlign=inputlist%year_align,                   &
                            offset=0,                                         &
                            domFilePath='',                                   &
                            domFileName=inputlist%file_name,                  &
                            domTvarName='time',                               &
                            domXvarName='TLONG', domYvarName='TLAT',          &
                            domZvarName='depth',                              &
                            domAreaName='TAREA', domMaskName='KMT',           &
                            FilePath='',                                      &
                            FileName=(/trim(inputlist%file_name)/),           &
                            fldListFile=inputlist%field_list,                 &
                            fldListModel=inputlist%field_list,                &
                            pio_subsystem=shr_pio_getiosys(inst_name),        &
                            pio_iotype=shr_pio_getiotype(inst_name),          &
                            fillalgo='none', mapalgo='none',                  &
                            tintalgo=inputlist%tintalgo,                      &
                            taxMode=inputlist%taxMode)
    else
       call shr_strdata_create(inputlist%sdat, name='not_used',               &
                            mpicom=POP_communicator,                          &
                            compid=POP_MCT_OCNID,                             &
                            gsmap=POP_MCT_gsMap_o, ggrid=POP_MCT_dom_o,       &
                            nxg=nx_global, nyg=ny_global,                     &
                            yearFirst=inputlist%year_first,                   &
                            yearLast=inputlist%year_last,                     &
                            yearAlign=inputlist%year_align,                   &
                            offset=0,                                         &
                            domFilePath='',                                   &
                            domFileName=inputlist%file_name,                  &
                            domTvarName='time',                               &
                            domXvarName='TLONG', domYvarName='TLAT',          &
                            domAreaName='TAREA', domMaskName='KMT',           &
                            FilePath='',                                      &
                            FileName=(/trim(inputlist%file_name)/),           &
                            fldListFile=inputlist%field_list,                 &
                            fldListModel=inputlist%field_list,                &
                            pio_subsystem=shr_pio_getiosys(inst_name),        &
                            pio_iotype=shr_pio_getiotype(inst_name),          &
                            fillalgo='none', mapalgo='none',                  &
                            tintalgo=inputlist%tintalgo,                      &
                            taxMode=inputlist%taxMode)
    endif

    if (my_task == master_task) then
      call shr_strdata_print(inputlist%sdat)
    endif

  end subroutine POP_strdata_create

  !***********************************************************************

  subroutine POP_strdata_advance_scalar(inputlist)

    use time_management, only : iyear, imonth, iday
    use time_management, only : ihour, iminute, isecond

    type(strdata_input_type), intent(inout) :: inputlist

    integer(int_kind) :: date
    integer(int_kind) :: time

    date = iyear*10000 + imonth*100 + iday
    time = isecond + 60 * (iminute + 60 * ihour)

    call shr_strdata_advance(inputlist%sdat,   &
                             date, time,       &
                             POP_communicator, &
                             trim(inputlist%timer_label))

  end subroutine POP_strdata_advance_scalar

  !***********************************************************************

  subroutine POP_strdata_advance_array(inputlist)

    use time_management, only : iyear, imonth, iday
    use time_management, only : ihour, iminute, isecond

    type(strdata_input_type), intent(inout) :: inputlist(:)

    integer(int_kind) :: date
    integer(int_kind) :: time
    integer(int_kind) :: n

    date = iyear*10000 + imonth*100 + iday
    time = isecond + 60 * (iminute + 60 * ihour)

    do n = 1, size(inputlist)
      call shr_strdata_advance(inputlist(n)%sdat, &
                               date, time,        &
                               POP_communicator,  &
                               trim(inputlist(n)%timer_label))
    end do

  end subroutine POP_strdata_advance_array

end module strdata_interface_mod
