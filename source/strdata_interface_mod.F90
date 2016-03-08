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
  use domain_size,      only : nx_global
  use domain_size,      only : ny_global
  use communicate,      only : my_task
  use communicate,      only : master_task
  use POP_IOUnitsMod,   only : inst_name
  use POP_CommMod,      only : POP_communicator
  use POP_MCT_vars_mod, only : POP_MCT_OCNID
  use POP_MCT_vars_mod, only : POP_MCT_gsMap_o
  use POP_MCT_vars_mod, only : POP_MCT_dom_o

  implicit none
  private

  public :: strdata_input_type
  public :: POP_strdata_create
  public :: POP_strdata_advance

  type strdata_input_type
    ! Contains arguments for shr_strdata_* that are unique to variable being
    ! read
    type(shr_strdata_type)  :: sdat
    character(len=char_len) :: timer_label
    integer(int_kind)       :: year_first
    integer(int_kind)       :: year_last
    integer(int_kind)       :: year_align
    character(len=char_len) :: file_name
    character(len=char_len) :: field_list
    integer(int_kind)       :: date
    integer(int_kind)       :: time
  end type strdata_input_type

contains

  subroutine POP_strdata_create(inputlist)

    type(strdata_input_type), intent(inout) :: inputlist

    call shr_strdata_create(inputlist%sdat, name='not_used',                  &
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
                            fillalgo='none', mapalgo='none')

    if (my_task == master_task) then
      call shr_strdata_print(inputlist%sdat)
    endif

  end subroutine POP_strdata_create

  subroutine POP_strdata_advance(inputlist)

    type(strdata_input_type), intent(inout) :: inputlist

    call shr_strdata_advance(inputlist%sdat,                                  &
                             inputlist%date,                                  &
                             inputlist%time,                                  &
                             POP_communicator,                                &
                             trim(inputlist%timer_label))

  end subroutine POP_strdata_advance

end module strdata_interface_mod
