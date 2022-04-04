module strdata_interface_mod

  ! !MODULE: strdata_interface_mod
  !
  ! !DESCRIPTION:
  !
  !  Provide an interface to shr_strdata to use to read input data sets such as
  !  initial data fields for ecosystem tracer fields use CDEPS functionality

  use ESMF
  use dshr_strdata_mod, only : shr_strdata_type
  use dshr_strdata_mod, only : shr_strdata_init_from_inline
  use dshr_strdata_mod, only : shr_strdata_advance
  use dshr_strdata_mod, only : shr_strdata_print
  use dshr_methods_mod, only : dshr_fldbun_getFldPtr

  use kinds_mod,        only : char_len
  use kinds_mod,        only : int_kind
  use kinds_mod,        only : log_kind
  use kinds_mod,        only : r8
  use communicate,      only : my_task
  use communicate,      only : master_task
  use io_types,         only : stdout
  use POP_CommMod,      only : POP_communicator

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
  public :: POP_strdata_get_streamdata
  public :: POP_strdata_set_n0

  interface POP_strdata_get_streamdata
    module procedure POP_strdata_get_streamdata_1d
    module procedure POP_strdata_get_streamdata_2d
  end interface

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
    character(len=char_len), allocatable :: field_names(:)
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

    use shr_string_mod, only: shr_string_listGetNum, shr_string_listGetName
    use ocn_comp_shr  , only: model_clock, model_meshfile, model_mesh

    type(strdata_input_type), intent(inout) :: inputlist

    character(char_len) :: stream_lev_dimname ! level dimension name
    integer(int_kind)   :: num_fields
    integer(int_kind)   :: n
    integer(int_kind)   :: rc

    ! create an array of names from the colon delimited fieldlist
    num_fields = shr_string_listGetNum(inputlist%field_list)
    allocate(inputlist%field_names(num_fields))
    do n = 1,num_fields
       call shr_string_listGetName(inputlist%field_list, n, inputlist%field_names(n))
    end do

    ! determine stream_lev_dimname
    if (inputlist%depth_flag) then
       stream_lev_dimname = 'depth'
    else
       stream_lev_dimname = 'null'
    end if

    call shr_strdata_init_from_inline(                        &
         inputlist%sdat,                                      &
         my_task             = my_task,                       &
         logunit             = stdout,                        &
         compname            = 'OCN',                         &
         model_clock         = model_clock,                   &
         model_mesh          = model_mesh,                    &
         stream_meshfile     = trim(model_meshfile),          &
         stream_lev_dimname  = trim(stream_lev_dimname),      &
         stream_mapalgo      = 'bilinear',                    &
         stream_filenames    = (/trim(inputlist%file_name)/), &
         stream_fldListFile  = inputlist%field_names,         &
         stream_fldListModel = inputlist%field_names,         &
         stream_yearFirst    = inputlist%year_first,          &
         stream_yearLast     = inputlist%year_last,           &
         stream_yearAlign    = inputlist%year_align,          &
         stream_offset       = 0,                             &
         stream_taxMode      = inputlist%taxMode,             &
         stream_dtlimit      = 1.e30_r8,                      &
         stream_tintalgo     = inputlist%tintalgo,            &
         stream_name         = 'unset',                       &
         rc                  = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

    if (my_task == master_task) then
      call shr_strdata_print(inputlist%sdat, 'pop_stream')
    endif

  end subroutine POP_strdata_create

  !***********************************************************************

  subroutine POP_strdata_advance_scalar(inputlist)

    use time_management, only : iyear, imonth, iday
    use time_management, only : ihour, iminute, isecond

    type(strdata_input_type), intent(inout) :: inputlist

    integer(int_kind) :: date
    integer(int_kind) :: time
    integer(int_kind) :: rc

    date = iyear*10000 + imonth*100 + iday
    time = isecond + 60 * (iminute + 60 * ihour)

    call shr_strdata_advance(inputlist%sdat, date, time, stdout, trim(inputlist%timer_label), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if

  end subroutine POP_strdata_advance_scalar

  !***********************************************************************

  subroutine POP_strdata_advance_array(inputlist)

    use time_management, only : iyear, imonth, iday
    use time_management, only : ihour, iminute, isecond

    type(strdata_input_type), intent(inout) :: inputlist(:)

    integer(int_kind) :: date
    integer(int_kind) :: time
    integer(int_kind) :: n
    integer(int_kind) :: rc

    date = iyear*10000 + imonth*100 + iday
    time = isecond + 60 * (iminute + 60 * ihour)

    do n = 1, size(inputlist)
       call shr_strdata_advance(inputlist(n)%sdat, date, time, stdout, trim(inputlist(n)%timer_label), rc=rc)
    end do

  end subroutine POP_strdata_advance_array

  !***********************************************************************

  subroutine POP_strdata_get_streamdata_1d(inputlist, var_index, stream_data1d)
    type(strdata_input_type) , intent(inout) :: inputlist
    integer                  , intent(in)    :: var_index
    real(r8)                 , pointer       :: stream_data1d(:)

    integer           :: lsize
    integer           :: rc
    real(r8), pointer :: dataptr1d(:)

    lsize = inputlist%sdat%model_lsize
    allocate(stream_data1d(lsize))

    call dshr_fldbun_getFldPtr(inputlist%sdat%pstrm(1)%fldbun_model, inputlist%field_names(var_index), &
         fldptr1=dataptr1d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    stream_data1d(:) = dataptr1d(:)

  end subroutine POP_strdata_get_streamdata_1d

  !***********************************************************************

  subroutine POP_strdata_get_streamdata_2d(inputlist, var_index, nlev, stream_data2d)
    type(strdata_input_type) , intent(inout) :: inputlist
    integer                  , intent(in)    :: var_index
    integer                  , intent(in)    :: nlev
    real(r8)                 , pointer       :: stream_data2d(:,:)

    integer           :: rc
    integer           :: n,k
    integer           :: lsize
    real(r8), pointer :: dataptr2d(:,:)

    lsize = inputlist%sdat%model_lsize
    allocate(stream_data2d(nlev,lsize))

    ! Note: in the stream_data output variable the level dimension is the second index
    call dshr_fldbun_getFldPtr(inputlist%sdat%pstrm(1)%fldbun_model, inputlist%field_names(var_index), &
         fldptr2=dataptr2d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) then
       call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end if
    stream_data2d(:,:) = dataptr2d(:,:)

  end subroutine POP_strdata_get_streamdata_2d

  !***********************************************************************

  subroutine POP_strdata_set_n0(nlev, n0)

    use domain, only : nblocks_clinic, blocks_clinic
    use blocks, only : block, get_block

    integer, intent(in)  :: nlev  ! not used here
    integer, intent(out) :: n0(nblocks_clinic)

    integer    :: iblock
    type(block):: this_block  ! block info for the current block

    n0(1) = 0
    do iblock = 1, nblocks_clinic-1
       this_block = get_block(blocks_clinic(iblock), iblock)
       n0(iblock+1) = n0(iblock) + (this_block%je-this_block%jb+1)*(this_block%ie-this_block%ib+1)
    enddo

  end subroutine POP_strdata_set_n0

end module strdata_interface_mod
