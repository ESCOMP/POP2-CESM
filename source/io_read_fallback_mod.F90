!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module io_read_fallback_mod

!BOP
! !MODULE: io_read_fallback_mod
! !DESCRIPTION:
!  Provide utilities for fallback conditions when attempting to read
!  a field from a file that is not present.
!

! !USES:

   use kinds_mod, only: char_len, r8, int_kind, log_kind
   use io_tools,  only: document
   use exit_mod,  only: exit_POP, sigAbort
   use io_types,  only: datafile, io_field_desc
   use constants, only: c1

   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: io_read_fallback_register_field
   public :: io_read_fallback_register_tracer
   public :: io_read_fallback_is_field_registered
   public :: define_field_fallback
   public :: read_field_fallback

!EOP
!BOC
!-----------------------------------------------------------------------
!  module private types and data
!-----------------------------------------------------------------------

   type :: fallback_type
      character (char_len) :: fieldname
      character (char_len) :: fallback_opt ! valid_values = ['const', 'alt_field']

      ! value for fallback_opt == 'const'
      real (r8) :: const_val

      ! fieldname and scaling value for fallback_opt == 'alt_field'
      character (char_len) :: alt_fieldname
      real (r8) :: scalefactor
   end type

   type (fallback_type), dimension(:), allocatable :: &
      fallback_array   ! array of registered fallbacks

   integer (int_kind) :: &
      fallback_cnt = 0 ! number of registered fallbacks

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: io_read_fallback_register_field
! !INTERFACE:

   subroutine io_read_fallback_register_field(fieldname, fallback_opt, const_val, alt_fieldname, scalefactor)

! !DESCRIPTION:
!  Register a fallback option.
!  It is a fatal error to register a fallback for a previously registered fieldname.
!
!  This should only be called once per task.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)           :: fieldname      ! name of field for which the fallback is being registered
   character (*), intent(in)           :: fallback_opt   ! fallback option being registered

   real (r8), optional, intent(in)     :: const_val      ! value for fallback_opt == 'const'

   character (*), optional, intent(in) :: alt_fieldname  ! fieldname value for fallback_opt == 'alt_field'
   real (r8), optional, intent(in)     :: scalefactor    ! scaling value for fallback_opt == 'alt_field'

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character (*), parameter :: subname = 'io_read_fallback_mod:io_read_fallback_register_field'
   character (char_len)     :: message

   type (fallback_type), dimension(:), allocatable :: &
      fallback_array_tmp ! Temporary array to store current fallback array contents. Used when growing the array.

!-----------------------------------------------------------------------
!  check to see if fieldname already has a registered fallback
!-----------------------------------------------------------------------

   if (io_read_fallback_is_field_registered(fieldname)) then
      write(message, '(3A)') subname, ': fallback already registered for ', trim(fieldname)
      call exit_POP(sigAbort, message)
   endif

   call document(subname, 'fieldname', fieldname)
   call document(subname, 'fallback_opt', fallback_opt)

!-----------------------------------------------------------------------
!  increase the length of the array fallback_array, preserving
!  the current contents if fallback_cnt > 0
!-----------------------------------------------------------------------

   if (fallback_cnt > 0) then
      allocate(fallback_array_tmp(fallback_cnt))
      fallback_array_tmp = fallback_array
      deallocate(fallback_array)
   endif

   allocate(fallback_array(fallback_cnt+1))

   if (fallback_cnt > 0) then
      fallback_array(1:fallback_cnt) = fallback_array_tmp
      deallocate(fallback_array_tmp)
   endif

   fallback_cnt = fallback_cnt + 1
   fallback_array(fallback_cnt)%fieldname    = fieldname
   fallback_array(fallback_cnt)%fallback_opt = fallback_opt

   select case(fallback_opt)

   case ('const')

      if (.not. present(const_val)) then
         write(message, '(2A)') subname, ': fallback_opt = const, but const_val not provided'
         call exit_POP(sigAbort, message)
      endif
      fallback_array(fallback_cnt)%const_val = const_val
      call document(subname, 'const_val', fallback_array(fallback_cnt)%const_val)

   case ('alt_field')

      if (.not. present(alt_fieldname)) then
         write(message, '(2A)') subname, ': fallback_opt = alt_field, but alt_fieldname not provided'
         call exit_POP(sigAbort, message)
      endif
      fallback_array(fallback_cnt)%alt_fieldname = alt_fieldname
      call document(subname, 'alt_fieldname', fallback_array(fallback_cnt)%alt_fieldname)

      if (present(scalefactor)) then
         fallback_array(fallback_cnt)%scalefactor = scalefactor
      else
         fallback_array(fallback_cnt)%scalefactor = c1
      endif
      call document(subname, 'scalefactor', fallback_array(fallback_cnt)%scalefactor)

   case default

      write(message, '(3A)') subname, ': unknown fallback_opt value, ', trim(fallback_opt)
      call exit_POP(sigAbort, message)

   end select

!-----------------------------------------------------------------------
!EOC

   end subroutine io_read_fallback_register_field

!***********************************************************************
!BOP
! !IROUTINE: io_read_fallback_register_tracer
! !INTERFACE:

   subroutine io_read_fallback_register_tracer(tracername, fallback_opt, const_val, alt_tracername, scalefactor)

! !DESCRIPTION:
!  Register a fallback option for a tracer.
!  Call io_read_fallback_register_field for the tracer, and
!  the generated names written to the restart file.
!
!  This should only be called once per task.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in)           :: tracername      ! name of tracer for which the fallback is being registered
   character (*), intent(in)           :: fallback_opt    ! fallback option being registered

   real (r8), optional, intent(in)     :: const_val       ! value for fallback_opt == 'const'

   character (*), optional, intent(in) :: alt_tracername  ! tracername value for fallback_opt == 'alt_field'
   real (r8), optional, intent(in)     :: scalefactor     ! scaling value for fallback_opt == 'alt_field'

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character (*), parameter :: subname = 'io_read_fallback_mod:io_read_fallback_register_tracer'
   character (char_len)     :: message

!-----------------------------------------------------------------------

   call document(subname, 'tracername', tracername)
   call document(subname, 'fallback_opt', fallback_opt)

   select case(fallback_opt)

   case ('const')

      if (.not. present(const_val)) then
         write(message, '(2A)') subname, ': fallback_opt = const, but const_val not provided'
         call exit_POP(sigAbort, message)
      endif
      call io_read_fallback_register_field(tracername, fallback_opt, const_val)
      call io_read_fallback_register_field(tracername//'_CUR', fallback_opt, const_val)
      call io_read_fallback_register_field(tracername//'_OLD', fallback_opt, const_val)

   case ('alt_field')

      if (.not. present(alt_tracername)) then
         write(message, '(2A)') subname, ': fallback_opt = alt_field, but alt_tracername not provided'
         call exit_POP(sigAbort, message)
      endif
      call io_read_fallback_register_field(tracername, fallback_opt, &
           alt_fieldname=alt_tracername, scalefactor=scalefactor)
      call io_read_fallback_register_field(tracername//'_CUR', fallback_opt, &
           alt_fieldname=alt_tracername//'_CUR', scalefactor=scalefactor)
      call io_read_fallback_register_field(tracername//'_OLD', fallback_opt, &
           alt_fieldname=alt_tracername//'_OLD', scalefactor=scalefactor)

   case default

      write(message, '(3A)') subname, ': unknown fallback_opt value, ', trim(fallback_opt)
      call exit_POP(sigAbort, message)

   end select

!-----------------------------------------------------------------------
!EOC

   end subroutine io_read_fallback_register_tracer

!***********************************************************************
!BOP
! !IROUTINE: io_read_fallback_is_field_registered
! !INTERFACE:

   function io_read_fallback_is_field_registered(fieldname)

! !DESCRIPTION:
!  Determine if a fallback is registered for fieldname
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: fieldname ! name of field for which fallback registery is being queried

! !OUTPUT PARAMETERS:

   logical (log_kind) :: io_read_fallback_is_field_registered ! result of this function

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: n ! index into array of registered fallbacks

!-----------------------------------------------------------------------

   io_read_fallback_is_field_registered = .false.

   do n = 1, fallback_cnt
      if (fallback_array(n)%fieldname == fieldname) then
         io_read_fallback_is_field_registered = .true.
         exit
      endif
   enddo

!-----------------------------------------------------------------------
!EOC

   end function io_read_fallback_is_field_registered

!***********************************************************************
!BOP
! !IROUTINE: define_field_fallback
! !INTERFACE:

   subroutine define_field_fallback(data_file, io_field)

! !DESCRIPTION:
!  Implement a fallback option to define a field.
!
! !REVISION HISTORY:
!  same as module

   use io_netcdf, only: define_field_netcdf
   use io_netcdf, only: field_exists_netcdf

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)      :: data_file ! file for reads (used for fallback_opt == 'alt_field')

   type (io_field_desc), intent (inout) :: io_field ! field descriptor for field being populated

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character (*), parameter :: subname = 'io_read_fallback_mod:define_field_fallback'
   character (char_len)     :: message
   integer (int_kind)       :: n       ! index into array of registered fallbacks
   logical (log_kind)       :: field_exists

!-----------------------------------------------------------------------
!  get index of registered fallback
!  ensure that field has a registered fallback
!-----------------------------------------------------------------------

   do n = 1, fallback_cnt
      if (fallback_array(n)%fieldname == io_field%short_name) exit
   enddo

   if (n > fallback_cnt) then
      write(message, '(3A)') subname, ': fallback not registered for ', trim(io_field%short_name)
      call exit_POP(sigAbort, message)
   endif

!-----------------------------------------------------------------------
!  implement define for registered fallback
!-----------------------------------------------------------------------

   select case(fallback_array(n)%fallback_opt)

!  case ('const') ! do nothing for fallback_opt = 'const'

   case ('alt_field')

      ! ensure that alt_fieldname is present in data_file
      call field_exists_netcdf(data_file,fallback_array(n)%alt_fieldname,field_exists)

      if (.not. field_exists) then
         write(message, '(6A)') subname, ': alt_fieldname ', trim(fallback_array(n)%alt_fieldname), &
            ' for fieldname ', trim(fallback_array(n)%fieldname), ' does not exist'
         call exit_POP(sigAbort, message)
      endif

      ! temporarily overwrite io_field%short_name for defining
      io_field%short_name = fallback_array(n)%alt_fieldname

      call define_field_netcdf(data_file, io_field)

      ! reset io_field%short_name after defining
      io_field%short_name = fallback_array(n)%fieldname

   end select

!-----------------------------------------------------------------------
!EOC

   end subroutine define_field_fallback

!***********************************************************************
!BOP
! !IROUTINE: read_field_fallback
! !INTERFACE:

   subroutine read_field_fallback(data_file, io_field)

! !DESCRIPTION:
!  Implement a fallback option to read set a field.
!
! !REVISION HISTORY:
!  same as module

   use io_netcdf, only: read_field_netcdf

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)      :: data_file ! file for reads (used for fallback_opt == 'alt_field')

   type (io_field_desc), intent (inout) :: io_field ! field descriptor for field being populated

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character (*), parameter :: subname = 'io_read_fallback_mod:read_field_fallback'
   character (char_len)     :: message
   integer (int_kind)       :: n       ! index into array of registered fallbacks

!-----------------------------------------------------------------------
!  get index of registered fallback
!  ensure that field has a registered fallback
!-----------------------------------------------------------------------

   do n = 1, fallback_cnt
      if (fallback_array(n)%fieldname == io_field%short_name) exit
   enddo

   if (n > fallback_cnt) then
      write(message, '(3A)') subname, ': fallback not registered for ', trim(io_field%short_name)
      call exit_POP(sigAbort, message)
   endif

!-----------------------------------------------------------------------
!  implement read for registered fallback
!-----------------------------------------------------------------------

   call document(subname, 'fieldname', fallback_array(n)%fieldname)
   call document(subname, 'fallback_opt', fallback_array(n)%fallback_opt)

   select case(fallback_array(n)%fallback_opt)

   case ('const')

      call document(subname, 'const_val', fallback_array(n)%const_val)

      if (associated(io_field%field_d_0d)) then
         io_field%field_d_0d          = fallback_array(n)%const_val
      else if (associated(io_field%field_d_1d)) then
         io_field%field_d_1d(:)       = fallback_array(n)%const_val
      else if (associated(io_field%field_d_2d)) then
         io_field%field_d_2d(:,:,:)   = fallback_array(n)%const_val
      else if (associated(io_field%field_d_3d)) then
         io_field%field_d_3d(:,:,:,:) = fallback_array(n)%const_val
      else
         write(message, '(2A)') subname, ': const fallback only supported for r8 fields'
         call exit_POP(sigAbort, message)
      endif

   case ('alt_field')

      call document(subname, 'alt_fieldname', fallback_array(n)%alt_fieldname)
      call document(subname, 'scalefactor', fallback_array(n)%scalefactor)

      ! temporarily overwrite io_field%short_name for reading
      io_field%short_name = fallback_array(n)%alt_fieldname

      call read_field_netcdf(data_file, io_field)

      ! reset io_field%short_name after reading
      io_field%short_name = fallback_array(n)%fieldname

      ! multiply read in value(s) by scalefactor, if it is not 1.0
      if (fallback_array(n)%scalefactor /= c1) then
         if (associated(io_field%field_d_0d)) then
            io_field%field_d_0d          = fallback_array(n)%scalefactor * io_field%field_d_0d
         else if (associated(io_field%field_d_1d)) then
            io_field%field_d_1d(:)       = fallback_array(n)%scalefactor * io_field%field_d_1d(:)
         else if (associated(io_field%field_d_2d)) then
            io_field%field_d_2d(:,:,:)   = fallback_array(n)%scalefactor * io_field%field_d_2d(:,:,:)
         else if (associated(io_field%field_d_3d)) then
            io_field%field_d_3d(:,:,:,:) = fallback_array(n)%scalefactor * io_field%field_d_3d(:,:,:,:)
         else
            write(message, '(2A)') subname, ': alt_field fallback with scaling only supported for r8 fields'
            call exit_POP(sigAbort, message)
         endif
      endif

   end select

!-----------------------------------------------------------------------
!EOC

   end subroutine read_field_fallback

 end module io_read_fallback_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
