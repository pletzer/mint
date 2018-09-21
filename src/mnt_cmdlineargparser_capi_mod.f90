 module mnt_cmdlineargparser_capi_mod

  integer, parameter :: mnt_string_size = 1024 ! May need to adjust!

  ! C function prototypes
  interface

    function mnt_cmdlineargparser_new(obj) bind(C)
      ! Constructor
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_cmdlineargparser_new
    end function mnt_cmdlineargparser_new

    function mnt_cmdlineargparser_del(obj) bind(C)
      ! Destructor
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_cmdlineargparser_del
    end function mnt_cmdlineargparser_del

    function mnt_cmdlineargparser_setdouble(obj, name, def_value, help) bind(C)
      ! Set double option
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @param def_value default value
      ! @param help description of the option
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(len=1), intent(in)             :: name
      real(c_double), value                    :: def_value
      character(len=1), intent(in)             :: help
      integer(c_int)                           :: mnt_cmdlineargparser_setdouble
    end function mnt_cmdlineargparser_setdouble

    function mnt_cmdlineargparser_setint(obj, name, def_value, help) bind(C)
      ! Set integer option
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @param def_value default value
      ! @param help description of the option
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(len=1), intent(in)             :: name
      integer(c_int), value                    :: def_value
      character(len=1), intent(in)             :: help
      integer(c_int)                           :: mnt_cmdlineargparser_setint
    end function mnt_cmdlineargparser_setint

    function mnt_cmdlineargparser_setstring(obj, name, def_value, help) bind(C)
      ! Set string option
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @param def_value default value
      ! @param help description of the option
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(len=1), intent(in)             :: name
      character(len=1), intent(in)             :: def_value
      character(len=1), intent(in)             :: help
      integer(c_int)                           :: mnt_cmdlineargparser_setstring
    end function mnt_cmdlineargparser_setstring

    function mnt_cmdlineargparser_setbool(obj, name, def_value, help) bind(C)
      ! Set boolean option
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @param def_value default value as an integer, 0 = false, 1 = true
      ! @param help description of the option
      ! @return 0 if successful
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void*
      character(len=1), intent(in)             :: name
      integer(c_int), value                    :: def_value
      character(len=1), intent(in)             :: help
      integer(c_int)                           :: mnt_cmdlineargparser_setbool
    end function mnt_cmdlineargparser_setbool

    function mnt_cmdlineargparser_parse(obj, nargs, n, args) bind(C)
      ! Parse command line options
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param nargs number of arguments
      ! @param n length of each command line argument string
      ! @param args list of command line options, each option must be \0 terminated
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_int), value                    :: nargs
      integer(c_int), value                    :: n
      character(len=1), intent(in)             :: args(*)
      integer(c_int)                           :: mnt_cmdlineargparser_parse
    end function mnt_cmdlineargparser_parse

    function mnt_cmdlineargparser_help(obj) bind(C)
      ! Print help message
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_int)                           :: mnt_cmdlineargparser_help
    end function mnt_cmdlineargparser_help

    function mnt_cmdlineargparser_getdouble(obj, name, value) bind(C)
      ! Get double command line argument value
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @param value return value
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void*
      character(len=1), intent(in)             :: name
      real(c_double), intent(out)              :: value
      integer(c_int)                           :: mnt_cmdlineargparser_getdouble
    end function mnt_cmdlineargparser_getdouble

    function mnt_cmdlineargparser_getint(obj, name, value) bind(C)
      ! Get integer command line argument value
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @param value return value
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(len=1), intent(in)             :: name
      integer(c_int), intent(out)              :: value
      integer(c_int)                           :: mnt_cmdlineargparser_getint
    end function mnt_cmdlineargparser_getint

    function mnt_cmdlineargparser_getstring(obj, name, value, n) bind(C)
      ! Get string command line argument value
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @param value return value
      ! @param n returned length of string
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(len=1), intent(in)             :: name(*)
      character(len=1), intent(out)            :: value(*)
      integer(c_int), intent(out)              :: n
      integer(c_int)                           :: mnt_cmdlineargparser_getstring
    end function mnt_cmdlineargparser_getstring

    function mnt_cmdlineargparser_getbool(obj, name, value) bind(C)
      ! Get boolean command line argument value
      ! @param obj instance of mntCmdLineArgParser_t (opaque handle)
      ! @param name name of the command line option (must be \0 terminated)
      ! @param value return value as an integer, 0 = false, 1 = true
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(len=1), intent(in)             :: name
      integer(c_int), intent(out)              :: value
      integer(c_int)                           :: mnt_cmdlineargparser_getbool
    end function mnt_cmdlineargparser_getbool

  end interface

contains

  subroutine mnt_f2c_string(f_string, c_string)
    ! Convert a Fortran to a C string
    ! @param f_string input Fortran string 
    ! @param c_string output C string, null terminated, trimmed and left-adjusted
    implicit none
    character(len=*), intent(in)  :: f_string
    character(len=1), intent(out) :: c_string(:)
    integer                       :: i, n
    n = len_trim(f_string)
    c_string(:) = char(0)
    do i = 1, n
      c_string(i) = f_string(i:i)
    enddo
  end subroutine mnt_f2c_string

  subroutine mnt_c2f_string(c_string, f_string)
    ! Convert a C to a Fortran string, detecting \0 (char(0)) 
    ! and setting any characted following \0 with a space
    ! @param c_string input C string 
    ! @param f_string output Fortran string, replacing 
    implicit none
    character(len=1), intent(in)  :: c_string(:)
    character(len=*), intent(out) :: f_string

    ! local variables
    integer                       :: i, n_f_string, n
    logical                       :: set_to_empty

    n_f_string = len(f_string)
    n = min(n_f_string, size(c_string))
    set_to_empty = .FALSE.
    ! replace '\0' with ' '
    do i = 1, n
      
      if (set_to_empty .OR. c_string(i) == char(0)) then
        ! end of C string detected
        f_string(i:i) = ' '
        set_to_empty = .TRUE.
      else
        f_string(i:i) = c_string(i)
      endif
    enddo

  end subroutine mnt_c2f_string

end module mnt_cmdlineargparser_capi_mod

