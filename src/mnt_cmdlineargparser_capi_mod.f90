module mnt_cmdlineargparser_capi_mod

  integer, parameter :: mnt_string_size = 64 ! May need to adjust!

  ! C function prototypes
  interface

    function mnt_cmdlineargparser_new(obj) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout) :: obj
      integer(c_int)                   :: mnt_cmdlineargparser_new
    end function mnt_cmdlineargparser_new

    function mnt_cmdlineargparser_del(obj) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout) :: obj
      integer(c_int)                   :: mnt_cmdlineargparser_del
    end function mnt_cmdlineargparser_del

    function mnt_cmdlineargparser_setdouble(obj, name, def_value, help) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: name
      real(c_double), value                    :: def_value
      character(len=1), intent(in)             :: help
      integer(c_int)                           :: mnt_cmdlineargparser_setdouble
    end function mnt_cmdlineargparser_setdouble

    function mnt_cmdlineargparser_setint(obj, name, def_value, help) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: name
      integer(c_int), value                    :: def_value
      character(len=1), intent(in)             :: help
      integer(c_int)                           :: mnt_cmdlineargparser_setint
    end function mnt_cmdlineargparser_setint

    function mnt_cmdlineargparser_setstring(obj, name, def_value, help) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: name
      character(len=1), intent(in)             :: def_value
      character(len=1), intent(in)             :: help
      integer(c_int)                           :: mnt_cmdlineargparser_setstring
    end function mnt_cmdlineargparser_setstring

    function mnt_cmdlineargparser_setbool(obj, name, def_value, help) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: name
      integer(c_int), value                    :: def_value
      character(len=1), intent(in)             :: help
      integer(c_int)                           :: mnt_cmdlineargparser_setint
    end function mnt_cmdlineargparser_setbool

    function mnt_cmdlineargparser_parse(obj, nargs, n, args) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      integer(c_int), value                    :: nargs
      integer(c_int), value                    :: n
      character(len=1), intent(in)             :: args(*)
      integer(c_int)                           :: mnt_cmdlineargparser_parse
    end function mnt_cmdlineargparser_parse

    function mnt_cmdlineargparser_help(obj) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      integer(c_int)                           :: mnt_cmdlineargparser_parse
    end function mnt_cmdlineargparser_help

    function mnt_cmdlineargparser_getdouble(obj, name, value) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: name
      real(c_double), intent(out)              :: value
      integer(c_int)                           :: mnt_cmdlineargparser_getdouble
    end function mnt_cmdlineargparser_getdouble

    function mnt_cmdlineargparser_getint(obj, name, value) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: name
      integer(c_int), intent(out)              :: value
      integer(c_int)                           :: mnt_cmdlineargparser_getint
    end function mnt_cmdlineargparser_getint

    function mnt_cmdlineargparser_getstring(obj, name, value, n) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: name
      character(len=1), intent(out)            :: value
      integer(c_int), intent(out)              :: n
      integer(c_int)                           :: mnt_cmdlineargparser_getstring
    end function mnt_cmdlineargparser_getstring

    function mnt_cmdlineargparser_getbool(obj, name, value) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: name
      integer(c_int), intent(out)              :: value
      integer(c_int)                           :: mnt_cmdlineargparser_getint
    end function mnt_cmdlineargparser_getbool

  end interface

contains

  subroutine mnt_make_c_string(fort_string, c_string)
    implicit none
    character(len=*), intent(in)  :: fort_string
    character(len=*), intent(out) :: c_string
    c_string = trim(adjustl(fort_string))//char(0)
  end subroutine mnt_make_c_string


end module mnt_cmdlineargparser_capi_mod

