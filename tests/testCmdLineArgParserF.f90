program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double

    implicit none
    integer(c_size_t)               :: prsr
    integer(c_int)                  :: ier, nargs, i, n
    character(len=512), allocatable :: argv(:)
    integer(c_int), allocatable     :: arg_lengths(:)
    character(len=512)              :: argv_full
    integer(c_int)                  :: ival
    real(c_double)                  :: dval
    character(len=512)              :: sval

    ier = mnt_cmdlineargparser_new(prsr)

    ! set 
    ier = mnt_cmdlineargparser_setint(prsr, "-i", -123, "some integer")
    ier = mnt_cmdlineargparser_setdouble(prsr, "-d", -1.23_c_double, "some double")
    ier = mnt_cmdlineargparser_setstring(prsr, "-s", "hello_there", "some string")

    ! parse
    nargs = command_argument_count()
    allocate(argv(nargs))
    allocate(arg_lengths(nargs))
    do i = 1, nargs
        call get_command_argument(i, argv_full)
        argv(i) = trim(argv_full)
        arg_lengths(i) = len(argv(i))
    enddo
    ier = mnt_cmdlineargparser_parse(prsr, nargs, arg_lengths, argv)

    ! extract
    ier = mnt_cmdlineargparser_getint(prsr, "-i", ival)
    ier = mnt_cmdlineargparser_getdouble(prsr, "-d", dval)
    ier = mnt_cmdlineargparser_getstring(prsr, "-s", sval, n)

    print *, "-i arg is ", ival
    print *, "-d arg is ", dval
    print *, "-s arg is ", sval

    ! clean up
    deallocate(arg_lengths)
    deallocate(argv)
    ier = mnt_cmdlineargparser_del(prsr)

end program 
