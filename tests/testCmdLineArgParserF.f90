program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod
    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_char

    implicit none
    integer(c_size_t)                                :: prsr
    integer(c_int)                                   :: ier, nargs, i, n, verbosity, nargs1
    character(len=1), allocatable                    :: args(:)
    character(len=mnt_string_size)                   :: argv_full
    integer(c_int)                                   :: ival
    real(c_double)                                   :: dval
    character(len=mnt_string_size)                   :: sval
    character(len=mnt_string_size)                   :: sval_f


    nargs = command_argument_count()
    nargs1 = nargs + 1

    ! args must be a contiguous string
    allocate(args(nargs1 * mnt_string_size))

    do i = 0, nargs 
        ! i = 0 is the executable
        call get_command_argument(i, argv_full)
        ! add termination character, trim...
        call mnt_f2c_string(argv_full, args(i*mnt_string_size + 1))
    enddo

    ! define the command line arguments
    ier = mnt_cmdlineargparser_new(prsr)

    ! set
    ier = mnt_cmdlineargparser_setint(prsr, "-i"//char(0), -123, &
                                      "some integer"//char(0))

    ier = mnt_cmdlineargparser_setdouble(prsr, "-d"//char(0), -1.23_c_double, &
                                         "some double"//char(0))

    ier = mnt_cmdlineargparser_setstring(prsr, "-s"//char(0), &
                                         "Hello there!"//char(0), &
                                         "some string"//char(0))

    ier = mnt_cmdlineargparser_setbool(prsr, "-v"//char(0), 0, &
                                         "verbose"//char(0))

    ier = mnt_cmdlineargparser_help(prsr)

    ! parse the command line arguments
    ier = mnt_cmdlineargparser_parse(prsr, nargs1, mnt_string_size, args(1))

    ! print help
    ier = mnt_cmdlineargparser_help(prsr)

    ! extract
    ier = mnt_cmdlineargparser_getint(prsr, "-i"//char(0), ival)
    ier = mnt_cmdlineargparser_getdouble(prsr, "-d"//char(0), dval)
    sval(:) = char(0)
    ier = mnt_cmdlineargparser_getstring(prsr, "-s"//char(0), sval, n)
    call mnt_c2f_string(sval, sval_f)
    ier = mnt_cmdlineargparser_getbool(prsr, "-v"//char(0), verbosity)

    print *, "-i arg is ", ival
    print *, "-d arg is ", dval
    print *, '-s arg is "'//trim(sval_f)//'"'
    print *, "-v arg is ", verbosity

    ! clean up
    deallocate(args)
    ier = mnt_cmdlineargparser_del(prsr)

end program 
