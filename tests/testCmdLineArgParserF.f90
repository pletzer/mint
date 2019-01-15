program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_ptr

    implicit none
    type(c_ptr)                                      :: prsr
    integer(c_int)                                   :: ier, nargs, i, n, verbosity, nargs1
    character(len=1), allocatable                    :: args(:)
    character(len=mnt_string_size)                   :: argv_full
    integer(c_int)                                   :: nval, vval
    real(c_double)                                   :: dval
    character(len=1), dimension(mnt_string_size)     :: ival, oval
    character(len=mnt_string_size)                   :: ival_f, oval_f


    nargs = command_argument_count()
    nargs1 = nargs + 1

    ! args must be a contiguous string
    allocate(args(nargs1 * mnt_string_size))
    args(:) = char(0)

    do i = 0, nargs 
        ! i = 0 is the executable
        call get_command_argument(i, argv_full)
        ! add termination character, trim...
        call mnt_f2c_string(argv_full, args(i*mnt_string_size + 1:(i+1)*mnt_string_size))
    enddo

    ! define the command line arguments
    ier = mnt_cmdlineargparser_new(prsr)

    ! set
    ier = mnt_cmdlineargparser_setint(prsr, "-n"//char(0), -123, &
                                      "some integer"//char(0))

    ier = mnt_cmdlineargparser_setdouble(prsr, "-d"//char(0), -1.23_c_double, &
                                         "some double"//char(0))

    ier = mnt_cmdlineargparser_setstring(prsr, "-i"//char(0), &
                                         "Hello input!"//char(0), &
                                         "some input string"//char(0))

    ier = mnt_cmdlineargparser_setstring(prsr, "-o"//char(0), &
                                         "Hello output!"//char(0), &
                                         "some output string"//char(0))

    ier = mnt_cmdlineargparser_setbool(prsr, "-v"//char(0), 0, &
                                         "verbose"//char(0))

    ! print help
    ier = mnt_cmdlineargparser_help(prsr)

    ! parse the command line arguments
    ier = mnt_cmdlineargparser_parse(prsr, nargs1, mnt_string_size, args(1))

    ! extract
    ier = mnt_cmdlineargparser_getint(prsr, "-n"//char(0), nval)
    ier = mnt_cmdlineargparser_getstring(prsr, "-i"//char(0), size(ival), ival)
    ier = mnt_cmdlineargparser_getstring(prsr, "-o"//char(0), size(oval), oval)
    ier = mnt_cmdlineargparser_getdouble(prsr, "-d"//char(0), dval)
    ier = mnt_cmdlineargparser_getbool(prsr, "-v"//char(0), verbosity)

    call mnt_c2f_string(ival, ival_f)
    call mnt_c2f_string(oval, oval_f)

    print *, "-n arg is ", nval
    print *, "-i arg is '", trim(ival_f), "'"
    print *, "-o arg is '", trim(oval_f), "'"
    print *, "-d arg is ", dval
    print *, "-v arg is ", verbosity

    ! clean up
    deallocate(args)
    ier = mnt_cmdlineargparser_del(prsr)

end program 
