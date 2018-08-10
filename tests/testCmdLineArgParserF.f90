program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod
    use, intrinsic :: iso_c_binding, only: c_loc

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double

    implicit none
    integer(c_size_t)               :: prsr
    integer(c_int)                  :: ier, nargs, i, n
    character(len=mnt_string_size), allocatable     :: args(:)
    character(len=mnt_string_size)                  :: argv_full
    integer(c_int)                  :: ival
    real(c_double)                  :: dval
    character(len=mnt_string_size)                    :: sval
    character(len=mnt_string_size)                    :: s1, s2, s3

    ier = mnt_cmdlineargparser_new(prsr)

    ! set
    call mnt_make_c_string("-i", s1)
    call mnt_make_c_string("some integer", s3)
    ier = mnt_cmdlineargparser_setint(prsr, s1, -123, s3)

    call mnt_make_c_string("-d", s1)
    call mnt_make_c_string("some double", s3)
    ier = mnt_cmdlineargparser_setdouble(prsr, s1, -1.23_c_double, s3)

    call mnt_make_c_string("-i", s1)
    call mnt_make_c_string("Hello there!", s2)
    call mnt_make_c_string("some integer", s3)
    ier = mnt_cmdlineargparser_setstring(prsr, s1, s2, s3)

    ! parse
    nargs = command_argument_count() + 1
    allocate(args(nargs))
    do i = 1, nargs
        call get_command_argument(i, argv_full)
        ! add termination character, trim...
        call mnt_make_c_string(argv_full, args(i))
    enddo
    n = len(args(1))
    write(0, *) '*** 4 nargs = ', nargs, ' n = ', n, ' args = ', args
    ier = mnt_cmdlineargparser_parse(prsr, nargs, n, args)

    ! extract
    write(0, *) '*** 5'
    call mnt_make_c_string("-i", s1)
    ier = mnt_cmdlineargparser_getint(prsr, s1, ival)
    call mnt_make_c_string("-d", s1)
    ier = mnt_cmdlineargparser_getdouble(prsr, s1, dval)
    call mnt_make_c_string("-s", s1)
    ier = mnt_cmdlineargparser_getstring(prsr, s1, sval, n)

    print *, "-i arg is ", ival
    print *, "-d arg is ", dval
    print *, "-s arg is ", sval

    ! clean up
    deallocate(args)
    ier = mnt_cmdlineargparser_del(prsr)

end program 
