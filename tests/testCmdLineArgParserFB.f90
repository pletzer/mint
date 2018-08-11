program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod
    use, intrinsic :: iso_c_binding, only: c_loc

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double

    implicit none
    integer(c_size_t)                                :: prsr
    integer(c_int)                                   :: ier, nargs, i, n, verbosity
    character(len=mnt_string_size), allocatable      :: args(:)
    character(len=mnt_string_size)                   :: argv_full, inp_filename, out_filename, inp_filename_f
    integer(c_int)                                   :: num_cells_per_bucket

    nargs = command_argument_count() + 1 ! must include the executable itself

    ier = mnt_cmdlineargparser_new(prsr)
    !                                     "Average number of cells per bucket"//char(0))
    ier = mnt_cmdlineargparser_setstring(prsr, "-i"//char(0), &
                                         "notvalid.vtk"//char(0), &
                                         "Input VTK file"//char(0))
    ! parse
    allocate(args(nargs))
    do i = 1, nargs
        ! make sure to include the name of the executable
        call get_command_argument(i - 1, argv_full)
        ! add termination character, trim...
        call mnt_f2c_string(argv_full, args(i))
    enddo
    ier = mnt_cmdlineargparser_help(prsr)

    ier = mnt_cmdlineargparser_parse(prsr, nargs, mnt_string_size, args)
    if (ier /= 0) then
        print *,'ERROR while parsing command line arguments'
    endif

    ier = mnt_cmdlineargparser_help(prsr)

    ! extract the arguments
    ier = mnt_cmdlineargparser_getstring(prsr, "-i"//char(0), inp_filename, n)
    call mnt_c2f_string(inp_filename, inp_filename_f)

    ! done
    ier = mnt_cmdlineargparser_del(prsr)

    print *, '-i arg is "'//trim(inp_filename_f)//'"'

    ! clean up
    deallocate(args)

end program 
