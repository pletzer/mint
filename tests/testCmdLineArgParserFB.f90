program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod
    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr

    implicit none
    type(c_ptr)                                      :: prsr
    integer(c_int)                                   :: ier, nargs, i, n, verbosity, nargs1
    character(len=1), allocatable                    :: args(:)
    character(len=mnt_string_size)                   :: argv_full, inp_filename, out_filename, inp_filename_f
    integer(c_int)                                   :: num_cells_per_bucket

    ! parse
    nargs = command_argument_count() ! excludes the executable
    nargs1 = nargs + 1               ! includes the executable

    allocate(args(nargs1*mnt_string_size))

    do i = 0, nargs
        ! make sure to include the name of the executable
        call get_command_argument(i, argv_full)
        ! add termination character, trim...
        call mnt_f2c_string(argv_full, args(i*mnt_string_size + 1:(i+1)*mnt_string_size))
    enddo

    ier = mnt_cmdlineargparser_new(prsr)
    if (ier /= 0) then
        print *,'ERROR after calling mnt_cmdlineargparser_new'
    endif

    ier = mnt_cmdlineargparser_setstring(prsr, "-i"//char(0), &
                                               "notvalid.vtk"//char(0), &
                                               "Input VTK file"//char(0))
    if (ier /= 0) then
        print *,'ERROR after calling mnt_cmdlineargparser_setstring'
    endif

    ier = mnt_cmdlineargparser_parse(prsr, nargs1, mnt_string_size, args(1))
    if (ier /= 0) then
        print *,'ERROR after calling mnt_cmdlineargparser_parse'
    endif

    ier = mnt_cmdlineargparser_help(prsr)
    if (ier /= 0) then
        print *,'ERROR after calling mnt_cmdlineargparser_help'
    endif

    ! extract the arguments
    ier = mnt_cmdlineargparser_getstring(prsr, "-i"//char(0), inp_filename, n)
    if (ier /= 0) then
        print *,'ERROR after calling mnt_cmdlineargparser_getstring'
    endif
    call mnt_c2f_string(inp_filename, inp_filename_f)

    ! done
    ier = mnt_cmdlineargparser_del(prsr)
    if (ier /= 0) then
        print *,'ERROR after calling mnt_cmdlineargparser_del'
    endif

    print *, '-i arg is "'//trim(inp_filename_f)//'"'

    ! clean up
    deallocate(args)

end program 
