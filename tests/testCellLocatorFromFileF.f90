program test
    ! defines the C API
    use mnt_celllocator_capi_mod
    use mnt_cmdlineargparser_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double

    implicit none
    integer(c_size_t)                           :: cloc, prsr
    integer                                     :: nargs, n, i, verbosity
    character(len=mnt_string_size)              :: argv_full, inp_filename, &
                                                   num_cells_per_bucket_str, out_filename, &
                                                   out_filename_f, inp_filename_f
    character(len=mnt_string_size), allocatable :: args(:)
    real(c_double)                              :: interp_point(3)
    real(c_double)                              :: pcoords(3)
    real(c_double)                              :: diff2
    integer(c_int)                              :: ier, num_cells_per_bucket
    integer(c_size_t)                           :: cell_id

    real(c_double)       :: target_point(3) = [ 4.3760449538518165_c_double, &
                                            &  -1.3135173587022644_c_double, &
                                            &   51392.815974060693_c_double]


    ier = mnt_cmdlineargparser_new(prsr)
    ier = mnt_cmdlineargparser_setint(prsr, "-n"//char(0), 100, &
                                         "Average number of cells per bucket"//char(0))
    ier = mnt_cmdlineargparser_setstring(prsr, "-i"//char(0), &
                                         "notvalid.vtk"//char(0), &
                                         "Input VTK file"//char(0))
    ier = mnt_cmdlineargparser_setstring(prsr, "-o"//char(0), &
                                         "testCellLocatorFromFile_grid.vtk"//char(0), &
                                         "Output VTK file"//char(0))
    ier = mnt_cmdlineargparser_setbool(prsr, "-v"//char(0), 0, &
                                         "Turn verbosity on"//char(0))
    ! parse
    nargs = command_argument_count() + 1 ! must include the executable's name
    allocate(args(nargs))
    do i = 1, nargs
        ! include the executable
        call get_command_argument(i - 1, argv_full)
        ! add termination character, trim...
        call mnt_f2c_string(argv_full, args(i))
    enddo
    ier = mnt_cmdlineargparser_help(prsr)

    ier = mnt_cmdlineargparser_parse(prsr, nargs, mnt_string_size, args)

    ier = mnt_cmdlineargparser_help(prsr)

    ! extract the arguments
    ier = mnt_cmdlineargparser_getint(prsr, "-n"//char(0), num_cells_per_bucket)
    ier = mnt_cmdlineargparser_getstring(prsr, "-i"//char(0), inp_filename, n)
    call mnt_c2f_string(inp_filename, inp_filename_f)
    ier = mnt_cmdlineargparser_getstring(prsr, "-o"//char(0), out_filename, n)
    call mnt_c2f_string(out_filename, out_filename_f)
    ier = mnt_cmdlineargparser_getbool(prsr, "-v"//char(0), verbosity)

    ! done
    ier = mnt_cmdlineargparser_del(prsr)

    if (inp_filename == 'not-valid.vtk') then
        stop 'ERROR: must provide VTK file name (-i <filename>)'
    endif

    print *,'Reading grid from file "'//trim(inp_filename_f)//'"'
    print *,'Number of cells per bucket = ', num_cells_per_bucket
    print *,'Output file name = "'//trim(out_filename_f)//'"'

    ier = mnt_celllocator_new(cloc)
    if(ier /= 0) print *,'ERROR after new ier = ', ier

    ier = mnt_celllocator_load(cloc, inp_filename, len(trim(inp_filename), c_size_t))
    if(ier /= 0) then 
        stop 'ERROR after load'
    endif

    ier = mnt_celllocator_build(cloc, num_cells_per_bucket)
    if(ier /= 0) print *,'ERROR ier = after build', ier

    if (verbosity /= 0) then
        ier = mnt_celllocator_rungriddiagnostics(cloc)
        if(ier /= 0) print *,'ERROR ier = after rungriddiagnostics', ier
    endif

    ! initialiaze
    pcoords = 0._c_double
    ier = mnt_celllocator_find(cloc, target_point(1), cell_id, pcoords(1))
    if(ier /= 0) print *,'ERROR ier = after find', ier

    if (cell_id >= 0) then
        ! the target point was found

        ier = mnt_celllocator_interp_point(cloc, cell_id, pcoords(1), interp_point(1))
        if(ier /= 0) print *,'ERROR ier = after find', ier

        ! check interpolation
        diff2 = dot_product(interp_point - target_point, interp_point - target_point)
        print *,'distance square error = ', diff2
    endif

    ier = mnt_celllocator_dumpgrid(cloc, trim(out_filename), len(trim(out_filename), c_size_t))
    if(ier /= 0) print *,'ERROR after dumpgrid ier = ', ier

    ! clean up
    deallocate(args)
    ier = mnt_celllocator_del(cloc)
    if(ier /= 0) print *,'ERROR after del ier = ', ier

    print *,'cell_id: ', cell_id, ' pcoords = ', pcoords

end program 
