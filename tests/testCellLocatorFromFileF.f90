program test
    ! defines the C API
    use mnt_celllocator_capi_mod
    use mnt_cmdlineargparser_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_long_long, c_int, c_double, c_ptr, c_loc

    implicit none
    type(c_ptr)                                 :: prsr ! void*
    type(c_ptr)                                 :: cloc ! void*
    integer                                     :: nargs, nargs1, n, i, verbosity
    character(len=mnt_string_size)              :: argv_full,  &
                                                   num_cells_per_bucket_str, &
                                                   out_filename_f, inp_filename_f
    character(len=1), dimension(mnt_string_size):: inp_filename, out_filename
    character(len=1), allocatable               :: args(:)
    real(c_double)                              :: interp_point(3)
    real(c_double)                              :: pcoords(3)
    real(c_double)                              :: diff2, tol
    integer(c_int)                              :: ier, num_cells_per_bucket, num_bad_cells
    integer(c_long_long)                        :: cell_id

    real(c_double)          :: target_point(3) = [ 4.3760449538518165_c_double, &
                                               &  -1.3135173587022644_c_double, &
                                               &   51392.815974060693_c_double]

    nargs = command_argument_count()
    nargs1 = nargs + 1

    ! args must be a contiguous string
    allocate(args(nargs1 * mnt_string_size))

    do i = 0, nargs 
        ! i = 0 is the executable
        call get_command_argument(i, argv_full)
        ! add termination character, trim...
        call mnt_f2c_string(argv_full, args(i*mnt_string_size + 1:(i+1)*mnt_string_size))
    enddo

    ier = mnt_cmdlineargparser_new(prsr)
    if (ier /= 0) write(0, *) 'ERROR after mnt_cmdlineargparser_new'

    ier = mnt_cmdlineargparser_setint(prsr, "-n"//char(0), 100, &
                                         "Average number of cells per bucket"//char(0))
    if (ier /= 0) write(0, *) 'ERROR after mnt_cmdlineargparser_new'

    ier = mnt_cmdlineargparser_setstring(prsr, "-i"//char(0), &
                                         "notvalid.vtk"//char(0), &
                                         "Input VTK file"//char(0))
    if (ier /= 0) write(0, *) 'ERROR after mnt_cmdlineargparser_setstring -i'

    ier = mnt_cmdlineargparser_setstring(prsr, "-o"//char(0), &
                                         "testCellLocatorFromFile_grid.vtk"//char(0), &
                                         "Output VTK file"//char(0))
    if (ier /= 0) write(0, *) 'ERROR after mnt_cmdlineargparser_setstring -o'

    ier = mnt_cmdlineargparser_setbool(prsr, "-v"//char(0), 0, &
                                         "Turn verbosity on"//char(0))
    if (ier /= 0) write(0, *) 'ERROR after mnt_cmdlineargparser_setbool'

    ! parse the command line arguments
    ier = mnt_cmdlineargparser_parse(prsr, nargs1, mnt_string_size, args(1))
    if (ier /= 0) stop 'ERROR after mnt_cmdlineargparser_parse'

    ier = mnt_cmdlineargparser_help(prsr)
    if (ier /= 0) stop 'ERROR after mnt_cmdlineargparser_help'

    ! extract the arguments
    ier = mnt_cmdlineargparser_getint(prsr, "-n"//char(0), num_cells_per_bucket)
    if (ier /= 0) stop 'ERROR after mnt_cmdlineargparser_getint'

    inp_filename(:) = char(0) ! must be initialized to \n
    ier = mnt_cmdlineargparser_getstring(prsr, "-i"//char(0), inp_filename, n)
    if (ier /= 0) stop 'ERROR after mnt_cmdlineargparser_getstring -i'
    call mnt_c2f_string(inp_filename, inp_filename_f)

    out_filename(:) = char(0) ! must be initialized to \n
    ier = mnt_cmdlineargparser_getstring(prsr, "-o"//char(0), out_filename, n)
    if (ier /= 0) stop 'ERROR after mnt_cmdlineargparser_getstring -o'
    call mnt_c2f_string(out_filename, out_filename_f)

    ier = mnt_cmdlineargparser_getbool(prsr, "-v"//char(0), verbosity)
    if (ier /= 0) stop 'ERROR after mnt_cmdlineargparser_getbool'

    ! done
    ier = mnt_cmdlineargparser_del(prsr)
    if (ier /= 0) stop 'ERROR after mnt_cmdlineargparser_del'

    if (inp_filename_f == 'not-valid.vtk') then
        print *, 'ERROR: must provide VTK file name (-i <filename>)'
        stop 1
    endif

    write(0, *) 'Reading grid from file "'//trim(inp_filename_f)//'"'
    write(0, *) 'Number of cells per bucket = ', num_cells_per_bucket
    write(0, *) 'Output file name = "'//trim(out_filename_f)//'"'

    ier = mnt_celllocator_new(cloc)
    if(ier /= 0) then 
        write(0, *) 'ERROR after mnt_celllocator_new ier = ', ier
        stop
    endif

    print *, '*** inp_filename = "', inp_filename, '"'
    print *, '*** out_filename = "', out_filename, '"'
    ier = mnt_celllocator_load(cloc, inp_filename, len(trim(inp_filename_f), c_size_t))
    if(ier /= 0) then 
        write(0, *) 'ERROR after mnt_celllocator_load ier = ', ier
        stop 1
    endif

    ier = mnt_celllocator_build(cloc, num_cells_per_bucket)
    if(ier /= 0) then
        write(0, *) 'ERROR after mnt_celllocator_build ier = ', ier
    endif

    if (verbosity /= 0) then
        ier = mnt_celllocator_rungriddiagnostics(cloc)
        if(ier /= 0) write(0, *) 'ERROR ier = after mnt_celllocator_rungriddiagnostics ier = ', ier

        tol = -1.e-12_c_double
        ier = mnt_celllocator_checkgrid(cloc, tol, num_bad_cells)
        if(ier /= 0) write(0, *) 'ERROR ier = after mnt_celllocator_checkgrid ier = ', ier
        if (num_bad_cells > 0) then
          write(0, *) 'Number of bad cells: ', num_bad_cells
        endif
    endif

    ! initialiaze
    pcoords = -1._c_double
    ier = mnt_celllocator_find(cloc, target_point(1), cell_id, pcoords(1))
    if(ier /= 0) write(0, *) 'ERROR ier = after mnt_celllocator_find ier = ', ier

    if (cell_id >= 0) then
        ! the target point was found

        ier = mnt_celllocator_interppoint(cloc, cell_id, pcoords(1), interp_point(1))
        if(ier /= 0) write(0, *) 'ERROR ier = after mnt_celllocator_interppoint ier = ', ier

        ! check interpolation
        diff2 = dot_product(interp_point - target_point, interp_point - target_point)
        write(0, *) 'distance square error = ', diff2
    endif

    ier = mnt_celllocator_dumpgrid(cloc, out_filename, len(trim(out_filename_f), c_size_t))
    if(ier /= 0) write(0, *) 'ERROR after mnt_celllocator_dumpgrid ier = ', ier

    ! clean up
    deallocate(args)
    ier = mnt_celllocator_del(cloc)
    if(ier /= 0) write(0, *) 'ERROR after del ier = ', ier

    write(0, *) 'cell_id: ', cell_id, ' pcoords = ', pcoords

end program 
