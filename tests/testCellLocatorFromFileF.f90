program test
    ! defines the C API
    use mnt_celllocator_capi_mod
    use mnt_cmdlineargparser_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double

    implicit none
    integer(c_size_t) :: cloc, prsr
    integer :: nargs, n, i
    character(len=mnt_string_size)  :: argv_full, filename, num_cells_per_bucket_str, out_filename
    character(len=mnt_string_size), allocatable :: args(:)
    integer(c_size_t), parameter :: num_cells = 1
    integer(c_int), parameter    :: num_verts_per_cell = 8
    integer(c_size_t), parameter :: num_points = num_cells*num_verts_per_cell
    real(c_double)       :: verts(num_points * 3)
    real(c_double)       :: target_point(3) = [ 4.3760449538518165_c_double, &
                                            &  -1.3135173587022644_c_double, &
                                            &   51392.815974060693_c_double]
    real(c_double)       :: interp_point(3)
    real(c_double)       :: pcoords(3)
    real(c_double)       :: diff2
    integer(c_int)    :: ier, num_cells_per_bucket
    integer(c_size_t) :: cell_id


    nargs = command_argument_count() + 1 ! must include the executable itself

    ier = mnt_cmdlineargparser_new(prsr)
    ier = mnt_cmdlineargparser_setstring(prsr, "-i"//char(0), "not-valid.vtk"//char(0), &
                                         "Input VTK file"//char(0))
    ier = mnt_cmdlineargparser_setint(prsr, "-n"//char(0), 100, &
                                         "Average number of cells per bucket"//char(0))
    ier = mnt_cmdlineargparser_setstring(prsr, "-o"//char(0), "testCellLocatorFromFile_grid.vtk"//char(0), &
                                         "Output VTK file"//char(0))

    ! parse
    allocate(args(nargs))
    do i = 1, nargs
        call get_command_argument(i, argv_full)
        ! add termination character, trim...
        call mnt_make_c_string(argv_full, args(i))
    enddo
    ier = mnt_cmdlineargparser_parse(prsr, nargs, len(args(1)), args)

    ! extract the arguments
    ier = mnt_cmdlineargparser_getint(prsr, "-n"//char(0), num_cells_per_bucket)
    ier = mnt_cmdlineargparser_getstring(prsr, "-i"//char(0), filename, n)
    ier = mnt_cmdlineargparser_getstring(prsr, "-o"//char(0), out_filename, n)

    ! done
    ier = mnt_cmdlineargparser_del(prsr)

    if (filename == 'not-valid.vtk') then
        stop 'ERROR: must provide VTK file name (-i <filename>)'
    endif

    print *,'reading grid from file >'//trim(filename)//'<'
    print *,'Number of cells per bucket = ', num_cells_per_bucket
    print *,'Output file name = >'//trim(out_filename)//'<'

    ier = mnt_celllocator_new(cloc)
    if(ier /= 0) print *,'ERROR after new ier = ', ier

    ier = mnt_celllocator_load(cloc, filename, len(trim(filename), c_size_t))
    if(ier /= 0) then 
        stop 'ERROR after load'
    endif

    ier = mnt_celllocator_build(cloc, num_cells_per_bucket)
    if(ier /= 0) print *,'ERROR ier = after build', ier

    ier = mnt_celllocator_rungriddiagnostics(cloc)
    if(ier /= 0) print *,'ERROR ier = after rungriddiagnostics', ier

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
    ier = mnt_celllocator_del(cloc)
    if(ier /= 0) print *,'ERROR after del ier = ', ier

    print *,'cell_id: ', cell_id, ' pcoords = ', pcoords

end program 
