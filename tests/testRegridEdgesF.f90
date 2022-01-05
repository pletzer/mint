program test
    ! defines the C API
    use mnt_regridedges_capi_mod
    use mnt_cmdlineargparser_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr

    implicit none
    type(c_ptr)                  :: prsr
    type(c_ptr)                  :: crg
    integer                      :: num_cells_per_bucket, ier, nargs, nargs1, i, debug, enablefolding
    character(len=1), allocatable                    :: args(:)
    character(len=mnt_string_size)                   :: argv_full
    character(len=1), dimension(mnt_string_size)     :: src_filename, dst_filename, output_filename
    character(len=mnt_string_size)                   :: src_filename_f, dst_filename_f, output_filename_f
    real(8)                                          :: periodx

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

    ier = mnt_cmdlineargparser_setstring(prsr, "-s"//char(0), &
                                         "cs_16.nc"//char(0), &
                                         "Set the source ugrid file"//char(0))

    ier = mnt_cmdlineargparser_setstring(prsr, "-d"//char(0), &
                                         "cs_4.nc"//char(0), &
                                         "Set the destination ugrid file"//char(0))

    ier = mnt_cmdlineargparser_setstring(prsr, "-o"//char(0), &
                                         "weights.nc"//char(0), &
                                         "Set the output netCDF file holding the weights"//char(0))
    ! parse the command line arguments
    ier = mnt_cmdlineargparser_parse(prsr, nargs1, mnt_string_size, args(1))

    ! extract
    ier = mnt_cmdlineargparser_getstring(prsr, "-s"//char(0), size(src_filename), src_filename)
    ier = mnt_cmdlineargparser_getstring(prsr, "-d"//char(0), size(dst_filename), dst_filename)
    ier = mnt_cmdlineargparser_getstring(prsr, "-o"//char(0), size(output_filename), output_filename)

    ! convert to fortran strings
    call mnt_c2f_string(src_filename, src_filename_f)
    call mnt_c2f_string(dst_filename, dst_filename_f)
    call mnt_c2f_string(output_filename, output_filename_f)

    print *,' src  ugrid: ', trim(src_filename_f)
    print *,' dst  ugrid: ', trim(dst_filename_f)
    print *,' out netcdf: ', trim(output_filename_f)

    ! now do the regridding

    ier = mnt_regridedges_new(crg)
    if(ier /= 0) print*,'ERROR after new ier = ', ier

    ier = mnt_regridedges_loadSrcGrid(crg, src_filename_f, len_trim(src_filename_f))
    if(ier /= 0) print*,'ERROR after loadsrc ier = ', ier

    ier = mnt_regridedges_loadDstGrid(crg, dst_filename_f, len_trim(dst_filename_f))
    if(ier /= 0) print*,'ERROR after loaddst ier = ', ier

    num_cells_per_bucket = 8
    periodx = 0._8
    enablefolding = 0
    ier = mnt_regridedges_buildLocator(crg, num_cells_per_bucket, periodx, enablefolding)
    if(ier /= 0) print*,'ERROR ier = after build', ier
    
    debug = 1
    ier = mnt_regridedges_computeWeights(crg, debug)
    if(ier /= 0) print*,'ERROR ier = after computeWeights', ier

    ier = mnt_regridedges_dumpWeights(crg, output_filename_f, len_trim(output_filename_f))
    if(ier /= 0) print*,'ERROR ier = after dump', ier

    ier = mnt_regridedges_del(crg)
    if(ier /= 0) print*,'ERROR after del ier = ', ier

end program test
