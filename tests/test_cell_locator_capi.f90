program test_cell_locator_capi
    use, intrinsic :: iso_c_binding, only: c_long_long, c_int, c_double, c_ptr, c_size_t
    use mnt_celllocator_capi_mod,       only: mnt_celllocator_new, &
                                              mnt_celllocator_del, &
                                              mnt_celllocator_setpointsptr, &
                                              mnt_celllocator_build, &
                                              mnt_celllocator_find, &
                                              mnt_celllocator_dumpgrid, &
                                              mnt_celllocator_interppoint
    implicit none
    integer              :: ier, nverts_per_cell, ndims, num_cells_per_bucket
    integer              :: ncells, nx, ny, i, j, icell, index0, index1, index2, index3, index4
    integer(c_long_long) :: cell_id
    real(8), allocatable :: verts(:)
    real(8)              :: point(3), pcoords(3), dx, dy, x0, y0, x1, y1
    type(c_ptr)          :: handle
    character(len=521)   :: output_file

    ! create a uniform grid
    nx = 20
    ny = 10
    ncells = nx * ny
    nverts_per_cell = 4
    ndims = 3
    allocate(verts(ndims * nverts_per_cell * ncells))
    dx = 360._8 / real(nx, 8)
    dy = 180._8 / real(ny, 8)
    icell = 0
    do j = 0, ny - 1
        do i = 0, nx - 1
            x0 = 0._8 + dx*i
            y0 = -90._8 + dy*j
            x1 = x0 + dx
            y1 = y0 + dy
            print*,'icell = ', icell, 'x0, y0 = ', x0, y0, ' x1, y1 = ', x1, y1
            index0 = ndims*nverts_per_cell*icell + 1
            index1 = index0 + ndims
            index2 = index1 + ndims
            index3 = index2 + ndims
            index4 = index3 + ndims
            verts(index0:index1 - 1) = (/x0, y0, 0._8/)
            verts(index1:index2 - 1) = (/x1, y0, 0._8/)
            verts(index2:index3 - 1) = (/x1, y1, 0._8/)
            verts(index3:index4 - 1) = (/x0, y1, 0._8/)
            icell = icell + 1
        enddo
    enddo

    ier = mnt_celllocator_new(handle)
    if (ier /= 0) print *,'ERROR: mnt_celllocator_new'

    ier = mnt_celllocator_setpointsptr(handle, nverts_per_cell, &
                                       int(ncells, kind=c_size_t), verts(1))
    if (ier /= 0) print *, 'ERROR: mnt_celllocator_setpointsptr'

    output_file = 'test_cell_locator_capi_grid.vtk'
    ier = mnt_celllocator_dumpgrid(handle, trim(output_file), int(len_trim(output_file), kind=8))
    if (ier /= 0) print *, 'ERROR: mnt_celllocator_dumpgrid'

    num_cells_per_bucket = 1
    ier = mnt_celllocator_build(handle, num_cells_per_bucket)
    if (ier /= 0) print *, 'ERROR: mnt_celllocator_build'

    point = (/10.0_8, 20.0_8, 0.0_8/)
    ier = mnt_celllocator_find(handle, point, cell_id, pcoords)
    if (ier /= 0) print *, 'ERROR: mnt_celllocator_find'

    print*,'cell_id = ', cell_id, ' for point ', point, ' (pcoords = ', pcoords, ')'
    if (cell_id < 0) then
        print*,'ERROR: could not find point'
    endif

    ier = mnt_celllocator_del(handle)
    if (ier /= 0) print *, 'ERROR: mnt_celllocator_del'

    deallocate(verts)

end program test_cell_locator_capi
