program test_cell_locator_capi
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_size_t
    use mnt_celllocator_capi_mod,       only: mnt_celllocator_new, &
                                              mnt_celllocator_del, &
                                              mnt_celllocator_setpointsptr, &
                                              mnt_celllocator_build, &
                                              mnt_celllocator_find, &
                                              mnt_celllocator_dumpgrid, &
                                              mnt_celllocator_interppoint
    implicit none
    integer              :: ier, nverts_per_cell, ndims, num_cells_per_bucket
    integer(8)           :: ncells, cell_id, nx, ny, i, j, icell, index0, index1, index2, index3
    real(8), allocatable :: verts(:)
    real(8)              :: point(3), pcoords(3), dx, dy, x0, y0, x1, y1
    type(c_ptr)          :: handle

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
    do j = 1, ny
        do i = 1, nx
            x0 = 0._8 + dx*i
            y0 = -90._8 + dy*j
            x1 = x0 + dx
            y1 = y0 + dy
            index0 = ndims * nverts_per_cell * icell
            index1 = index0 + ndims
            index2 = index1 + ndims
            index3 = index2 + ndims
            verts(index0:index0 + ndims - 1) = (/x0, y0, 0._8/)
            verts(index1:index1 + ndims - 1) = (/x1, y0, 0._8/)
            verts(index2:index2 + ndims - 1) = (/x1, y1, 0._8/)
            verts(index3:index3 + ndims - 1) = (/x0, y1, 0._8/)
            icell = icell + 1
        enddo
    enddo

    ier = mnt_celllocator_new(handle)
    if (ier /= 0) stop 1

    ier = mnt_celllocator_setpointsptr(handle, nverts_per_cell, ncells, verts(1))
    if (ier /= 0) stop 3

    num_cells_per_bucket = 2
    ier = mnt_celllocator_build(handle, num_cells_per_bucket)
    if (ier /= 0) stop 4

    point = (/0.9_8, 0.4_8, 0.0_8/)
    ier = mnt_celllocator_find(handle, point, cell_id, pcoords)
    if (ier /= 0) stop 5

    print*,'cell_id = ', cell_id, ' for point ', point, ' (pcoords = ', pcoords, ')'

    ier = mnt_celllocator_del(handle)
    if (ier /= 0) stop 2

    deallocate(verts)


end program test_cell_locator_capi