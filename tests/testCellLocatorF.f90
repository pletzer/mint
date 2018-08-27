program test
    ! defines the C API
    use mnt_celllocator_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr

    implicit none
    type(c_ptr)                  :: cloc 
    integer(c_size_t), parameter :: num_cells = 1
    integer(c_int), parameter    :: num_verts_per_cell = 8
    integer(c_size_t), parameter :: num_points = num_cells*num_verts_per_cell
    real(c_double)               :: verts(num_points * 3)
    real(c_double), allocatable  :: target_point(:) 
    real(c_double)               :: interp_point(3)
    real(c_double)               :: pcoords(3)
    real(c_double)       :: diff2
    integer(c_int)       :: ier, num_cells_per_bucket
    integer(c_size_t)    :: cell_id

    verts = [0., 0., 0., &
             1., 0., 0., &
             1., 1., 0., &
             0., 1., 0., &
             0., 0., 1., &
             1., 0., 1., &
             1., 1., 1., &
             0., 1., 1.]

    ier = mnt_celllocator_new(cloc)
    if(ier /= 0) print*,'ERROR after new ier = ', ier

    ier = mnt_celllocator_setpoints(cloc, num_verts_per_cell, num_cells, verts(1))
    if(ier /= 0) print*,'ERROR after setpoints ier = ', ier

    num_cells_per_bucket = 512
    ier = mnt_celllocator_build(cloc, num_cells_per_bucket)
    if(ier /= 0) print*,'ERROR ier = after build', ier

    allocate(target_point(3))
    target_point = [0.2_c_double, 0.3_c_double, 0.4_c_double]
    ier = mnt_celllocator_find(cloc, target_point(1), cell_id, pcoords(1))
    if(ier /= 0) print*,'ERROR ier = after find', ier

    ier = mnt_celllocator_interp_point(cloc, cell_id, pcoords(1), interp_point(1))
    if(ier /= 0) print*,'ERROR ier = after find', ier
    ! check interpolation
    diff2 = dot_product(interp_point - target_point, interp_point - target_point)
    print *,'distance square error = ', diff2
    deallocate(target_point)

    ier = mnt_celllocator_del(cloc)
    if(ier /= 0) print*,'ERROR after del ier = ', ier

    print *,'cell_id: ', cell_id, ' pcoords = ', pcoords

end program 
