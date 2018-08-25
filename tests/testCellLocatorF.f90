program test
    ! defines the C API
    use mnt_celllocator_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr, c_loc

    implicit none
    integer(c_size_t) :: cloc
    integer(c_size_t), parameter :: num_cells = 1
    integer(c_int), parameter    :: num_verts_per_cell = 8
    integer(c_size_t), parameter :: num_points = num_cells*num_verts_per_cell
    real(c_double), target       :: verts(num_points * 3)
    real(c_double), target       :: target_point(3) = [0.2_c_double, 0.3_c_double, 0.4_c_double]
    real(c_double), target       :: interp_point(3)
    real(c_double), target       :: pcoords(3)
    real(c_double)       :: diff2
    integer(c_int)       :: ier, num_cells_per_bucket
    integer(c_size_t)    :: cell_id
    type(c_ptr)          :: pcoords_ptr, interp_point_ptr

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

    ier = mnt_celllocator_setpoints(cloc, num_verts_per_cell, num_cells, c_loc(verts(1)))
    if(ier /= 0) print*,'ERROR after setpoints ier = ', ier

    num_cells_per_bucket = 512
    ier = mnt_celllocator_build(cloc, num_cells_per_bucket)
    if(ier /= 0) print*,'ERROR ier = after build', ier

    pcoords_ptr = c_loc(pcoords(1))
    ier = mnt_celllocator_find(cloc, c_loc(target_point(1)), cell_id, pcoords_ptr)
    if(ier /= 0) print*,'ERROR ier = after find', ier

    interp_point_ptr = c_loc(interp_point(1))
    ier = mnt_celllocator_interp_point(cloc, cell_id, c_loc(pcoords(1)), interp_point_ptr)
    if(ier /= 0) print*,'ERROR ier = after find', ier
    ! check interpolation
    diff2 = dot_product(interp_point - target_point, interp_point - target_point)
    print *,'distance square error = ', diff2

    ier = mnt_celllocator_del(cloc)
    if(ier /= 0) print*,'ERROR after del ier = ', ier

    print *,'cell_id: ', cell_id, ' pcoords = ', pcoords

end program 
