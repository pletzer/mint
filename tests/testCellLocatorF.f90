program test
    ! defines the C API
    use mnt_celllocator_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double

    implicit none
    integer(c_size_t) :: cloc
    integer(c_size_t), parameter :: num_cells = 1
    integer(c_int), parameter    :: num_verts_per_cell = 8
    integer(c_size_t), parameter :: num_points = num_cells*num_verts_per_cell
    real(c_double)       :: verts(num_points * 3)
    real(c_double)       :: target_point(3) = [0.2_c_double, 0.3_c_double, 0.4_c_double]
    real(c_double)       :: pcoords(3)
    integer(c_int)    :: ier
    integer(c_size_t) :: cell_id

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
    ier = mnt_celllocator_build(cloc)
    if(ier /= 0) print*,'ERROR ier = after build', ier
    ier = mnt_celllocator_find(cloc, target_point(1), cell_id, pcoords(1))
    if(ier /= 0) print*,'ERROR ier = after find', ier
    ier = mnt_celllocator_del(cloc)
    if(ier /= 0) print*,'ERROR after del ier = ', ier

    print *,'cell_id: ', cell_id, ' pcoords = ', pcoords

end program 
