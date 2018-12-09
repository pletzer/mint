module mnt_regridedges_capi_mod
  ! C function prototypes
  interface

    function mnt_regridedges_new(obj) &
                              &  bind(C, name='mnt_regridedges_new')
      ! Constructor
      ! @param obj opaque handle
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_regridedges_new
    end function mnt_regridedges_new

    function mnt_regridedges_del(obj) &
                               & bind(C, name='mnt_regridedges_del')
      ! Destructor
      ! @param obj opaque handle
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_regridedges_del
    end function mnt_regridedges_del

    function mnt_regridedges_load(obj, filename, n) &
                                & bind(C, name='mnt_regridedges_load')
      ! Load weights from file 
      ! @param obj opaque handle
      ! @param filename file name
      ! @return 0 if sucecssful
      ! @note must invoke constructor prior to this call
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_size_t), value                 :: n
      integer(c_int)                           :: mnt_regridedges_load
    end function mnt_regridedges_load

    function mnt_regridedges_setSrcPointsPtr(obj, nverts_per_cell, ncells, verts) &
                                          &  bind(C, name='mnt_regridedges_setSrcPointsPtr')
      ! Set the source grid points
      ! @param obj opaque handle
      ! @param nverts_per_cell number of vertices per cell
      ! @param ncells number of cells
      ! @param verts flat array of vertices [x0, y0, z0, x1, y1, z1, ...]
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_int), value                    :: nverts_per_cell
      integer(c_size_t), value                 :: ncells
      real(c_double), intent(in)               :: verts(*) ! const double*
      integer(c_int)                           :: mnt_regridedges_setsrcpointsptr
    end function mnt_regridedges_setSrcPointsPtr

    function mnt_regridedges_setDstPointsPtr(obj, nverts_per_cell, ncells, verts) &
                                          &  bind(C, name='mnt_regridedges_setDstPointsPtr')
      ! Set the destination grid points
      ! @param obj opaque handle
      ! @param nverts_per_cell number of vertices per cell
      ! @param ncells number of cells
      ! @param verts flat array of vertices [x0, y0, z0, x1, y1, z1, ...]
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_int), value                    :: nverts_per_cell
      integer(c_size_t), value                 :: ncells
      real(c_double), intent(in)               :: verts(*) ! const double*
      integer(c_int)                           :: mnt_regridedges_setdstpointsptr
    end function mnt_regridedges_setDstPointsPtr

    function mnt_regridedges_build(obj, num_cells_per_bucket) &
                                 & bind(C, name='mnt_regridedges_build')
      ! Build the regridder
      ! @param obj opaque handle
      ! @param num_cells_per_bucket number of cells per bucket
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: num_cells_per_bucket
      integer(c_int)                   :: mnt_regridedges_build
    end function mnt_regridedges_build

    function mnt_regridedges_applyWeights(obj, src_data, dst_data) &
                                & bind(C, name='mnt_regridedges_applyWeights')
      ! Apply the interpolation weights
      ! @param obj opaque handle
      ! @param src_data edge centred data on the source grid 
      ! @param dst_data edge centred data on the destination grid 
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_long_long, c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      real(c_double), intent(in)               :: src_data(*) ! const double*
      real(c_double), intent(out)              :: dst_data(*) ! double*
      integer(c_int)                           :: mnt_regridedges_applyWeights
    end function mnt_regridedges_applyWeights

    function mnt_regridedges_dump(obj, filename, n) & 
                                & bind(C, name='mnt_regridedges_dump')
      ! Dump the weights to a file
      ! @param obj instance of mntregridedges_t (opaque handle)
      ! @param filename file name, \0 (char(0)) terminated
      ! @param n length of filename 
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_size_t), value                 :: n
      integer(c_int)                           :: mnt_regridedges_dump  
    end function mnt_regridedges_dump

  end interface

end module mnt_regridedges_capi_mod

