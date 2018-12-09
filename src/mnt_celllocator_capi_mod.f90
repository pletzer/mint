module mnt_celllocator_capi_mod
  ! C function prototypes
  interface

    function mnt_celllocator_new(obj) &
                              &  bind(C, name='mnt_celllocator_new')
      ! Constructor
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_celllocator_new
    end function mnt_celllocator_new

    function mnt_celllocator_del(obj) &
                               & bind(C, name='mnt_celllocator_del')
      ! Destructor
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_celllocator_del
    end function mnt_celllocator_del

    function mnt_celllocator_load(obj, filename, n) &
                                & bind(C, name='mnt_celllocator_load')
      ! Load cell locator object from file 
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @param filename file name
      ! @return 0 if sucecssful
      ! @note must invoke constructor prior to this call
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_size_t), value                 :: n
      integer(c_int)                           :: mnt_celllocator_load
    end function mnt_celllocator_load

    function mnt_celllocator_setPointsPtr(obj, nverts_per_cell, ncells, verts) &
                                       &  bind(C, name='mnt_celllocator_setPointsPtr')
      ! Set the grid points
      ! @param obj instance of mntcellLocator_t (opaque handle)
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
      integer(c_int)                           :: mnt_celllocator_setpointsptr
    end function mnt_celllocator_setPointsPtr

    function mnt_celllocator_build(obj, num_cells_per_bucket) &
                                 & bind(C, name='mnt_celllocator_build')
      ! Build locator object
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @param num_cells_per_bucket number of cells per bucket
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: num_cells_per_bucket
      integer(c_int)                   :: mnt_celllocator_build
    end function mnt_celllocator_build

    function mnt_celllocator_find(obj, point, cell_id, pcoords) &
                                & bind(C, name='mnt_celllocator_find')
      ! Find point
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @param point 3d point
      ! @param cell_id output cell Id (zero-based)
      ! @param output parametric coordinates
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_long_long, c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      real(c_double), intent(in)               :: point(3) ! const double*
      integer(c_long_long), intent(out)        :: cell_id
      real(c_double), intent(out)              :: pcoords(3) ! double*
      integer(c_int)                           :: mnt_celllocator_find
    end function mnt_celllocator_find

    function mnt_celllocator_interpPoint(obj, cell_id, pcoords, point) &
                                      &  bind(C, name='mnt_celllocator_interpPoint')
      ! Interpolate point
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @param cell_id 0-based cell Id
      ! @param pcoords array of parametric coordinates
      ! @param point output point
      ! @return 0 if successful
      ! @note call mnt_celllocator_find to compute pcoords and get cell_id
      use, intrinsic :: iso_c_binding, only: c_long_long, c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_long_long), value              :: cell_id
      real(c_double), intent(in)               :: pcoords(3) ! const double*
      real(c_double), intent(out)              :: point(3)   ! double*
      integer(c_int)                           :: mnt_celllocator_interppoint
    end function mnt_celllocator_interpPoint

    function mnt_celllocator_dumpGrid(obj, filename, n) & 
                                    & bind(C, name='mnt_celllocator_dumpGrid')
      ! Dump grid to VTK file
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @param filename file name, \0 (char(0)) terminated
      ! @param n length of filename 
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_size_t), value                 :: n
      integer(c_int)                           :: mnt_celllocator_dumpgrid      
    end function mnt_celllocator_dumpGrid

    function mnt_celllocator_checkGrid(obj, tol, num_bad_cells) &
                                     & bind(C, name='mnt_celllocator_checkGrid')
      ! Check/detect grid cells with bad topology
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @param tol min area/volume tolerance (typically small and negative)
      ! @param num_bad_cells number of bad cells, ie cells with area/volume < tol
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      real(c_double), value                    :: tol
      integer(c_int), intent(out)              :: num_bad_cells
      integer(c_int)                           :: mnt_celllocator_checkgrid   
    end function mnt_celllocator_checkGrid

    function mnt_celllocator_runGridDiagnostics(obj) &
                                              & bind(C, name='mnt_celllocator_runGridDiagnostics')
      ! Run grid diagnostics
      ! @param obj instance of mntcellLocator_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_int)                           :: mnt_celllocator_rungriddiagnostics    
    end function mnt_celllocator_runGridDiagnostics

    subroutine mnt_celllocator_printAddress(something) &
                                          & bind(C, name='mnt_celllocator_printAddress')
      ! Print the address of something
      use, intrinsic :: iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: something
    end subroutine mnt_celllocator_printAddress

  end interface

end module mnt_celllocator_capi_mod

