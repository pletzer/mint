module mnt_celllocator_capi_mod
  ! C function prototypes
  interface

    function mnt_celllocator_new(obj) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout) :: obj
      integer(c_int)                   :: mnt_celllocator_new
    end function mnt_celllocator_new

    function mnt_celllocator_del(obj) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout) :: obj
      integer(c_int)                   :: mnt_celllocator_del
    end function mnt_celllocator_del

    function mnt_celllocator_load(obj, filename, n) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: filename
      integer(c_size_t), value                 :: n
      integer(c_int)                           :: mnt_celllocator_load
    end function mnt_celllocator_load

    function mnt_celllocator_setpoints(obj, nverts_per_cell, ncells, verts) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr
      integer(c_size_t), intent(inout)         :: obj
      integer(c_int), value                    :: nverts_per_cell
      integer(c_size_t), value                 :: ncells
      type(c_ptr), intent(in)                  :: verts ! const double*
      integer(c_int)                           :: mnt_celllocator_setpoints
    end function mnt_celllocator_setpoints

    function mnt_celllocator_build(obj, num_cells_per_bucket) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout) :: obj
      integer(c_int), intent(in)       :: num_cells_per_bucket
      integer(c_int)                   :: mnt_celllocator_build
    end function mnt_celllocator_build

    function mnt_celllocator_find(obj, point, cell_id, pcoords) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr
      integer(c_size_t), intent(inout)         :: obj
      type(c_ptr), intent(in)                  :: point   ! const double*
      integer(c_size_t), intent(out)           :: cell_id
      type(c_ptr), intent(out)                 :: pcoords ! double*
      integer(c_int)                           :: mnt_celllocator_find
    end function mnt_celllocator_find

    function mnt_celllocator_interp_point(obj, cell_id, pcoords, point) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr
      integer(c_size_t), intent(inout)         :: obj
      integer(c_size_t), value                 :: cell_id
      type(c_ptr), intent(in)                  :: pcoords ! const double*
      type(c_ptr), intent(out)                 :: point   ! double*
      integer(c_int)                           :: mnt_celllocator_find
    end function mnt_celllocator_interp_point

    function mnt_celllocator_dumpgrid(obj, filename, n) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      character(len=1), intent(in)             :: filename
      integer(c_size_t), value                 :: n
      integer(c_int)                           :: mnt_celllocator_dumpgrid      
    end function mnt_celllocator_dumpgrid

    function mnt_celllocator_rungriddiagnostics(obj) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      integer(c_int)                           :: mnt_celllocator_rungriddiagnostics    
    end function mnt_celllocator_rungriddiagnostics

  end interface

end module mnt_celllocator_capi_mod

