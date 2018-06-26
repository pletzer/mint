module constants_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64, int32, int64

  implicit none

  integer, parameter :: r_def     = real64 !< Default real kind for application.
  integer, parameter :: r_single  = real32 !< Default single precision real kind for application.
  integer, parameter :: r_double  = real64 !< Default double precision real kind for application.
  integer, parameter :: i_def        = int32       !< Default integer kind for application.
  integer, parameter :: i_addrss_def = int64       !< For integers that hold addresses.

end module constants_mod


module test_mod
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

    function mnt_celllocator_setpoints(obj, nverts_per_cell, ncells, verts) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      integer(c_int), intent(in)               :: nverts_per_cell
      integer(c_size_t), intent(in)            :: ncells
      real(c_double), intent(in)               :: verts ! const double*
      integer(c_int)                           :: mnt_celllocator_setpoints
    end function mnt_celllocator_setpoints

    function mnt_celllocator_build(obj) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout) :: obj
      integer(c_int)                   :: mnt_celllocator_build
    end function mnt_celllocator_build

    function mnt_celllocator_find(obj, point, cell_id, pcoords) bind(C)
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double
      integer(c_size_t), intent(inout)         :: obj
      real(c_double), intent(in)               :: point   ! const double*
      integer(c_size_t), intent(out)           :: cell_id
      real(c_double), intent(out)              :: pcoords ! double*
      integer(c_int)                           :: mnt_celllocator_find
    end function mnt_celllocator_find
  end interface
end module test_mod


program test
    use constants_mod
    use test_mod
    implicit none
    integer(c_size_t) :: cloc
    integer(c_size_t) :: num_cells = 1
    integer(idef)     :: num_verts_per_cell = 8
    integer(c_size_t) :: num_points = num_cells*num_verts_per_cell
    real(r_def)       :: verts(num_points * 3)
    real(r_def)       :: target_point(3) = [0.2, 0.3, 0.4]
    real(r_def)       :: pcoords(3)
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
    ier = mnt_celllocator_setpoints(cloc, num_verts_per_cell, num_cells, verts)
    ier = mnt_celllocator_build(cloc)
    ier = mnt_celllocator_find(cloc, target_point, cell_id, pcoords)
    ier = mnt_celllocator_del(cloc)

    print *,'cell_id: ', cell_id, ' pcoords = ', pcoords

end program 