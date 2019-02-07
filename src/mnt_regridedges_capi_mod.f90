module mnt_regridedges_capi_mod
  ! C function prototypes
  interface

    function mnt_regridedges_new(obj) &
                                 bind(C, name='mnt_regridedges_new')
      ! Constructor
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_regridedges_new
    end function mnt_regridedges_new

    function mnt_regridedges_del(obj) &
                                 bind(C, name='mnt_regridedges_del')
      ! Destructor
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_regridedges_del
    end function mnt_regridedges_del

    function mnt_regridedges_loadsrc(obj, filename, n) &
                                     bind(C, name='mnt_regridedges_loadSrc')
      ! Load source grid from 2d UGRID file 
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @param filename file name
      ! @param n length of string filename
      ! @return 0 if sucecssful
      ! @note must invoke constructor prior to this call
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_loadsrc
    end function mnt_regridedges_loadsrc

    function mnt_regridedges_loaddst(obj, filename, n) &
                                     bind(C, name='mnt_regridedges_loadDst')
      ! Load destination grid from 2d UGRID file 
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @param filename file name
      ! @param n length of string filename
      ! @return 0 if sucecssful
      ! @note must invoke constructor prior to this call
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_loaddst
    end function mnt_regridedges_loaddst

    function mnt_regridedges_build(obj, num_cells_per_bucket) &
                                   bind(C, name='mnt_regridedges_build')
      ! Build locator object
      ! @param obj instance of mntregridedges_t (opaque handle)
      ! @param num_cells_per_bucket number of cells per bucket
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: num_cells_per_bucket
      integer(c_int)                   :: mnt_regridedges_build
    end function mnt_regridedges_build

    function mnt_regridedges_dump(obj, filename, n) & 
                                  bind(C, name='mnt_regridedges_dump')
      ! Dump interpolation weights in NetCDF file
      ! @param obj instance of mntregridedges_t (opaque handle)
      ! @param filename file name, \0 (char(0)) terminated
      ! @param n length of filename 
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_dump    
    end function mnt_regridedges_dump

    function mnt_regridedges_print(obj) & 
                                  bind(C, name='mnt_regridedges_print')
      ! Print interpolation weights
      ! @param obj instance of mntregridedges_t (1opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_int)                           :: mnt_regridedges_print   
    end function mnt_regridedges_print
  end interface

end module mnt_regridedges_capi_mod

