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

    function mnt_regridedges_loadSrc(obj, filename, n) &
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
    end function mnt_regridedges_loadSrc

    function mnt_regridedges_loadDst(obj, filename, n) &
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
    end function mnt_regridedges_loadDst

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

    function mnt_regridedges_load(obj, filename, n) & 
                                  bind(C, name='mnt_regridedges_load')
      ! Load interpolation weights from NetCDF file
      ! @param obj instance of mntregridedges_t (opaque handle)
      ! @param filename file name, \0 (char(0)) terminated
      ! @param n length of filename 
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_load  
    end function mnt_regridedges_load

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
      ! @param obj instance of mntregridedges_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_int)                           :: mnt_regridedges_print   
    end function mnt_regridedges_print

    function mnt_regridedges_applyWeightsToCellEdgeField(obj, src_data, &
                                                         dst_data) &
                    bind(C, name='mnt_regridedges_applyWeightsToCellEdgeField')
      ! Apply weights to cell by cell edge field
      ! @param obj instance of mntregridedges_t (opaque handle)
      ! @param src_data edge field defined cell by cell on the source grid (array of size 4*num_src_cells)
      ! @param dst_data edge field defined cell by cell on the destination grid (array of size 4*num_dst_cells)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      real(c_double), intent(in)               :: src_data(*)
      real(c_double), intent(out)              :: dst_data(*)
      integer(c_int)                           :: mnt_regridedges_applyWeightsToCellEdgeField 
    end	function mnt_regridedges_applyWeightsToCellEdgeField

    function mnt_regridedges_applyWeightsToEdgeIdField(obj, src_data, &
                                                       num_dst_edges, dst_data) &
                    bind(C, name='mnt_regridedges_applyWeightsToEdgeIdField')
      ! Apply weights to unique cell Id edge field
      ! @param obj instance of mntregridedges_t (opaque handle)
      ! @param src_data source field defined for each unique edge Id
      ! @param num_dst_edges number of edges on destination grid
      ! @param dst_data destination field defined for each for each unique edge Id
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      real(c_double), intent(in)               :: src_data(*)
      integer(c_size_t), value                 :: num_dst_edges
      real(c_double), intent(out)              :: dst_data(*)
      integer(c_int)                           :: mnt_regridedges_applyWeightsToEdgeIdField 
    end	function mnt_regridedges_applyWeightsToEdgeIdField

  end interface

end module mnt_regridedges_capi_mod

