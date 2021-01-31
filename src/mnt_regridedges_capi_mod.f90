module mnt_regridedges_capi_mod
  !> C function prototypes
  interface

    !> Constructor
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @return 0 if successful
    function mnt_regridedges_new(obj) &
                                 bind(C, name='mnt_regridedges_new')
        use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
        implicit none
        type(c_ptr), intent(inout)       :: obj !> void**
        integer(c_int)                   :: mnt_regridedges_new
    end function mnt_regridedges_new

    !> Destructor
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @return 0 if successful
    function mnt_regridedges_del(obj) &
                                 bind(C, name='mnt_regridedges_del')
        use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
        implicit none
        type(c_ptr), intent(inout)       :: obj !> void**
        integer(c_int)                   :: mnt_regridedges_del
    end function mnt_regridedges_del

    !> Set the source grid flags that control periodicity and regularization
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param fixLonAcrossDateline set this to 1 to ensure that cells that are cut by the dateline don't wrap around
    !> @param averageLonAtPole set this to 1 to reset the longitude at the pole to match the average longitude of the cell, this
    !>                         reduces the size of the cell in lon-lat space
    !> @return 0 if successful
    function mnt_regridedges_setSrcGridFlags(obj, fixLonAcrossDateline, averageLonAtPole) &
                                             bind(C, name='mnt_regridedges_setSrcGridFlags')
        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj !> void**
        integer(c_int), value                    :: fixLonAcrossDateline
        integer(c_int), value                    :: averageLonAtPole
        integer(c_int)                           :: mnt_regridedges_setSrcGridFlags
    end function mnt_regridedges_setSrcGridFlags

    !> Set the destination grid flags that control periodicity and regularization
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param fixLonAcrossDateline set this to 1 to ensure that cells that are cut by the dateline don't wrap around
    !> @param averageLonAtPole set this to 1 to reset the longitude at the pole to match the average longitude of the cell, this
    !>                         reduces the size of the cell in lon-lat space
    !> @return 0 if successful
    function mnt_regridedges_setDstGridFlags(obj, fixLonAcrossDateline, averageLonAtPole) &
                                             bind(C, name='mnt_regridedges_setDstGridFlags')
        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj !> void**
        integer(c_int), value                    :: fixLonAcrossDateline
        integer(c_int), value                    :: averageLonAtPole
        integer(c_int)                           :: mnt_regridedges_setDstGridFlags
    end function mnt_regridedges_setDstGridFlags

    !> Initialize the slice iterator
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param src_fort_filename file name (does not require termination character)
    !> @param src_nFilenameLength length of filename string (excluding '\0' if present)
    !> @param dst_fort_filename file name (does not require termination character)
    !> @param dst_nFilenameLength length of filename string (excluding '\0' if present)
    !> @param append set to 1 to append to exisiting netcdf file, 0 otherwise
    !> @param field_name name of the field
    !> @param nFieldNameLength length of field_name string (excluding '\0' if present)
    !> @param numSlices number of slices (output)
    !> @return 0 if successful
    function mnt_regridedges_initSliceIter(obj, &
                                           src_fort_filename, src_nFilenameLength, &
                                           dst_fort_filename, dst_nFilenameLength, &
                                           append, field_name, nFieldNameLength, numSlices) &
                                           bind(C, name='mnt_regridedges_initSliceIter')
        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
        implicit none
        type(c_ptr), intent(inout)               :: obj !> void**
        character(kind=c_char), intent(in)       :: src_fort_filename(*)
        integer(c_int), value                    :: src_nFilenameLength
        character(kind=c_char), intent(in)       :: dst_fort_filename(*)
        integer(c_int), value                    :: dst_nFilenameLength
        integer(c_int), value                    :: append
        character(kind=c_char), intent(in)       :: field_name(*)
        integer(c_int), value                    :: nFieldNameLength
        integer(c_size_t), intent(out)           :: numSlices
        integer(c_int)                           :: mnt_regridedges_initSliceIter
    end function mnt_regridedges_initSliceIter

    !> Load/read a source field slice from a netcdf file
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param data slice of edge data to be loaded
    !> @return 0 if successful
    function mnt_regridedges_loadSrcSlice(obj, data) bind(C, name='mnt_regridedges_loadSrcSlice')
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj !> void**
        real(c_double), intent(out)              :: data(*)
        integer(c_int)                           :: mnt_regridedges_loadSrcSlice
    end function mnt_regridedges_loadSrcSlice

    !> Dump/write a destination field slice to a netcdf file
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param data slice of edge data to be written
    !> @return 0 if successful
    function mnt_regridedges_dumpDstSlice(obj, data) bind(C, name='mnt_regridedges_dumpDstSlice')
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj !> void**
        real(c_double), intent(out)              :: data(*)
        integer(c_int)                           :: mnt_regridedges_dumpDstSlice
    end function mnt_regridedges_dumpDstSlice

    !> Increment the slice iterator
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @return 0 if successful
    function mnt_regridedges_nextSlice(obj) bind(C, name='mnt_regridedges_nextSlice')
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none
        type(c_ptr), intent(inout)               :: obj !> void**
        integer(c_int)                           :: mnt_regridedges_nextSlice
    end function mnt_regridedges_nextSlice

   !> Load a field from 2D UGRID file
   !> @param obj instance of mntRegridEdges_t (opaque handle)
   !> @param filename file name (does not require termination character)
   !> @param nFilenameLength length of filename string (excluding '\0' if present)
   !> @param field_name name of the field
   !> @param nFieldNameLength length of field_name string (excluding '\0' if present)
   !> @param ndata number of edges
   !> @param data array of size number of unique edges (output)
   !> @return 0 if successful
    function mnt_regridedges_loadEdgeField(obj, filename, nFilenameLength, &
                                           field_name, nFieldNameLength, &
                                           ndata, data) &
                                           bind(C, name='mnt_regridedges_loadEdgeField')
        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj !> void**
        character(kind=c_char), intent(in)       :: filename(*)
        integer(c_int), value                    :: nFilenameLength
        character(kind=c_char), intent(in)       :: field_name(*)
        integer(c_int), value                    :: nFieldNameLength
        integer(c_size_t), value                 :: ndata
        real(c_double), intent(out)              :: data(*)
        integer(c_int)                           :: mnt_regridedges_loadEdgeField
    end function mnt_regridedges_loadEdgeField

   !> Dump a field to a 2D UGRID file
   !> @param obj instance of mntRegridEdges_t (opaque handle)
   !> @param filename file name (does not require termination character)
   !> @param nFilenameLength length of filename string (excluding '\0' if present)
   !> @param field_name name of the field
   !> @param nFieldNameLength length of field_name string (excluding '\0' if present)
   !> @param ndata number of edges
   !> @param data array of size number of unique edges (output)
   !> @return 0 if successful
    function mnt_regridedges_dumpEdgeField(obj, filename, nFilenameLength, &
                                           field_name, nFieldNameLength, &
                                           ndata, data) &
                                           bind(C, name='mnt_regridedges_dumpEdgeField')
        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj !> void**
        character(kind=c_char), intent(in)       :: filename(*)
        integer(c_int), value                    :: nFilenameLength
        character(kind=c_char), intent(in)       :: field_name(*)
        integer(c_int), value                    :: nFieldNameLength
        integer(c_size_t), value                 :: ndata
        real(c_double), intent(in)               :: data(*)
        integer(c_int)                           :: mnt_regridedges_dumpEdgeField
    end function mnt_regridedges_dumpEdgeField

    !> Load a source grid from a 2d UGRID file 
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param filename file name
    !> @param n length of string filename
    !> @return 0 if sucecssful
    !> @note must invoke constructor prior to this call
    function mnt_regridedges_loadSrcGrid(obj, filename, n) &
                                         bind(C, name='mnt_regridedges_loadSrcGrid')
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj !> void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_loadSrcGrid
    end function mnt_regridedges_loadSrcGrid

    !> Load a destination grid from a 2d UGRID file 
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param filename file name
    !> @param n length of string filename
    !> @return 0 if sucecssful
    !> @note must invoke constructor prior to this call
    function mnt_regridedges_loadDstGrid(obj, filename, n) &
                                         bind(C, name='mnt_regridedges_loadDstGrid')
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj !> void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_loadDstGrid
    end function mnt_regridedges_loadDstGrid

    !> Build the locator object
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param num_cells_per_bucket number of cells per bucket
    !> @param periodX period in longitudes (0 if non-periodic)
    !> @param debug 0=no debug info, 1=print debug info, 2=save bad edges in VTK file
    !> @return 0 if successful
    function mnt_regridedges_build(obj, num_cells_per_bucket, periodX, debug) &
                                   bind(C, name='mnt_regridedges_build')
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj !> void**
      integer(c_int), value            :: num_cells_per_bucket
      real(c_double), value            :: periodX
      integer(c_int), value            :: debug
      integer(c_int)                   :: mnt_regridedges_build
    end function mnt_regridedges_build

    !> Get the number of unique edges in the source grid
    !> @param n number (output)
    function mnt_regridedges_getNumSrcEdges(obj, n) &
                            bind(C, name='mnt_regridedges_getNumSrcEdges')
      use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_size_t
      type(c_ptr), intent(inout)       :: obj !> void**
      integer(c_size_t), intent(out)   :: n 
      integer(c_int)                   :: mnt_regridedges_getNumSrcEdges

    end function mnt_regridedges_getNumSrcEdges

    !> Get the number of unique edges in the destination grid
    !> @param n number (output)
    function mnt_regridedges_getNumDstEdges(obj, n) &
                            bind(C, name='mnt_regridedges_getNumDstEdges')
      use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_size_t
      type(c_ptr), intent(inout)       :: obj !> void**
      integer(c_size_t), intent(out)   :: n 
      integer(c_int)                   :: mnt_regridedges_getNumDstEdges

    end function mnt_regridedges_getNumDstEdges

    !> Load the interpolation weights from a netCDF file
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param filename file name, \0 (char(0)) terminated
    !> @param n length of filename 
    !> @return 0 if successful
    function mnt_regridedges_loadWeights(obj, filename, n) & 
                                         bind(C, name='mnt_regridedges_loadWeights')
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj !> void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_loadWeights
    end function mnt_regridedges_loadWeights

    !> Dump/write the interpolation weights to a NetCDF file
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param filename file name, \0 (char(0)) terminated
    !> @param n length of filename 
    !> @return 0 if successful
    function mnt_regridedges_dumpWeights(obj, filename, n) & 
                                         bind(C, name='mnt_regridedges_dumpWeights')
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj !> void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_dumpWeights
    end function mnt_regridedges_dumpWeights

    !> Print the interpolation weights
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @return 0 if successful
    function mnt_regridedges_print(obj) & 
                                  bind(C, name='mnt_regridedges_print')
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj !> void**
      integer(c_int)                           :: mnt_regridedges_print   
    end function mnt_regridedges_print

    !> Apply the weights to the cell by cell edge field
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param src_data edge field defined cell by cell on the source grid (array of size 4*num_src_cells)
    !> @param dst_data edge field defined cell by cell on the destination grid (array of size 4*num_dst_cells)
    !> @return 0 if successful
    function mnt_regridedges_applyCellEdge(obj, src_data, dst_data) &
                                           bind(C, name='mnt_regridedges_applyCellEdge')
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj !> void**
      real(c_double), intent(in)               :: src_data(*)
      real(c_double), intent(out)              :: dst_data(*)
      integer(c_int)                           :: mnt_regridedges_applyCellEdge
    end function mnt_regridedges_applyCellEdge

    !> Apply the weights to an edge field with unique edge Ids
    !> @param obj instance of mntRegridEdges_t (opaque handle)
    !> @param src_data source field defined for each unique edge Id
    !> @param dst_data destination field defined for each for each unique edge Id
    !> @return 0 if successful
    function mnt_regridedges_apply(obj, src_data, dst_data) &
                                   bind(C, name='mnt_regridedges_apply')
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj !> void**
      real(c_double), intent(in)               :: src_data(*)
      real(c_double), intent(out)              :: dst_data(*)
      integer(c_int)                           :: mnt_regridedges_apply
    end function mnt_regridedges_apply

  end interface

end module mnt_regridedges_capi_mod

