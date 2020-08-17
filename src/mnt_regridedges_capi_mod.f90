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

    function mnt_regridedges_setSrcGridFlags(obj, fixLonAcrossDateline, averageLonAtPole) &
                                             bind(C, name='mnt_regridedges_setSrcGridFlags')
        ! Set source grid flags that control periodicity and regularization
        ! @param obj instance of mntRegridEdges_t (opaque handle)
        ! @param fixLonAcrossDateline set this to 1 to ensure that cells that are cut by the dateline don't wrap around
        ! @param averageLonAtPole set this to 1 to reset the longitude at the pole to match the average longitude of the cell, this
        !                         reduces the size of the cell in lon-lat space
        ! @return 0 if successful

        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj ! void**
        integer(c_int), value                    :: fixLonAcrossDateline
        integer(c_int), value                    :: averageLonAtPole
        integer(c_int)                           :: mnt_regridedges_setSrcGridFlags
    end function mnt_regridedges_setSrcGridFlags

    function mnt_regridedges_setDstGridFlags(obj, fixLonAcrossDateline, averageLonAtPole) &
                                             bind(C, name='mnt_regridedges_setDstGridFlags')
        ! Set destination grid flags that control periodicity and regularization
        ! @param obj instance of mntRegridEdges_t (opaque handle)
        ! @param fixLonAcrossDateline set this to 1 to ensure that cells that are cut by the dateline don't wrap around
        ! @param averageLonAtPole set this to 1 to reset the longitude at the pole to match the average longitude of the cell, this
        !                         reduces the size of the cell in lon-lat space
        ! @return 0 if successful

        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj ! void**
        integer(c_int), value                    :: fixLonAcrossDateline
        integer(c_int), value                    :: averageLonAtPole
        integer(c_int)                           :: mnt_regridedges_setDstGridFlags
    end function mnt_regridedges_setDstGridFlags

    function mnt_regridedges_initSliceIter(obj, &
                                           src_fort_filename, src_nFilenameLength, &
                                           dst_fort_filename, dst_nFilenameLength, &
                                           append, field_name, nFieldNameLength, numSlices) &
                                           bind(C, name='mnt_regridedges_initSliceIter')
        ! Read field metadata and initialize the slice iterator
        ! @param obj instance of mntRegridEdges_t (opaque handle)
        ! @param src_fort_filename file name (does not require termination character)
        ! @param src_nFilenameLength length of filename string (excluding '\0' if present)
        ! @param dst_fort_filename file name (does not require termination character)
        ! @param dst_nFilenameLength length of filename string (excluding '\0' if present)
        ! @param append set to 1 to append to exisiting netcdf file, 0 otherwise
        ! @param field_name name of the field
        ! @param nFieldNameLength length of field_name string (excluding '\0' if present)
        ! @param numSlices number of slices (output)
        ! @return 0 if successful
        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
        implicit none
        type(c_ptr), intent(inout)               :: obj ! void**
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

    function mnt_regridedges_loadSrcSlice(obj, data) bind(C, name='mnt_regridedges_loadSrcSlice')
        ! Load a source field slice from netcdf file
        ! @param obj instance of mntRegridEdges_t (opaque handle)
        ! @param data slice of edge data to be loaded
        ! @return 0 if successful
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj ! void**
        real(c_double), intent(out)              :: data(*)
        integer(c_int)                           :: mnt_regridedges_loadSrcSlice
    end function mnt_regridedges_loadSrcSlice

    function mnt_regridedges_dumpDstSlice(obj, data) bind(C, name='mnt_regridedges_dumpDstSlice')
        ! Dump a destination field slice to netcdf file
        ! @param obj instance of mntRegridEdges_t (opaque handle)
        ! @param data slice of edge data to be written
        ! @return 0 if successful
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj ! void**
        real(c_double), intent(out)              :: data(*)
        integer(c_int)                           :: mnt_regridedges_dumpDstSlice
    end function mnt_regridedges_dumpDstSlice

    function mnt_regridedges_nextSlice(obj) bind(C, name='mnt_regridedges_nextSlice')
        ! Incfrement the slice iterator
        ! @param obj instance of mntRegridEdges_t (opaque handle)
        ! @return 0 if successful
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        implicit none
        type(c_ptr), intent(inout)               :: obj ! void**
        integer(c_int)                           :: mnt_regridedges_nextSlice
    end function mnt_regridedges_nextSlice

    function mnt_regridedges_loadEdgeField(obj, filename, nFilenameLength, &
                                           field_name, nFieldNameLength, &
                                           ndata, data) &
                                           bind(C, name='mnt_regridedges_loadEdgeField')

       ! Load field from 2D UGRID file
       ! @param obj instance of mntRegridEdges_t (opaque handle)
       ! @param filename file name (does not require termination character)
       ! @param nFilenameLength length of filename string (excluding '\0' if present)
       ! @param field_name name of the field
       ! @param nFieldNameLength length of field_name string (excluding '\0' if present)
       ! @param ndata number of edges
       ! @param data array of size number of unique edges (output)
       ! @return 0 if successful

        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj ! void**
        character(kind=c_char), intent(in)       :: filename(*)
        integer(c_int), value                    :: nFilenameLength
        character(kind=c_char), intent(in)       :: field_name(*)
        integer(c_int), value                    :: nFieldNameLength
        integer(c_size_t), value                 :: ndata
        real(c_double), intent(out)              :: data(*)
        integer(c_int)                           :: mnt_regridedges_loadEdgeField
    end function mnt_regridedges_loadEdgeField

    function mnt_regridedges_dumpEdgeField(obj, filename, nFilenameLength, &
                                           field_name, nFieldNameLength, &
                                           ndata, data) &
                                           bind(C, name='mnt_regridedges_dumpEdgeField')
       ! Dump field to 2D UGRID file
       ! @param obj instance of mntRegridEdges_t (opaque handle)
       ! @param filename file name (does not require termination character)
       ! @param nFilenameLength length of filename string (excluding '\0' if present)
       ! @param field_name name of the field
       ! @param nFieldNameLength length of field_name string (excluding '\0' if present)
       ! @param ndata number of edges
       ! @param data array of size number of unique edges (output)
       ! @return 0 if successful

        use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char, c_double
        implicit none
        type(c_ptr), intent(inout)               :: obj ! void**
        character(kind=c_char), intent(in)       :: filename(*)
        integer(c_int), value                    :: nFilenameLength
        character(kind=c_char), intent(in)       :: field_name(*)
        integer(c_int), value                    :: nFieldNameLength
        integer(c_size_t), value                 :: ndata
        real(c_double), intent(in)               :: data(*)
        integer(c_int)                           :: mnt_regridedges_dumpEdgeField
    end function mnt_regridedges_dumpEdgeField

    function mnt_regridedges_loadSrcGrid(obj, filename, n) &
                                         bind(C, name='mnt_regridedges_loadSrcGrid')
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
      integer(c_int)                           :: mnt_regridedges_loadSrcGrid
    end function mnt_regridedges_loadSrcGrid

    function mnt_regridedges_loadDstGrid(obj, filename, n) &
                                         bind(C, name='mnt_regridedges_loadDstGrid')
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
      integer(c_int)                           :: mnt_regridedges_loadDstGrid
    end function mnt_regridedges_loadDstGrid

    function mnt_regridedges_build(obj, num_cells_per_bucket, debug) &
                                   bind(C, name='mnt_regridedges_build')
      ! Build locator object
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @param num_cells_per_bucket number of cells per bucket
      ! @param debug 0=no debug info, 1=print debug info, 2=save bad edges in VTK file
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: num_cells_per_bucket
      integer(c_int), value            :: debug
      integer(c_int)                   :: mnt_regridedges_build
    end function mnt_regridedges_build

    function mnt_regridedges_getNumSrcEdges(obj, n) &
                            bind(C, name='mnt_regridedges_getNumSrcEdges')
      ! Get the number of unique edges in the source grid
      ! @param n number (output)
      use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_size_t
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_size_t), intent(out)   :: n 
      integer(c_int)                   :: mnt_regridedges_getNumSrcEdges

    end function mnt_regridedges_getNumSrcEdges

    function mnt_regridedges_getNumDstEdges(obj, n) &
                            bind(C, name='mnt_regridedges_getNumDstEdges')
      ! Get the number of unique edges in the destination grid
      ! @param n number (output)
      use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_size_t
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_size_t), intent(out)   :: n 
      integer(c_int)                   :: mnt_regridedges_getNumDstEdges

    end function mnt_regridedges_getNumDstEdges

    function mnt_regridedges_loadWeights(obj, filename, n) & 
                                         bind(C, name='mnt_regridedges_loadWeights')
      ! Load interpolation weights from NetCDF file
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @param filename file name, \0 (char(0)) terminated
      ! @param n length of filename 
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_loadWeights
    end function mnt_regridedges_loadWeights

    function mnt_regridedges_dumpWeights(obj, filename, n) & 
                                         bind(C, name='mnt_regridedges_dumpWeights')
      ! Dump interpolation weights in NetCDF file
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @param filename file name, \0 (char(0)) terminated
      ! @param n length of filename 
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_int), value                    :: n
      integer(c_int)                           :: mnt_regridedges_dumpWeights
    end function mnt_regridedges_dumpWeights

    function mnt_regridedges_print(obj) & 
                                  bind(C, name='mnt_regridedges_print')
      ! Print interpolation weights
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      integer(c_int)                           :: mnt_regridedges_print   
    end function mnt_regridedges_print

    function mnt_regridedges_apply(obj, src_data, dst_data) &
                                   bind(C, name='mnt_regridedges_apply')
      ! Apply weights to unique cell Id edge field
      ! @param obj instance of mntRegridEdges_t (opaque handle)
      ! @param src_data source field defined for each unique edge Id
      ! @param dst_data destination field defined for each for each unique edge Id
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)               :: obj ! void**
      real(c_double), intent(in)               :: src_data(*)
      real(c_double), intent(out)              :: dst_data(*)
      integer(c_int)                           :: mnt_regridedges_apply
    end function mnt_regridedges_apply

  end interface

end module mnt_regridedges_capi_mod

