module mnt_ncdimensions_capi_mod
  ! C function prototypes
  interface

    function mnt_ncdimensions_new(obj, ncid, varid) &
                                 bind(C, name='mnt_ncdimensions_new')
      ! Constructor
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)             :: obj ! void**
      integer(c_int), value                  :: ncid
      integer(c_int), value                  :: varid
      integer(c_int)                         :: mnt_ncdimensions_new
    end function mnt_ncdimensions_new

    function mnt_ncdimensions_del(obj) &
                                 bind(C, name='mnt_ncdimensions_del')
      ! Destructor
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_ncdimensions_del
    end function mnt_ncdimensions_del

    function mnt_ncdimensions_read(obj, ncid, varid) &
                                        bind(C, name='mnt_ncdimensions_read')
      ! Read the dimensions of the netcdf file variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param ncid netcdf file Id
      ! @param varid netcdf variable Id
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: ncid
      integer(c_int), value            :: varid
      integer(c_int)                   :: mnt_ncdimensions_read
    end function mnt_ncdimensions_read
    
    function mnt_ncdimensions_getNumDims(obj, ndims) &
                                        bind(C, name='mnt_ncdimensions_getNumDims')
      ! Get the number of dimenions of this netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param ndims number
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), intent(out)      :: ndims
      integer(c_int)                   :: mnt_ncdimensions_getNumDims
    end function mnt_ncdimensions_getNumDims

    function mnt_ncdimensions_get(obj, i, size) &
                                        bind(C, name='mnt_ncdimensions_get')
      ! Get the number of dimenions of this netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param ndims number
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_size_t
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: i
      integer(c_size_t), intent(out)   :: size
      integer(c_int)                   :: mnt_ncdimensions_get
    end function mnt_ncdimensions_get



    function mnt_ncdimensions_getDimName(obj, iAxis, dimName, dimNameLen) &
                                        bind(C, name='mnt_ncdimensions_getDimName')
      ! Get the dimension name of this netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param iAxis dimension/axis number (0 based)
      ! @param dimName dimension name
      ! @param dimNameLen max number of characters in dimName 
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)             :: obj ! void**
      integer(c_int), value                  :: iAxis
      character(kind=c_char), intent(out)    :: dimName(*)
      integer(c_int), value                  :: dimNameLen
      integer(c_int)                         :: mnt_ncdimensions_getDimName
    end function mnt_ncdimensions_getDimName

    function mnt_ncdimensions_getDim(obj, iAxis, dim) &
                                        bind(C, name='mnt_ncdimensions_getDim')
      ! Get the dimension of this netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param iAxis dimension/axis number (0 based)
      ! @param dim dimension size
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_size_t
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: iAxis
      integer(c_size_t), intent(out)   :: dim
      integer(c_int)                   :: mnt_ncdimensions_getDim
    end function mnt_ncdimensions_getDim

    function mnt_ncdimensions_data(obj, data) &
                                      bind(C, name='mnt_ncdimensions_data')
      ! Read the entire netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param data array, will be filled in
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      real(c_double), intent(out)      :: data(*)
      integer(c_int)                   :: mnt_ncdimensions_data
    end function mnt_ncdimensions_data

    function mnt_ncdimensions_dataSlice(obj, startInds0, counts, data) &
                                       bind(C, name='mnt_ncdimensions_dataSlice')
      ! Read a slice of the netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param startInds0 starting indices (0 based and assuming C ordering)
      ! @param counts number of values for each axis/dimension (assumed C ordering)
      ! @param data array, will be filled in
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double, c_size_t
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_size_t), intent(in)    :: startInds0(*)
      integer(c_size_t), intent(in)    :: counts(*)
      real(c_double), intent(out)      :: data(*)
      integer(c_int)                   :: mnt_ncdimensions_dataSlice
    end function mnt_ncdimensions_dataSlice   

  end interface

end module mnt_ncdimensions_capi_mod

