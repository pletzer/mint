module mnt_ncfieldread_capi_mod
  ! C function prototypes
  interface

    function mnt_ncfieldread_new(obj, ncid, varid) &
                                 bind(C, name='mnt_ncfieldread_new')
      ! Constructor
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)             :: obj ! void**
      integer(c_int), value                  :: ncid
      integer(c_int), value                  :: varid
      integer(c_int)                         :: mnt_ncfieldread_new
    end function mnt_ncfieldread_new

    function mnt_ncfieldread_del(obj) &
                                 bind(C, name='mnt_ncfieldread_del')
      ! Destructor
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_ncfieldread_del
    end function mnt_ncfieldread_del

    
    function mnt_ncfieldread_getNumDims(obj, ndims) &
                                        bind(C, name='mnt_ncfieldread_getNumDims')
      ! Get the number of dimenions of this netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param ndims number
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: ndims
      integer(c_int)                   :: mnt_ncfieldread_getNumDims
    end function mnt_ncfieldread_getNumDims

    function mnt_ncfieldread_getDimName(obj, iAxis, dimName, dimNameLen) &
                                        bind(C, name='mnt_ncfieldread_getDimName')
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
      integer(c_int)                         :: mnt_ncfieldread_getDimName
    end function mnt_ncfieldread_getDimName

    function mnt_ncfieldread_getDim(obj, iAxis, dim) &
                                        bind(C, name='mnt_ncfieldread_getDim')
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
      integer(c_int)                   :: mnt_ncfieldread_getDim
    end function mnt_ncfieldread_getDim

    function mnt_ncfieldread_data(obj, data) &
                                      bind(C, name='mnt_ncfieldread_data')
      ! Read the entire netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param data array, will be filled in
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      real(c_double), intent(out)      :: data(*)
      integer(c_int)                   :: mnt_ncfieldread_data
    end function mnt_ncfieldread_data

    function mnt_ncfieldread_dataSlice(obj, startInds0, counts, data) &
                                       bind(C, name='mnt_ncfieldread_dataSlice')
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
      integer(c_int)                   :: mnt_ncfieldread_dataSlice
    end function mnt_ncfieldread_dataSlice   

  end interface

end module mnt_ncfieldread_capi_mod

