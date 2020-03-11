module mnt_ncfieldwrite_capi_mod
  ! C function prototypes
  interface

    function mnt_ncfieldwrite_new(obj, fileName, fileNameLen, &
    	                            varName, varNameLen, append) &
                                  bind(C, name='mnt_ncfieldwrite_new')
      ! Constructor
      ! @param obj instance of mntNcFieldwrite_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)             :: obj ! void**
      character(kind=c_char), intent(in)     :: fileName(*)
      integer(c_int), value                  :: fileNameLen
      character(kind=c_char), intent(in)     :: varName(*)
      integer(c_int), value                  :: varNameLen
      integer(c_int), value                  :: append
      integer(c_int)                         :: mnt_ncfieldwrite_new
    end function mnt_ncfieldwrite_new

    function mnt_ncfieldwrite_del(obj) &
                                 bind(C, name='mnt_ncfieldwrite_del')
      ! Destructor
      ! @param obj instance of mntNcFieldwrite_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_ncfieldwrite_del
    end function mnt_ncfieldwrite_del

    
    function mnt_ncfieldwrite_setNumDims(obj, ndims) &
                                         bind(C, name='mnt_ncfieldwrite_setNumDims')
      ! Set the number of dimenions of this netcdf variable
      ! @param obj instance of mntNcFieldwrite_t (opaque handle)
      ! @param ndims number
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: ndims
      integer(c_int)                   :: mnt_ncfieldwrite_setNumDims
    end function mnt_ncfieldwrite_setNumDims

    function mnt_ncfieldwrite_setDim(obj, iAxis, dimName, dimNameLen, dim) &
                                        bind(C, name='mnt_ncfieldwrite_setDim')
      ! set the dimension name of this netcdf variable
      ! @param obj instance of mntNcFieldwrite_t (opaque handle)
      ! @param iAxis dimension/axis number (0 based)
      ! @param dimName dimension name
      ! @param dimNameLen max number of characters in dimName
      ! @param dim dimension size
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char, c_size_t
      implicit none
      type(c_ptr), intent(inout)             :: obj ! void**
      integer(c_int), value                  :: iAxis
      character(kind=c_char), intent(in)     :: dimName(*)
      integer(c_int), value                  :: dimNameLen
      integer(c_size_t), value               :: dim
      integer(c_int)                         :: mnt_ncfieldwrite_setDim
    end function mnt_ncfieldwrite_setDim

    function mnt_ncfieldwrite_data(obj, data) &
                                   bind(C, name='mnt_ncfieldwrite_data')
      ! write the entire netcdf variable
      ! @param obj instance of mntNcFieldwrite_t (opaque handle)
      ! @param data array, will be filled in
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      real(c_double), intent(in)       :: data(*)
      integer(c_int)                   :: mnt_ncfieldwrite_data
    end function mnt_ncfieldwrite_data

    function mnt_ncfieldwrite_dataSlice(obj, startInds0, counts, data) &
                                        bind(C, name='mnt_ncfieldwrite_dataSlice')
      ! write a slice of the netcdf variable
      ! @param obj instance of mntNcFieldwrite_t (opaque handle)
      ! @param startInds0 starting indices (0 based and assuming C ordering)
      ! @param counts number of values for each axis/dimension (assumed C ordering)
      ! @param data array, will be filled in
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double, c_size_t
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_size_t), intent(in)    :: startInds0(*)
      integer(c_size_t), intent(in)    :: counts(*)
      real(c_double), intent(in)       :: data(*)
      integer(c_int)                   :: mnt_ncfieldwrite_dataSlice
    end function mnt_ncfieldwrite_dataSlice   

  end interface

end module mnt_ncfieldwrite_capi_mod

