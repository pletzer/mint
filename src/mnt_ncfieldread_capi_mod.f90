module mnt_ncfieldread_capi_mod
  ! C function prototypes
  interface

    function mnt_ncfieldread_new(obj, fileName, fileNameLen, varName, varNameLen) &
                                 bind(C, name='mnt_ncfieldread_new')
      ! Constructor
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      character(len=1), intent(in)     :: fileName
      integer(c_int), value            :: fileNameLen
      character(len=1), intent(in)     :: varName
      integer(c_int), value            :: varNameLen
      integer(c_int)                   :: mnt_ncfieldread_new
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
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: ndims
      integer(c_int)                   :: mnt_ncfieldread_getNumDims
    end function mnt_ncfieldread_getNumDims

    function mnt_ncfieldread_getDimName(obj, iAxis, dimName, dimNameLen) &
                                        bind(C, name='mnt_ncfieldread_getDimName')
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), intent(in)       :: iAxis
      character(len=1), intent(out)    :: dimName
      integer(c_int), value            :: dimNameLen
      integer(c_int)                   :: mnt_ncfieldread_getDimName
    end function mnt_ncfieldread_getDimName

    function mnt_ncfieldread_getDim(obj, iAxis, dim) &
                                        bind(C, name='mnt_ncfieldread_getDim')
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: iAxis
      integer(c_int), intent(out)      :: dim
      integer(c_int)                   :: mnt_ncfieldread_getDim
    end function mnt_ncfieldread_getDim

    function mnt_ncfieldread_readData(obj, data) &
                                        bind(C, name='mnt_ncfieldread_readData')
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      real(c_double), intent(out)      :: data
      integer(c_int)                   :: mnt_ncfieldread_readData
    end function mnt_ncfieldread_readData

    function mnt_ncfieldread_readDataSlice(obj, startInds0, counts, data) &
                                        bind(C, name='mnt_ncfieldread_readDataSlice')
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_double, c_size_t
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_size_t), intent(in)    :: startInds0(*)
      integer(c_size_t), intent(in)    :: counts(*)
      real(c_double), intent(out)      :: data
      integer(c_int)                   :: mnt_ncfieldread_readDataSlice
    end function mnt_ncfieldread_readDataSlice   

  end interface

end module mnt_ncfieldread_capi_mod

