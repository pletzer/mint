module mnt_ncdimensions_capi_mod
  ! C function prototypes
  interface

    function mnt_ncdimensions_new(obj) &
                                 bind(C, name='mnt_ncdimensions_new')
      ! Constructor
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)             :: obj ! void**
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

    function mnt_ncdimensions_print(obj) &
                                       bind(C, name='mnt_ncdimensions_print')
      ! Print the object
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_ncdimensions_print
    end function mnt_ncdimensions_print


  end interface

end module mnt_ncdimensions_capi_mod

