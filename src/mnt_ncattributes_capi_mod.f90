module mnt_ncattributes_capi_mod
  ! C function prototypes
  interface

    function mnt_ncattributes_new(obj) &
                                  bind(C, name='mnt_ncattributes_new')
      ! Constructor
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char
      implicit none
      type(c_ptr), intent(inout)             :: obj ! void**
      integer(c_int)                         :: mnt_ncattributes_new
    end function mnt_ncattributes_new

    function mnt_ncattributes_del(obj) &
                                  bind(C, name='mnt_ncattributes_del')
      ! Destructor
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_ncattributes_del
    end function mnt_ncattributes_del

    function mnt_ncattributes_read(obj, ncid, varid) &
                                        bind(C, name='mnt_ncattributes_read')
      ! Read the attributes of a variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param ncid netcdf file Id
      ! @param varid netcdf variable Id
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: ncid
      integer(c_int), value            :: varid
      integer(c_int)                   :: mnt_ncattributes_read
    end function mnt_ncattributes_read
    
    function mnt_ncattributes_write(obj, ncid, varid) &
                                        bind(C, name='mnt_ncattributes_write')
      ! Write the attributes of a variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param ncid netcdf file Id
      ! @param varid netcdf variable Id
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: ncid
      integer(c_int), value            :: varid
      integer(c_int)                   :: mnt_ncattributes_write
    end function mnt_ncattributes_write

    function mnt_ncattributes_isIntensive(obj) &
                                        bind(C, name='mnt_ncattributes_isIntensive')
      ! Get the number of dimenions of this netcdf variable
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @param ndims number
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_ncattributes_isIntensive
    end function mnt_ncattributes_isIntensive

    function mnt_ncattributes_print(obj) &
                                       bind(C, name='mnt_ncattributes_print')
      ! Print the object
      ! @param obj instance of mntNcFieldRead_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_ncattributes_print
    end function mnt_ncattributes_print

  end interface

end module mnt_ncattributes_capi_mod

