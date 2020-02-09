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

  end interface

end module mnt_ncfieldread_capi_mod

