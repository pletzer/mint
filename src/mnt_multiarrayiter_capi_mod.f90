module mnt_multiarrayiter_capi_mod
  ! C function prototypes
  interface

    function mnt_multiarrayiter_new(obj, ndims, dims) &
                                 bind(C, name='mnt_multiarrayiter_new')
      ! Constructor
      ! @param obj instance of mntmultiarrayiter_t (opaque handle)
      ! @param ndims number of axes/dimensions
      ! @param dims dimension along each axis
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int), value            :: ndims
      integer(c_size_t), intent(in)    :: dims(*)
      integer(c_int)                   :: mnt_multiarrayiter_new
    end function mnt_multiarrayiter_new

    function mnt_multiarrayiter_del(obj) &
                                 bind(C, name='mnt_multiarrayiter_del')
      ! Destructor
      ! @param obj instance of mntmultiarrayiter_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_multiarrayiter_del
    end function mnt_multiarrayiter_del

    function mnt_multiarrayiter_begin(obj) &
                                 bind(C, name='mnt_multiarrayiter_begin')
      ! Set the iterator to the beginning
      ! @param obj instance of mntmultiarrayiter_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_multiarrayiter_begin
    end function mnt_multiarrayiter_begin

    function mnt_multiarrayiter_next(obj) &
                                 bind(C, name='mnt_multiarrayiter_next')
      ! Increment the iterator
      ! @param obj instance of mntmultiarrayiter_t (opaque handle)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_int)                   :: mnt_multiarrayiter_next
    end function mnt_multiarrayiter_next

    function mnt_multiarrayiter_getNumIters(obj, n) &
                                 bind(C, name='mnt_multiarrayiter_getNumIters')
      ! Get the number of iterations
      ! @param obj instance of mntmultiarrayiter_t (opaque handle)
      ! @param n number (output)
      ! @return 0 if successful
      use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_size_t), intent(out)   :: n
      integer(c_int)                   :: mnt_multiarrayiter_getNumIters
    end function mnt_multiarrayiter_getNumIters

    function mnt_multiarrayiter_getIndices(obj, indices) &
                                 bind(C, name='mnt_multiarrayiter_getIndices')
      ! Get the indices for this iteration
      ! @param obj instance of mntmultiarrayiter_t (opaque handle)
      ! @param indices array of indices (output)
      ! @return 0 if successful
      ! @note indices are zero based
      use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_ptr
      implicit none
      type(c_ptr), intent(inout)       :: obj ! void**
      integer(c_size_t), intent(out)   :: indices(*)
      integer(c_int)                   :: mnt_multiarrayiter_getIndices
    end function mnt_multiarrayiter_getIndices


  end interface

end module mnt_multiarrayiter_capi_mod

