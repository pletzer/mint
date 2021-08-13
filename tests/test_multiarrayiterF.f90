subroutine test(ndims)
    use mnt_multiarrayiter_capi_mod
    use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_ptr
    implicit none
    integer, intent(in) :: ndims

    integer           :: ier, i
    integer(c_size_t) :: nelems
    integer(c_size_t), allocatable :: dims(:)
    integer(c_size_t), allocatable  :: indices(:)
    type(c_ptr)                    :: handle

    allocate(dims(ndims), indices(ndims))

    do i = 1, ndims
        dims(i) = 1 + i
    enddo
    ier = mnt_multiarrayiter_new(handle, ndims, dims)
    if (ier /= 0) print*,'ERROR'

    ier = mnt_multiarrayiter_getNumIters(handle, nelems)
    if (ier /= 0) print*,'ERROR'    

    ier = mnt_multiarrayiter_begin(handle)
    if (ier /= 0) print*,'ERROR'

    do i = 1, nelems

        ier = mnt_multiarrayiter_getIndices(handle, indices)
        if (ier /= 0) print*,'ERROR'

        print*,' ndims = ', ndims, ' iteration: ', i, ' indices = ', indices

        ier = mnt_multiarrayiter_next(handle)
        if (ier /= 0) print*,'ERROR'
    enddo


    ier = mnt_multiarrayiter_del(handle)
    if (ier /= 0) print*,'ERROR'

    deallocate(dims, indices)

end subroutine test

program test_multiarrayiterF

    call test(1)
    call test(2)
    call test(3)

end program