module write_data_mod

    implicit none

contains

    subroutine write_data(ncfile, varname, mode, u1)

        use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_ptr, c_size_t
        use mnt_ncfieldwrite_capi_mod 
        implicit none

        character(len=*), intent(in)   :: ncfile, varname, mode
        real(8), intent(inout)         :: u1(:, :, :)

        ! local variables
        integer                        :: ier
        type(c_ptr)                    :: handle
        character(len=32)              :: xName, yName, zName, unitStr, ndimsStr, fooStr
        character(len=32)              :: unitVal
        integer(c_int)                 :: ndimsVal
        real(c_double)                 :: fooVal
        integer(c_size_t)              :: start(3), nsizes(3), nx, ny, nz
        integer(c_int)                 :: append
 
        nx = size(u1, 1)
        ny = size(u1, 2)
        nz = size(u1, 3)

        append = 0
        if (mode(1:1) == 'a') append = 1

        ier = mnt_ncfieldwrite_new(handle, trim(ncfile), len_trim(ncfile), &
                                           trim(varname), len_trim(varname), &
                                           append)
        if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_new'

        if (append /= 1) then

            ! only need to define the dimensions and valiable when not in append mode

            ier = mnt_ncfieldwrite_setNumDims(handle, rank(u1))
            if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_setNumDims'

            ! axes
            xName = 'x'
            ier = mnt_ncfieldwrite_setDim(handle, 0, trim(xName), len_trim(xName), nx)
            if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_setDim of axis 0'

            yName = 'y'
            ier = mnt_ncfieldwrite_setDim(handle, 1, trim(yName), len_trim(yName), ny)
            if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_setDim of axis 1'

            zName = 'z'
            ier = mnt_ncfieldwrite_setDim(handle, 2, trim(zName), len_trim(zName), nz)
            if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_setDim of axis 2'

        endif

        ! write all the data
        ier = mnt_ncfieldwrite_data(handle, u1)
        if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_data ier = ', ier

        ! write a slice of the data. Note indexing starts at 0 and is column major (a la C)!
        start = [1, 1, 1]
        nsizes = [nx-1, ny-1, nz-1]
        u1 = -1
        ier = mnt_ncfieldwrite_dataSlice(handle, start, nsizes, u1)
        if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_data ier = ', ier


        ier = mnt_ncfieldwrite_del(handle)
        if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_del'

end subroutine write_data

end module write_data_mod

program test

    use write_data_mod
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_ptr, c_size_t

    implicit none

    integer(c_size_t), parameter   :: nx = 3, ny = 4, nz = 5
    integer(c_size_t)              :: i, j, k
    real(8)                        :: u1(nx, ny, nz)

    ! fill in
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                u1(i, j, k) = i-1 + (j-1)*nx + (k-1)*nx*ny
            enddo
        enddo
    enddo

    call write_data('test_ncfieldwrite.nc', 'u1', 'new', u1)

    ! uodate the data
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                u1(i, j, k) = 1000 + i-1 + (j-1)*nx + (k-1)*nx*ny
            enddo
        enddo
    enddo
    call write_data('test_ncfieldwrite.nc', 'u1', 'append', u1)


end program