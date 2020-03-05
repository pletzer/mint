program test
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_ptr, c_size_t
    use mnt_ncfieldwrite_capi_mod 
    implicit none
    integer                        :: ier
    type(c_ptr)                    :: handle
    character(len=512)             :: ncfile = 'test_ncfieldwrite.nc'
    character(len=32)              :: varname = 'u1'
    character(len=32)              :: xName, yName, zName, unitStr, ndimsStr, fooStr
    integer(c_size_t)              :: i, j, k
    integer(c_size_t), parameter   :: nx = 3, ny = 4, nz = 5
    real(8)                        :: u1(nx, ny, nz)
    character(len=32)              :: unitVal
    integer(c_int)                 :: ndimsVal
    real(c_double)                 :: fooVal
    integer(c_size_t)              :: start(3), nsizes(3)

    ier = mnt_ncfieldwrite_new(handle, trim(ncfile), len_trim(ncfile), &
                                       trim(varname), len_trim(varname))
    if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_new'

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

    ! attribues
    unitStr = 'units'
    unitVal = 'm/s'
    ier = mnt_ncfieldwrite_setAttStr(handle, trim(unitStr), len_trim(unitStr), trim(unitVal), len_trim(unitVal))
    if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_setAttStr'


    ndimsStr = 'ndims'
    ndimsVal = 3
    ier = mnt_ncfieldwrite_setAttInt(handle, trim(ndimsStr), len_trim(ndimsStr), ndimsVal)
    if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_setAttInt'


    fooStr = 'foo'
    fooVal = 3.4_c_double
    ier = mnt_ncfieldwrite_setAttDbl(handle, trim(fooStr), len_trim(fooStr), fooVal)
    if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_setAttDbl'

    ! fill in
    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                u1(i, j, k) = i-1 + (j-1)*nx + (k-1)*nx*ny
            enddo
        enddo
    enddo

    ! write all the data
    ier = mnt_ncfieldwrite_data(handle, u1)
    if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_data ier = ', ier

    ! write a slice of the data
    start = [1, 1, 1]
    nsizes = [nx-1, ny-1, nz-1]
    u1 = -1
    ier = mnt_ncfieldwrite_dataSlice(handle, start, nsizes, u1)
    if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_data ier = ', ier



    ier = mnt_ncfieldwrite_del(handle)
    if (ier /= 0) print*, 'ERROR after mnt_ncfieldwrite_del'


end program