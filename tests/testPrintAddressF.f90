program test
    ! defines the C API
    use mnt_celllocator_capi_mod
    use, intrinsic :: iso_c_binding
    implicit none
    real*8, target :: point(3)
    point(1) = 1; point(2) = 2; point(3) = 3;
    call mnt_celllocator_printaddress(c_loc(point))

end program 
