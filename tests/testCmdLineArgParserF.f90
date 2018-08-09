program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod

    use, intrinsic :: iso_c_binding, only: c_size_t, c_int, c_double

    implicit none
    integer(c_size_t) :: prsr
    integer(c_int)    :: ier

    ier = mnt_cmdlineargparser_new(prsr)

    ier = mnt_cmdlineargparser_del(prsr)

end program 
