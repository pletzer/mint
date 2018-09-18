program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod

    implicit none
    character(len=32) :: f_string
    character(len=32) :: c_string

    f_string = '0123456789'
    print *,'f_string = ', f_string
    call mnt_f2c_string(f_string, c_string)
    print *,'c_string = ', c_string

end program
