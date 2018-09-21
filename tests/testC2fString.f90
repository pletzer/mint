program test
    ! defines the C API
    use mnt_cmdlineargparser_capi_mod

    implicit none
    character(len=20) :: f_string, f_string2
    character(len=1)  :: c_string(32)

    f_string = '0123456789'
    print *,'f_string = "', f_string, '"'
    call mnt_f2c_string(f_string, c_string)
    print *,'c_string = "', c_string, '"'
    call mnt_c2f_string(c_string, f_string2)
    print *,'f_string2 = "', f_string2, '"'

    if (f_string2 /= f_string) then
        print *,'ERROR: f_string and f_string2 do not match'
        stop 1
    endif

end program
