program test
    use mnt_logger_capi_mod
    use, intrinsic :: iso_c_binding, only: c_size_t
    implicit none
    character(len=32) :: filename = 'test_loggerF.log'
    integer(c_size_t) :: n
    ! empty
    call mnt_printLogMessages
    n = len_trim(filename)
    call mnt_writeLogMessages(filename, n)
end program test