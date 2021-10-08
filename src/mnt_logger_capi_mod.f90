module mnt_logger_capi_mod

  ! C function prototypes
  interface

    subroutine mnt_printLogMessages() bind(C, name='mnt_printLogMessages')
      ! Print log messages
    end subroutine mnt_printLogMessages

    subroutine mnt_writeLogMessages(filename, filename_len) bind(C, name='mnt_writeLogMessages')
      ! Write log messages to file
      ! @param filename file name
      ! @param filename_len number of characters in filename, excluding '\0'
      use, intrinsic :: iso_c_binding, only: c_char, c_size_t
      implicit none
      character(kind=c_char), intent(in)       :: filename(*)
      integer(c_size_t), value                 :: filename_len
    end subroutine mnt_writeLogMessages

  end interface

end module mnt_logger_capi_mod

