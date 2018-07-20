!>
!! @brief The function pointer module.
!!
!! @details The routines for calling the specified functions.
!!          This module contains no subroutines or functions it
!!          only provies an interface to the C counterparts.
!
module ccpp_memory

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_int32_t, c_char, c_null_char

    implicit none

    private
    public :: ccpp_memory_usage

    interface
        integer(c_int32_t)                                             &
        function ccpp_memory_usage_c                                   &
                   (mpicomm, str, lstr)                                &
                   bind(c, name='ccpp_memory_usage_c')
            import :: c_char, c_int32_t
            integer(c_int32_t), value, intent(in) :: mpicomm
            character(kind=c_char), dimension(*)  :: str
            integer(c_int32_t), value, intent(in) :: lstr
        end function ccpp_memory_usage_c
    end interface

    contains

    function ccpp_memory_usage(mpicomm, memory_usage) result(ierr)

        implicit none

        ! Interface variables
        integer, intent(in)                          :: mpicomm
        character(len=*), intent(out)                :: memory_usage
        ! Function return value
        integer                                      :: ierr
        ! Local variables
        character(len=len(memory_usage),kind=c_char) :: memory_usage_c
        integer                                      :: i

        ierr = ccpp_memory_usage_c(mpicomm, memory_usage_c, len(memory_usage_c))
        if (ierr /= 0) then
            write(memory_usage,fmt='(a)') "An error occurred in the call to ccpp_memory_usage_c in ccpp_memory_usage"
            return
        end if

        memory_usage = memory_usage_c(1:index(memory_usage_c, c_null_char)-1)

    end function ccpp_memory_usage

end module ccpp_memory
