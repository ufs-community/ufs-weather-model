!>
!! @brief The function pointer module.
!!
!! @details The routines for calling the specified functions.
!!          This module contains no subroutines or functions it
!!          only provies an interface to the C counterparts.
!
module ccpp_dl

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_int32_t, c_char, c_ptr

    implicit none

    private
    public :: ccpp_dl_open,                                            &
              ccpp_dl_close,                                           &
              ccpp_dl_call

    interface
        integer(c_int32_t)                                             &
        function ccpp_dl_open                                          &
                   (name, library, version, shdl, lhdl)                &
                   bind(c, name='ccpp_dl_open')
            import :: c_char, c_int32_t, c_ptr
            character(kind=c_char), dimension(*)  :: name
            character(kind=c_char), dimension(*)  :: library
            character(kind=c_char), dimension(*)  :: version
            type(c_ptr)                           :: shdl
            type(c_ptr)                           :: lhdl
        end function ccpp_dl_open

        integer(c_int32_t)                                             &
        function ccpp_dl_close                                         &
                   (lhdl)                                              &
                   bind(c, name='ccpp_dl_close')
            import :: c_int32_t, c_ptr
            type(c_ptr)                           :: lhdl
        end function ccpp_dl_close

        integer(c_int32_t)                                             &
        function ccpp_dl_call                                          &
                   (shdl, cdata)                                       &
                   bind(c, name='ccpp_dl_call')
            import :: c_int32_t, c_ptr
            type(c_ptr)                           :: shdl
            type(c_ptr)                           :: cdata
        end function ccpp_dl_call
    end interface

    contains

end module ccpp_dl
