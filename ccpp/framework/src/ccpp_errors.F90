!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Error/Warning reporting module.
!!
!! @details Subroutines for reporting warnings.
!
module ccpp_errors

    use, intrinsic :: iso_fortran_env,                                 &
                      only: error_unit, output_unit

    implicit none

    private
    public :: ccpp_error,                                              &
              ccpp_warn,                                               &
              ccpp_info,                                               &
              ccpp_debug,                                              &
              ccpp_if_error,                                           &
              ccpp_if_warn

    contains

    !>
    !! Fatal error reporting.
    !!
    !! Write an error message to error_unit/stderr.
    !!
    !! @param[in] message   The error message to write.
    !
    subroutine ccpp_error(message)
        character(len=*),        intent(in) :: message

        write(error_unit, *) 'ERROR: ', trim(message)
    end subroutine ccpp_error

    !>
    !! Non-fatal warning reporting.
    !!
    !! Write an warning message to error_unit/stderr.
    !!
    !! @param[in] message   The warning message to write.
    !
    subroutine ccpp_warn(message)
        character(len=*),        intent(in) :: message

        write(error_unit, *) 'WARN: ', trim(message)
    end subroutine ccpp_warn

    !>
    !! Reporting on info level
    !!
    !! Write an info message to output_unit/stdout.
    !!
    !! @param[in] message   The info message to write.
    !
    subroutine ccpp_info(message)
        character(len=*),        intent(in) :: message

        write(output_unit, *) 'INFO: ', trim(message)
    end subroutine ccpp_info

    !>
    !! Reporting on debug level
    !!
    !! Write a debug message to output_unit/stdout.
    !!
    !! @param[in] message   The debug message to write.
    !
    subroutine ccpp_debug(message)
        character(len=*),        intent(in) :: message

#ifdef DEBUG
        write(output_unit, *) 'DEBUG: ', trim(message)
#endif
    end subroutine ccpp_debug

    !>
    !! Fatal error checking and reporting.
    !!
    !! Check to see if ierr is non-zero. If it is
    !! write an error message to error_unit/stderr.
    !!
    !! @param[in] ierr      The exit code.
    !! @param[in] message   The error message to write.
    !
    subroutine ccpp_if_error(ierr, message)
        integer,                 intent(in) :: ierr
        character(len=*),        intent(in) :: message

        if (ierr /= 0) then
            write(error_unit, *) 'ERROR: ', trim(message)
        end if

    end subroutine ccpp_if_error

    !>
    !! Non-fatal warning checking and reporting.
    !!
    !! Check to see if ierr is non-zero. If it is
    !! write an warning message to error_unit/stderr.
    !!
    !! @param[in] ierr      The exit code.
    !! @param[in] message   The warning message to write.
    !
    subroutine ccpp_if_warn(ierr, message)
        integer,                 intent(in) :: ierr
        character(len=*),        intent(in) :: message

        if (ierr /= 0) then
            write(error_unit, *) 'WARN: ', trim(message)
        end if

    end subroutine ccpp_if_warn

end module ccpp_errors
