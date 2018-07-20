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
!! @brief The CCPP library main entry and exit points.
!!
!
module ccpp

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_ptr
    use            :: ccpp_types,                                      &
                      only: ccpp_t, ccpp_suite_t
    use            :: ccpp_suite,                                      &
                      only: ccpp_suite_init, ccpp_suite_finalize
    use            :: ccpp_fields,                                     &
                      only: ccpp_fields_init, ccpp_fields_finalize
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug

    implicit none

    private
    public :: ccpp_init,                                               &
              ccpp_finalize,                                           &
              ccpp_initialized

    contains

    !>
    !! CCPP initialization subroutine.
    !!
    !! @param[in]     filename The file name of the XML scheme file to load.
    !! @param[in,out] cdata    The ccpp_t type data.
    !! @param[  out]  ierr     Integer error flag.
    !
    subroutine ccpp_init(filename, cdata, ierr, suite)
        character(len=*),   intent(in)           :: filename
        type(ccpp_t),       intent(inout)        :: cdata
        integer,            intent(  out)        :: ierr
        type(ccpp_suite_t), intent(in), optional :: suite

        ierr = 0

        call ccpp_debug('Called ccpp_init')

        if (present(suite)) then
            ! Makes a copy of the suite to avoid multiple
            ! reads/parses of the suite definiton file
            cdata%suite = suite
            cdata%suite%iscopy = .True.
        else
            ! Initialize the suite
            call ccpp_suite_init(filename, cdata%suite, ierr)
            if (ierr /= 0) then
                call ccpp_error('In initializing the CCPP suite')
                return
            end if
        end if

        ! Initialize the fields
        call ccpp_fields_init(cdata, ierr)
        if (ierr /= 0) then
            call ccpp_error('In initializing the CCPP fields')
            return
        end if

        ! Set flag indicating initialization state of cdata
        cdata%initialized = .true.

    end subroutine ccpp_init

    !>
    !! CCPP finalization subroutine.
    !!
    !! @param[in,out] cdata    The ccpp_t type data.
    !! @param[  out]  ierr     Integer error flag.
    !
    subroutine ccpp_finalize(cdata, ierr)
        type(ccpp_t),           intent(inout) :: cdata
        integer,                intent(  out) :: ierr

        ierr = 0

        call ccpp_debug('Called ccpp_finalize')

        ! Finalize the suite
        call ccpp_suite_finalize(cdata%suite, ierr)
        if (ierr /= 0) then
                call ccpp_error('In finalizing the CCPP suite')
                return
        end if

        ! Finalize the fields
        call ccpp_fields_finalize(cdata, ierr)
        if (ierr /= 0) then
                call ccpp_error('In finalizing the CCPP fields')
                return
        end if

        ! Set flag indicating initialization state of cdata
        cdata%initialized = .false.

    end subroutine ccpp_finalize

    !>
    !! CCPP test initialization routine
    !!
    !! @param[in]     cdata        The ccpp_t type data
    !! @return        initialized  .true. or .false.
    !
    function ccpp_initialized(cdata) result(initialized)
        type(ccpp_t), intent(in) :: cdata
        logical                  :: initialized

        call ccpp_debug('Called ccpp_initialized')

        initialized = cdata%initialized

    end function ccpp_initialized

end module ccpp
