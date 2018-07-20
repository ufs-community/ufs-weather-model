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
!! @brief A test program to test the CCPP.
!!
!! @details This will test the
!!          - initialization and finalization subroutines of
!!            -- CCPP
!!            -- Suite
!!            -- Fields
!!          It can be used multipile times to test the parsing
!!          of various suite XML files.
!
program test_init_finalize

    use            :: ccpp_types,                          &
                      only: CCPP_STR_LEN, ccpp_t
    use            :: ccpp,                                &
                      only: ccpp_init, ccpp_finalize

    implicit none

    integer                :: ierr
    integer                :: len
    character(len=CCPP_STR_LEN) :: filename
    type(ccpp_t), target   :: cdata


    ierr = 0

    call get_command_argument(1, filename, len, ierr)
    if (ierr /= 0) then
        print *, 'Error: no suite XML file specified.'
        call exit(ierr)
    end if

    call ccpp_init(filename, cdata, ierr)
    if (ierr /= 0) then
        call exit(ierr)
    end if

    call ccpp_finalize(cdata, ierr)
    if (ierr /= 0) then
        call exit(ierr)
    end if

end program test_init_finalize
