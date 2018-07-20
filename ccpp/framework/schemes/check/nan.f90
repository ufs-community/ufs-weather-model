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
!! @brief A physics module to check for NaNs.
!!
!
module check_nans

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr
    use            :: ccpp_types,                                      &
                      only: ccpp_t
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    implicit none

    private
    public :: nans_cap

    contains

    subroutine nans_cap(ptr) bind(c)
        implicit none
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t),     pointer  :: cdata
        real,             pointer  :: v(:,:,:)
        integer                    :: i
        integer                    :: ierr

        call c_f_pointer(ptr, cdata)

        call ccpp_field_get(cdata, 'northward_wind', v, ierr)

        call test_run(gravity, u, v, surf_t)

    end subroutine nans_cap

    subroutine run(gravity, u, v, surf_t)
        implicit none
        real, pointer, intent(inout) :: gravity
        real, pointer, intent(inout) :: surf_t(:)
        real, pointer, intent(inout) :: u(:,:,:)
        real, pointer, intent(inout) :: v(:,:,:)

        print *, 'In physics check nans run'

    end subroutine test_run

end module check_test
