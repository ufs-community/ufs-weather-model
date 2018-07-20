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
!! @brief A NO-OP physics modules.
!!
!
module check_noop

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr
    use            :: ccpp_types,                                      &
                      only: ccpp_t
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    implicit none

    private
    public :: noop_cap

    contains

    subroutine noop_cap(ptr) bind(c)
        implicit none
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer      :: cdata

        call c_f_pointer(ptr, cdata)

        print *, 'In noop_cap'
        print *, cdata%suite%groups(1)%subcycles(1)%schemes(1)%name

    end subroutine noop_cap

end module check_noop
