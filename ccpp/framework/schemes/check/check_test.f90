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
!! @brief A checking physics modules.
!!
!
module check_test

    implicit none

    private
    public :: test

    contains

    subroutine test(gravity, u, v, surf_t)
        implicit none
        real, intent(inout) :: gravity
        real, intent(inout) :: surf_t(:)
        real, intent(inout) :: u(:,:,:)
        real, intent(inout) :: v(:,:,:)

        print *, 'In physics test_run'
        print *, 'gravity: ', gravity
        print *, 'surf_t:  ', surf_t
        print *, 'updating u to be 10m/s'
        u = 10.0
        print *, 'updating v to be -10m/s'
        v = -10.0

    end subroutine test

end module check_test
