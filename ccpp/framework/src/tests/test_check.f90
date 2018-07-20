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
!! @brief A Test Atmospheric Driver Program.
!!
program test_check

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_loc, c_f_pointer
    use            :: ccpp_api,                                        &
                      only: CCPP_STR_LEN,                              &
                            ccpp_t,                                    &
                            ccpp_init,                                 &
                            ccpp_finalize,                             &
                            ccpp_physics_run,                          &
                            ccpp_field_add

    implicit none

    type(ccpp_t), target                         :: cdata
    character(len=CCPP_STR_LEN)                  :: filename
    integer                                      :: len
    integer                                      :: ierr
    integer                                      :: asize
    real, target                                 :: gravity
    real, target, allocatable, dimension(:)      :: surf_t
    real, target, allocatable, dimension(:,:,:)  :: u
    real, target, allocatable, dimension(:,:,:)  :: v

    ierr = 0

    call get_command_argument(1, filename, len, ierr)
    if (ierr /= 0) then
            call exit(1)
    end if

    ! Allocate the data
    asize = 5
    allocate(surf_t(asize), stat=ierr)
    if (ierr /= 0) then
            print *, 'Unable to allocate surface temperature array'
            call exit(1)
    end if

    allocate(u(asize,asize,asize), stat=ierr)
    if (ierr /= 0) then
            print *, 'Unable to allocate U array'
            call exit(1)
    end if

    allocate(v(asize,asize,asize), stat=ierr)
    if (ierr /= 0) then
            print *, 'Unable to allocate U array'
            call exit(1)
    end if

    ! Generate data to pass into a physics driver
    gravity = 9.80665
    surf_t = [290.0, 291.0, 292.0, 293.0, 294.0]
    u = 0.0
    v = 10.0

    ! Initalize the CCPP (with the filename of the suite to load).
    call ccpp_init(filename, cdata, ierr)
    if (ierr /= 0) then
        call exit(1)
    end if

    ! Add all the fields we want to expose to the physics driver.
    call ccpp_field_add(cdata, 'gravity', gravity, ierr, 'm s-2')
    if (ierr /= 0) then
            call exit(1)
    end if

    call ccpp_field_add(cdata, 'surface_temperature', surf_t, ierr, 'K')

    call ccpp_field_add(cdata, 'eastward_wind', u, ierr, 'm s-1')

    call ccpp_field_add(cdata, 'northward_wind', v, ierr, 'm s-1')

    call ccpp_physics_run(cdata, scheme_name="test", ierr=ierr)

    print *, 'In test dummy main'
    print *, 'gravity: ', gravity
    print *, 'surf_t:  ', surf_t(1:2)
    print *, 'u: ', u(1,1,1)
    print *, 'v: ', v(1,1,1)

    call ccpp_finalize(cdata, ierr)

    if (allocated(surf_t)) then
        deallocate(surf_t)
    end if

    if (allocated(u)) then
        deallocate(u)
    end if

    if (allocated(v)) then
        deallocate(v)
    end if

end program test_check
