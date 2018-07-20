module example_ccpp_host_cap

  use ccpp_types,         only: ccpp_t
  use ccpp,               only: ccpp_init, ccpp_finalize
  use ccpp_fcall,         only: ccpp_run
  use ccpp_fields,        only: ccpp_field_add
  use iso_c_binding,      only: c_loc

! Include auto-generated list of modules for ccpp
#include "ccpp_modules.inc"

  implicit none

!  CCPP data structure
  type(ccpp_t), save, target :: cdata

  public :: physics_init, physics_run, physics_finalize

contains

  subroutine physics_init(ccpp_suite_name)
    character(len=*), intent(in) :: ccpp_suite_name
    integer :: ierr
    ierr = 0

    call ccpp_init(ccpp_suite_name, cdata, ierr)
    if (ierr/=0) then
      write(*,'(a)') "An error occurred in ccpp_init"
      stop
    end if

! Include auto-generated list of calls to ccpp_field_add
#include "ccpp_fields.inc"

  end subroutine physics_init

  subroutine physics_run(step)
    ! the step (currently called IPD step) to run as
    ! defined in the runtime suite definition file
    integer, intent(in) :: step
    integer :: ierr
    ierr = 0

    call ccpp_run(cdata%suite%ipds(step), cdata, ierr)
    ! future versions: call ccpp_run(cdata, step=step, ierr=ierr)
    if (ierr/=0) then
      ! errmsg is known because of #include ccpp_modules.inc
      write(*,'(a,i0,a)') "An error occurred in physics step ", step, &
                          "; error message: '" // trim(errmsg) // "'"
      stop
    end if

  end subroutine physics_run

  subroutine physics_finalize()
    integer :: ierr
    ierr = 0

    call ccpp_finalize(cdata, ierr)
    if (ierr/=0) then
      write(*,'(a)') "An error occurred in ccpp_finalize"
      stop
    end if

  end subroutine physics_finalize

end module example_ccpp_host_cap
