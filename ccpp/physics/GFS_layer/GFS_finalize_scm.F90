module GFS_finalize_scm

  implicit none

  private

!----------------
! Public entities
!----------------
  public  GFS_finalize_scm_init, GFS_finalize_scm_run, GFS_finalize_scm_finalize

  CONTAINS
!*******************************************************************************************

!--------------
! GFS initialze
!--------------

  subroutine GFS_finalize_scm_init()
  end subroutine GFS_finalize_scm_init

  subroutine GFS_finalize_scm_finalize()
  end subroutine GFS_finalize_scm_finalize

!> \section arg_table_GFS_finalize_scm_run Argument Table
!! | local_name           | standard_name                                               | long_name                                                               | units         | rank | type                          |    kind   | intent | optional |
!! |----------------------|-------------------------------------------------------------|-------------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | errmsg               | error_message                                               | error message for error handling in CCPP                                | none          |    0 | character                     | len=*     | out    | F        |
!! | errflg               | error_flag                                                  | error flag for error handling in CCPP                                   | flag          |    0 | integer                       |           | out    | F        |
!!
  subroutine GFS_finalize_scm_run (errmsg, errflg)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine GFS_finalize_scm_run

end module GFS_finalize_scm
