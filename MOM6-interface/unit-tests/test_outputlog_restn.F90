!> @file test_outputlog_restn.F90
!> @brief Test per-part completion for track_restn
!!
!! Calls track_restn directly,  against real restart-fixture files to verify
!! both per-part allDone(n) and all(allDone). Also tests filenames constructed
!! by track_restn
!!
!> @date 09-15-2026

!> Main program for testing outputlog_restn tracking
program test_outputlog_restn

  use ESMF
  use mpi_f08,                only : MPI_Init, MPI_Comm, MPI_Comm_rank, MPI_COMM_WORLD, MPI_Barrier
  use test_utils
  use mom_outputlog_methods,  only : track_restn
  use nc_fixture_mod,         only : restart_part_fname, make_restart_fixture
  use test_helpers,           only : base_yy, base_mm, base_dd

  implicit none

  integer, parameter :: maxtests = 20

  type(MPI_Comm) :: comm
  integer        :: rank, ierr, rootpe
  logical        :: isroot

  character(len=128) :: testname
  character(len=256) :: assertmsg
  character(len=20)  :: subname = 'test_track_restn'
  character(len=256) :: restartdir = './'

  type(testsummary) :: restntests

  logical :: assertrc
  integer :: n, nt

  comm = MPI_COMM_WORLD
  rootpe = 0

  call restntests%init(maxtests)

  call MPI_Init(ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  isroot = (rank == rootpe)
  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=ierr)
  call esmf_err(ierr, subname, "ESMF_Initialize")

  if (isroot) call execute_command_line('rm -f '//trim(restartdir)//'*.MOM.res*.nc', wait=.true.)
  call MPI_Barrier(comm, ierr)
  if (ierr /= 0) then
     write(0,'(A)') "FATAL ("//trim(subname)//"): MPI_Barrier (post-cleanup) failed"
     stop 99
  endif

  nt = 0
  ! ===========================================================================
  ! Test cases
  ! ===========================================================================

  nt = nt + 1
  write(testname,'(A,I2.2,A)')'test ',nt,' test track_restn: 3 parts, 2nd incomplete '
  call check_partial(trim(testname),hour=0)

  ! ------------------
  nt = nt + 1
  write(testname,'(A,I2.2,A)')'test ',nt,' test track_restn: 3 parts, all complete '
  call check_all_complete(trim(testname),hour=6)

  ! ------------------
  nt = nt + 1
  write(testname,'(A,I2.2,A)')'test ',nt,' test track_restn: 1 part, complete '
  call check_single_file(trim(testname), hour=12, complete=.true.)

  ! ------------------
  nt = nt + 1
  write(testname,'(A,I2.2,A)')'test ',nt,' test track_restn: 1 part, incomplete '
  call check_single_file(trim(testname),hour=18, complete=.false.)

  ! ------------------
  ! Test results
  ! ------------------
  if (isroot) then
  if (restntests%nfail > 0) then
     print '(A)', 'FAIL: At least one test failed '
     do n = 1,restntests%count
        if (.not. restntests%teststatus(n)) print '(A)', trim(restntests%testmessage(n)%str)
     enddo
  else
     do n = 1,restntests%count
        print '(A)', trim(restntests%testmessage(n)%str)
     enddo
  endif
  print '(3(A,I0))','Total tests = ',restntests%count,' Passing = ',restntests%npass, &
       ' Failing = ',restntests%nfail
  endif

  call MPI_Barrier(comm, ierr)
  if (isroot) call execute_command_line('rm -f '//trim(restartdir)//'*.MOM.res*.nc', wait=.true.)

  call ESMF_Finalize(rc=ierr)
  call esmf_err(ierr, subname, "ESMF_Finalize")

  if (restntests%nfail > 0) then
     if (isroot) print '(A)','Test failures! '
     stop 1
  endif

contains
  !> Construct filename independently of track_restn
  !!
  !! @param[in]     hour          filename hour
  !! @param[in]     part_index    filename part
  !! @return        fname         constructed filename
  function expected_fname(hour, part_index) result(fname)
    integer, intent(in) :: hour
    integer, intent(in) :: part_index

    character(len=256) :: fname
    type(ESMF_Time)    :: nextTime
    integer :: yr, mon, day, hr, minute, sec, rc

    character(len=15)  :: base_timestr
    character(len=256) :: base

    call ESMF_TimeSet(nextTime, yy=base_yy, mm=base_mm, dd=base_dd, h=hour, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeSet(nextTime)")
    call ESMF_TimeGet(nextTime, yy=yr, mm=mon, dd=day, h=hr, m=minute, s=sec, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeGet(nextTime)")

    write(base_timestr,'(I4.4,2(I2.2),A,3(I2.2))') yr, mon, day, ".", hr, minute, sec
    base = trim(restartdir)//trim(base_timestr)//'.MOM.res'
    fname = restart_part_fname(base, part_index)
  end function expected_fname
  !> Check partial completion case
  !!
  !! @param[in]   test    test identifier string
  !! @param[in]   hour    filename hour
  subroutine check_partial(test,hour)
    character(len=*), intent(in) :: test
    integer,          intent(in) :: hour

    integer, parameter :: num_rest_files = 3
    logical, parameter :: part_complete(num_rest_files) = [.true., .false., .true.]

    integer :: n, rc
    type(ESMF_Time) :: nextTime
    logical, allocatable :: allDone(:)
    character(len=256), allocatable :: fnames(:)
    character(len=256) :: fname

    call ESMF_TimeSet(nextTime, yy=base_yy, mm=base_mm, dd=base_dd, h=hour, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeSet(nextTime)")

    ! create three restart files, with 2nd file incomplete
    do n = 1, num_rest_files
       fname = expected_fname(hour, n-1)
       if (isroot) call make_restart_fixture(fname, complete=part_complete(n))
    enddo
    call MPI_Barrier(comm, ierr)

    call track_restn(nextTime, num_rest_files, comm, isroot, rootpe, restartdir, allDone, fnames, rc)
    call esmf_err(rc, subname, "track_restn (partial completion)")

    do n = 1, num_rest_files
       call assert_equal(trim(fnames(n))==trim(expected_fname(hour,n-1)), .true., &
            test//', check fname(n)', assertrc, assertmsg)
       call addresult(restntests, assertrc, trim(assertmsg), '')
    enddo
    do n = 1, num_rest_files
       call assert_equal(allDone(n), part_complete(n), test//', check alldone(n)', assertrc, assertmsg)
       call addresult(restntests, assertrc, trim(assertmsg), '')
    enddo

    call assert_equal(all(allDone), .false., test//', check all(allDone)', assertrc, assertmsg)
    call addresult(restntests, assertrc, trim(assertmsg), '')
  end subroutine check_partial
  !> Check full completion case
  !!
  !! @param[in]   test    test identifier string
  !! @param[in]   hour    filename hour
  subroutine check_all_complete(test,hour)
    character(len=*), intent(in) :: test
    integer,          intent(in) :: hour

    integer, parameter :: num_rest_files = 3

    type(ESMF_Time) :: nextTime
    logical, allocatable :: allDone(:)
    character(len=256), allocatable :: fnames(:)
    integer :: n, rc
    character(len=256) :: fname

    call ESMF_TimeSet(nextTime, yy=base_yy, mm=base_mm, dd=base_dd, h=hour, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeSet(nextTime)")

    ! create three restart files, all complete
    do n = 1, num_rest_files
       fname = expected_fname(hour, n-1)
       if (isroot) call make_restart_fixture(fname, complete=.true.)
    enddo
    call MPI_Barrier(comm, ierr)

    call track_restn(nextTime, num_rest_files, comm, isroot, rootpe, restartdir, allDone, fnames, rc)
    call esmf_err(rc, subname, "track_restn (all complete)")

    do n = 1, num_rest_files
       call assert_equal(allDone(n), .true., test//', check alldone(n)', assertrc, assertmsg)
       call addresult(restntests, assertrc, trim(assertmsg), '')
    enddo

    call assert_equal(all(allDone), .true., test//', check all(allDone)', assertrc, assertmsg)
    call addresult(restntests, assertrc, trim(assertmsg), '')

  end subroutine check_all_complete
  !> Check single file, either complete or not
  !!
  !! @param[in]   test     test identifier string
  !! @param[in]   hour     filename hour
  !! @param[in]   complete logical to stage either complete or incomplete file
  subroutine check_single_file(test,hour,complete)
    character(len=*), intent(in) :: test
    integer,          intent(in) :: hour
    logical,          intent(in) :: complete

    integer, parameter :: num_rest_files = 1

    type(ESMF_Time) :: nextTime
    logical, allocatable :: allDone(:)
    character(len=256), allocatable :: fnames(:)
    integer :: rc
    character(len=256) :: fname
    character(len=16) :: tag

    call ESMF_TimeSet(nextTime, yy=base_yy, mm=base_mm, dd=base_dd, h=hour, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeSet(nextTime)")

    fname = expected_fname(hour, 0)
    ! create single restart files, either complete or not
    if (isroot) call make_restart_fixture(fname, complete=complete)
    call MPI_Barrier(comm, ierr)

    call track_restn(nextTime, num_rest_files, comm, isroot, rootpe, restartdir, allDone, fnames, rc)
    call esmf_err(rc, subname, "track_restn (single file)")

    if (complete) then
       tag = ' complete'
    else
       tag = ' incomplete'
    endif

    call assert_equal(allDone(1), complete, test//trim(tag)//', check alldone(1)', assertrc, assertmsg)
    call addresult(restntests, assertrc, trim(assertmsg), '')

    call assert_equal(all(allDone), complete, test//trim(tag)//', check all(allDone)', assertrc, assertmsg)
    call addresult(restntests, assertrc, trim(assertmsg), '')
  end subroutine check_single_file

end program test_outputlog_restn
