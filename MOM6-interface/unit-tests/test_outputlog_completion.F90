!> @file test_outputlog_completion.F90
!> @brief Tests file completion contract in isolation
!!
!! Probes get_file_state and file_is_complete directly against on-disk netCDF-4 fixture
!! files (see nc_fixture_mod.F90), independent of the alarm/clock machinery in outputlog_freqn.
!!
!! @date 09-01-2026

!> Main program for testing outputlog file completion logic
program test_outputlog_completion

  use ESMF
  use mpi_f08,               only : MPI_Init, MPI_Comm, MPI_Comm_rank, MPI_COMM_WORLD, MPI_Barrier
  use mom_outputlog_methods, only : get_file_state, file_is_complete, set_restfname
  use nc_fixture_mod,        only : make_datm_incomplete, make_datm_complete
  use nc_fixture_mod,        only : make_atm_incomplete,  make_atm_complete
  use nc_fixture_mod,        only : make_restart_fixture
  use test_helpers,          only : base_yy, base_mm, base_dd

  implicit none

  type(MPI_Comm) :: comm
  integer        :: rank, ierr, rootpe
  logical        :: isroot
  integer        :: total_errors
  logical        :: verbose = .false.

  comm   = MPI_COMM_WORLD
  rootpe = 0
  total_errors = 0

  call MPI_Init(ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  isroot = (rank == rootpe)
  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=ierr)
  if (ierr /= ESMF_SUCCESS) then
    write(0,'(A)') "FATAL (test_outputlog_completion): ESMF_Initialize failed"
    stop 99
  end if
  ! cleanup state files (if run outside of CI)
  if (isroot) call execute_command_line('rm -f test_*.nc *.MOM.res*.nc', wait=.true.)
  call MPI_Barrier(comm, ierr)

  call check_datm_incomplete()
  call check_datm_complete()
  call check_atm_incomplete()
  call check_atm_complete()
  call check_file_does_not_exist()
  call check_restart_single_file_complete()
  call check_restart_single_file_incomplete()
  call check_restart_multiple_files_all_complete()
  call check_restart_multiple_files_partial()

  if (isroot) then
    print *, "========================================================"
    if (total_errors == 0) then
      print *, "SUCCESS: all completion-contract cases passed"
    else
      print *, "FAILURE: ", total_errors, " assertions failed"
    end if
    print *, "========================================================"
  end if
  call MPI_Barrier(comm, ierr)
  if (isroot) call execute_command_line('rm -f test_*.nc *.MOM.res*.nc', wait=.true.)

  call ESMF_Finalize(rc=ierr)
  if (ierr /= ESMF_SUCCESS .and. isroot) then
    write(0,'(A)') "WARNING (test_outputlog_completion): ESMF_Finalize returned an error"
  end if
  if (total_errors == 0) then
    stop 0
  else
    stop 1
  end if

contains
  !> Checks logic for an incomplete DATM file
  !!
  subroutine check_datm_incomplete()
    character(len=*), parameter :: fname = "test_datm_incomplete.nc"
    integer :: rc, nlen, fsize
    logical :: complete

    if (isroot) call make_datm_incomplete(fname)
    call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, rc=rc)

    call assert_equal(0, rc, "DATM incomplete: get_file_state rc")
    call assert_equal(0, nlen, "DATM incomplete: nlen should be 0")

    complete = file_is_complete(comm, isroot, rootpe, fname, .false., 0, rc)
    if (isroot .and. verbose) print '(A,i6,A,L2)',fname//',  nlen = ',nlen,', fsize value ignored, complete = ',complete

    call assert_equal(0, rc, "DATM incomplete: file_is_complete rc")
    call assert_false(complete, "DATM incomplete: should NOT be complete (nlen=0)")
  end subroutine check_datm_incomplete
  !> Checks logic for a complete DATM file
  !!
  subroutine check_datm_complete()
    character(len=*), parameter :: fname = "test_datm_complete.nc"
    integer :: rc, nlen, fsize
    logical :: complete

    if (isroot) call make_datm_complete(fname)

    call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, fsize=fsize, rc=rc)

    call assert_equal(0, rc, "DATM complete: get_file_state rc")
    call assert_equal(1, nlen, "DATM complete: nlen should be 1")

    complete = file_is_complete(comm, isroot, rootpe, fname, .false., 0, rc)
    if (isroot .and. verbose) print '(A,i6,A,L2)',fname//',  nlen = ',nlen,', fsize value ignored, complete = ',complete

    call assert_equal(0, rc, "DATM complete: file_is_complete rc")
    call assert_true(complete, "DATM complete: SHOULD be complete (nlen=1, use_filesize=F)")
  end subroutine check_datm_complete
  !> Checks logic for an incomplete ATM file
  !!
  subroutine check_atm_incomplete()
    character(len=*), parameter :: fname = "test_atm_incomplete.nc"
    integer :: rc, nlen, fsize, createsize
    logical :: complete

    createsize = 0
    if (isroot) call make_atm_incomplete(fname, createsize)

    call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, fsize=fsize, rc=rc)

    call assert_equal(0, rc, "ATM incomplete: get_file_state rc")
    call assert_equal(1, nlen, "ATM incomplete: nlen should be 1 (record written)")
    call assert_equal(createsize, fsize, "ATM incomplete: fsize should equal createsize (no bulk data yet)")

    complete = file_is_complete(comm, isroot, rootpe, fname, .true., createsize, rc)
    if (isroot .and. verbose) print '(2(A,i6),A,L2)',fname//',  nlen = ',nlen,', fsize = ',fsize,', complete = ',complete

    call assert_equal(0, rc, "ATM incomplete: file_is_complete rc")
    call assert_false(complete, "ATM incomplete: should NOT be complete (nlen>0 but size==createsize)")
  end subroutine check_atm_incomplete
  !> Checks logic for a complete ATM file
  !!
  subroutine check_atm_complete()
    character(len=*), parameter :: fname = "test_atm_complete.nc"
    integer :: rc, nlen, fsize, createsize
    logical :: complete

    createsize = 0
    if (isroot) call make_atm_complete(fname, createsize)

    call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, fsize=fsize, rc=rc)

    call assert_equal(0, rc, "ATM complete: get_file_state rc")
    call assert_equal(1, nlen, "ATM complete: nlen should be 1")
    call assert_true(fsize > createsize, "ATM complete: fsize should exceed createsize (bulk data written)")

    complete = file_is_complete(comm, isroot, rootpe, fname, .true., createsize, rc)
    if (isroot .and. verbose) print '(2(A,i6),A,L2)',fname//',  nlen = ',nlen,', fsize = ',fsize,', complete = ',complete

    call assert_equal(0, rc, "ATM complete: file_is_complete rc")
    call assert_true(complete, "ATM complete: SHOULD be complete (nlen>0 and size>createsize)")
  end subroutine check_atm_complete
  !> Edge case: get_file_state/file_is_complete asked about a file that was
  !> never created. Confirm that propagates correctly and that file_is_complete
  !> treats it as incomplete either way.
  !!
  subroutine check_file_does_not_exist()
    character(len=*), parameter :: fname = "test_does_not_exist.nc"
    integer :: rc, nlen
    logical :: complete_datm, complete_atm

    call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, rc=rc)
    call assert_equal(0, rc, "Nonexistent file: get_file_state rc")
    call assert_true(nlen < 0, "Nonexistent file: nlen should be the fill-value sentinel, not a real count")

    complete_datm = file_is_complete(comm, isroot, rootpe, fname, .false., 0, rc)
    call assert_false(complete_datm, "Nonexistent file: use_filesize=F must not report complete")

    complete_atm = file_is_complete(comm, isroot, rootpe, fname, .true., 0, rc)
    call assert_false(complete_atm, "Nonexistent file: use_filesize=T must not report complete")
  end subroutine check_file_does_not_exist
  !> Restart, num_rest_files=1: the one part is complete (nlen=1).
  !!
  subroutine check_restart_single_file_complete()
    type(ESMF_Time) :: nextTime
    character(len=256) :: fname
    integer :: rc, nlen
    logical :: alldone

    call ESMF_TimeSet(nextTime, yy=base_yy, mm=base_mm, dd=base_dd, h=0, rc=rc)
    call assert_equal(0, rc, "Restart single file complete: ESMF_TimeSet rc")

    fname = set_restfname(nextTime, 1, './', rc)
    call assert_equal(0, rc, "Restart single file complete: set_restfname rc")
    if (isroot) call make_restart_fixture(fname, complete=.true.)
    call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, rc=rc)
    call assert_equal(0, rc, "Restart single file complete: get_file_state rc")
    call assert_equal(1, nlen, "Restart single file complete: nlen should be 1")

    alldone = file_is_complete(comm, isroot, rootpe, fname, .false., 0, rc)
    if (isroot .and. verbose) print '(A,i6,A,L2)',trim(fname)//',  nlen = ',nlen,', complete = ',alldone

    call assert_equal(0, rc, "Restart single file complete: file_is_complete rc")
    call assert_true(alldone, "Restart single file complete: allDone should be true")
  end subroutine check_restart_single_file_complete

  !> Restart, num_rest_files=1: the one part is NOT complete (nlen=0).
  !!
  subroutine check_restart_single_file_incomplete()
    type(ESMF_Time) :: nextTime
    character(len=256) :: fname
    integer :: rc, nlen
    logical :: alldone

    call ESMF_TimeSet(nextTime, yy=base_yy, mm=base_mm, dd=base_dd, h=6, rc=rc)
    call assert_equal(0, rc, "Restart single file incomplete: ESMF_TimeSet rc")

    fname = set_restfname(nextTime, 1, './', rc)
    call assert_equal(0, rc, "Restart single file incomplete: set_restfname rc")
    if (isroot) call make_restart_fixture(fname, complete=.false.)
    call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, rc=rc)

    call assert_equal(0, rc, "Restart single file incomplete: get_file_state rc")
    call assert_equal(0, nlen, "Restart single file incomplete: nlen should be 0")

    alldone = file_is_complete(comm, isroot, rootpe, fname, .false., 0, rc)
    if (isroot .and. verbose) print '(A,i6,A,L2)',trim(fname)//',  nlen = ',nlen,', complete = ',alldone

    call assert_equal(0, rc, "Restart single file incomplete: file_is_complete rc")
    call assert_false(alldone, "Restart single file incomplete: allDone should be false")
  end subroutine check_restart_single_file_incomplete
  !> Restart, num_rest_files=3: all three parts complete
  !!
  subroutine check_restart_multiple_files_all_complete()
    integer, parameter :: num_rest_files = 3
    type(ESMF_Time) :: nextTime
    character(len=256) :: fname
    integer :: rc, nlen, n
    logical :: alldone(num_rest_files)

    call ESMF_TimeSet(nextTime, yy=base_yy, mm=base_mm, dd=base_dd, h=12, rc=rc)
    call assert_equal(0, rc, "Restart multi-file all complete: ESMF_TimeSet rc")

    do n = 1, num_rest_files
      fname = set_restfname(nextTime, n, './', rc)
      call assert_equal(0, rc, "Restart multi-file all complete: set_restfname rc")
      if (isroot) call make_restart_fixture(fname, complete=.true.)
      call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, rc=rc)
      call assert_equal(0, rc, "Restart multi-file all complete: get_file_state rc")
      call assert_equal(1, nlen, "Restart multi-file all complete: part nlen should be 1")

      alldone(n) = file_is_complete(comm, isroot, rootpe, fname, .false., 0, rc)
      if (isroot .and. verbose) print '(A,i6,A,L2)',trim(fname)//',  nlen = ',nlen,', complete = ',alldone(n)
      call assert_equal(0, rc, "Restart multi-file all complete: file_is_complete rc")
    end do

    call assert_true(all(alldone), "Restart multi-file all complete: allDone should be true")
  end subroutine check_restart_multiple_files_all_complete
  !> Restart, num_rest_files=3: two parts complete, one not
  !!
  subroutine check_restart_multiple_files_partial()
    integer, parameter :: num_rest_files = 3
    logical, parameter :: part_complete(num_rest_files) = [.true., .false., .true.]
    type(ESMF_Time) :: nextTime
    character(len=256) :: fname
    integer :: rc, nlen, n
    logical :: alldone(num_rest_files)

    call ESMF_TimeSet(nextTime, yy=base_yy, mm=base_mm, dd=base_dd, h=18, rc=rc)
    call assert_equal(0, rc, "Restart multi-file partial: ESMF_TimeSet rc")

    do n = 1, num_rest_files
      fname = set_restfname(nextTime, n, './', rc)
      call assert_equal(0, rc, "Restart multi-file partial: set_restfname rc")
      if (isroot) call make_restart_fixture(fname, complete=part_complete(n))
      call get_file_state(comm, isroot, rootpe, fname, nlen=nlen, rc=rc)
      call assert_equal(0, rc, "Restart multi-file partial: get_file_state rc")
      if (part_complete(n)) then
        call assert_equal(1, nlen, "Restart multi-file partial: expected-complete part nlen should be 1")
      else
        call assert_equal(0, nlen, "Restart multi-file partial: expected-incomplete part nlen should be 0")
      end if

      alldone(n) = file_is_complete(comm, isroot, rootpe, fname, .false., 0, rc)
      if (isroot .and. verbose) print '(A,i6,A,L2)',trim(fname)//',  nlen = ',nlen,', complete = ',alldone(n)
      call assert_equal(0, rc, "Restart multi-file partial: file_is_complete rc")
    end do

    call assert_false(all(alldone), "Restart multi-file partial: allDone should be false (one part still incomplete)")
  end subroutine check_restart_multiple_files_partial
  ! --- Assertion helpers ---

  !> Asserts a logical condition is true
  !!
  !! @param[in] condition  The boolean condition to evaluate
  !! @param[in] msg        Message printed if assertion fails
  subroutine assert_true(condition, msg)
    logical,          intent(in) :: condition
    character(len=*), intent(in) :: msg
    if (.not. condition .and. isroot) then
      print *, "  -> ASSERTION FAILED: ", trim(msg)
      total_errors = total_errors + 1
    end if
  end subroutine assert_true

  !> Asserts a logical condition is false
  !!
  !! @param[in] condition  The boolean condition to evaluate
  !! @param[in] msg        Message printed if assertion fails
  subroutine assert_false(condition, msg)
    logical,          intent(in) :: condition
    character(len=*), intent(in) :: msg
    if (condition .and. isroot) then
      print *, "  -> ASSERTION FAILED: ", trim(msg)
      total_errors = total_errors + 1
    end if
  end subroutine assert_false

  !> Asserts equality between two integers
  !!
  !! @param[in] expected  The expected integer value
  !! @param[in] actual    The actual computed integer value
  !! @param[in] msg       Message printed if assertion fails
  subroutine assert_equal(expected, actual, msg)
    integer,          intent(in) :: expected, actual
    character(len=*), intent(in) :: msg
    if (expected /= actual .and. isroot) then
      print *, "  -> ASSERTION FAILED: ", trim(msg), " (expected ", expected, ", got ", actual, ")"
      total_errors = total_errors + 1
    end if
  end subroutine assert_equal

end program test_outputlog_completion
