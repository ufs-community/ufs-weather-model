!> @file test_outputlog_freqn.F90
!> @brief Orchestration test for the track_freqn machinery
!!
!! Set up a modelClock and ring state for a logging frequency; calls into track_freqn
!! with fixed curr/next Time, mimicking actual ModelAdvance where the clock itself never
!! advances. Creates and completes files in same sequence as actual production code.
!!
!> @date 08-12-2026

!> Main program for testing outputlog_freqn file tracking
program test_outputlog_freqn

  use ESMF
  use mpi_f08,                only : MPI_Init, MPI_Comm, MPI_Comm_rank, MPI_COMM_WORLD, MPI_Barrier
  use test_utils
  use mom_cap_outputlog,      only : track_freqn
  use mom_outputlog_methods,  only : outputlog_config_type, outputlog_state_type, outputlog_modeltime_type
  use mom_outputlog_methods,  only : get_timestr, get_importexport, set_toffset, get_file_state, debug_info
  use test_helpers,           only : base_yy, base_mm, base_dd, setup_case, handlefiles, setup_restarttimes
  use test_helpers,           only : setup_expected_lastrestart_times

  implicit none

  integer, parameter :: maxtests = 10

  type(MPI_Comm) :: comm
  integer        :: rank, ierr, rootpe
  logical        :: isroot

  type(ESMF_Time), allocatable :: lastrestart_times(:), expected_lastrestarts(:)
  integer,         allocatable :: expected_lastrestart_hours(:)

  character(len=128) :: testname
  character(len=256) :: assertmsg
  character(len=20)  :: subname = 'test_track_freqn'
  character(len=256) :: cmdstr = ''
  character(len=256) :: outputdir = ''

  type(testsummary) :: freqntests

  logical :: debug_onroot
  logical :: assertrc
  integer :: n,nt
  integer :: expected, completions, nmatches
  logical :: verbose = .false.

  character(len=16) :: timestr
  integer :: rc

  comm = MPI_COMM_WORLD
  rootpe = 0
  outputdir = "./"

  call freqntests%init(maxtests)

  call MPI_Init(ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  isroot = (rank == rootpe)
  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=ierr)
  call esmf_err(ierr, subname, "ESMF_Initialize")

  debug_onroot = verbose .and. isroot
  nt = 0
  ! ===========================================================================
  ! Test cases
  ! ===========================================================================

  !------------------
  nt = nt + 1
  expected = 1 ! ring at hour=18 completes 09 file
  write(testname,'(A,I2.2,A)')'test ',nt,' two rings, one tracked averaging window completes '

  call run_case(trim(testname),           &
       freq=6, start_hour=6, runhours=13, &
       use_filesize=.true.,               &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  !------------------
  nt = nt + 1
  expected = 5 ! rings at hour=8,9,10 complete 07,08,09 files; first finalize completes file 10, second complete 11 file
  write(testname,'(A,I2.2,A)')'test ',nt,' snapshots, multiple rings complete correctly; 2 at finalize'

  call run_case(trim(testname),                      &
       freq=1, start_hour=6, runhours=5,             &
       use_filesize=.true., timereduce='none',       &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  ! ------------------
  nt = nt + 1
  expected = 2  ! stop at hour=18 completes 09 and 15 file
  write(testname,'(A,I2.2,A)')'test ',nt,' two tracked averaging windows, both complete at finalize'

  call run_case(trim(testname),           &
       freq=6, start_hour=6, runhours=12, &
       use_filesize=.true.,               &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  ! ------------------
  nt = nt + 1
  expected = 1 ! ring at hour 18 tracks 09 file for DATM
  write(testname,'(A,I2.2,A)')'test ',nt,' two rings, one tracked averaging window completes (DATM)'

  call run_case(trim(testname),               &
       freq=6, start_hour=6, runhours=13,     &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  ! ------------------
  nt = nt + 1
  expected = 0 ! ring at hour 12 tracks 03 phantom file
  write(testname,'(A,I2.2,A)')'test ',nt,' single ring tracks phantom 03 file'

  call run_case(trim(testname),          &
       freq=6, start_hour=6, runhours=7, &
       use_filesize=.true.,              &
       completions = completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  ! ------------------
  nt = nt + 1
  expected = 2 ! rings at hours 12,18 track phantom files 03,09; hour=24 tracks 15 and 21 file at finalize
  write(testname,'(A,I2.2,A)')'test ',nt,' three rings, but only two real windows, both complete at finalize'

  call run_case(trim(testname),                     &
       freq=6, start_hour=9, runhours=15,           &
       use_filesize=.true.,                         &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  ! ------------------
  nt = nt + 1
  expected = 5 ! rings at hour 15,18,21 track files 10.5,13.5,16.5; hour=24 tracks 19.5 and 22.5 file at finalize
  write(testname,'(A,I2.2,A)')'test ',nt,' three rings, but only two real windows, both complete at finalize'

  call run_case(trim(testname),                     &
       freq=3, start_hour=9, runhours=15,           &
       use_filesize=.true.,                         &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')
  ! ------------------
  nt = nt + 1
  expected = 1
  write(testname,'(A,I2.2,A)')'test ',nt,' two ring, one real window, completes correctly'

  call run_case(trim(testname),               &
       freq=24, start_hour=6, runhours=56,    &
       use_filesize=.true.,                   &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  ! ------------------
  ! io_layout
  ! ------------------
  nt = nt + 1
  expected = 1 ! ring at hour 18 tracks 09 file for DATM, nfiles=4
  write(testname,'(A,I2.2,A)')'test ',nt,' io_layout: nfiles>1 uses .0000 suffix, completes correctly'
  call run_case(trim(testname),               &
       freq=6, start_hour=6, runhours=13,     &
       nfiles=4,                              &
       completions=completions)
  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  ! ===========================================================================
  ! Test cases with restart pairing
  !   - restart_hours are a frequency if a single value is provided
  !   - restart_hours are a set value if specific hours are provided
  ! ===========================================================================

  ! ------------------
  nt = nt + 1
  expected = 6
  write(testname,'(A,I2.2,A)')'test ',nt,' starthour=6, runhours=36, restart freq=15 '
  allocate(lastrestart_times(expected))
  allocate(expected_lastrestarts(expected))
  allocate(expected_lastrestart_hours(expected))

  call run_case(trim(testname), freq=6, start_hour=6, runhours=36, &
       use_filesize              = .true.,                         &
       restart_hours             =[15],                            &
       lastrestart_times         =lastrestart_times,               &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  if (completions == expected) then
     testname = trim(testname)//', restart pairing at file completion'
     nmatches = 0
     expected_lastrestart_hours=[ 6,21,21,36,36,36]
     call setup_expected_lastrestart_times(expected_lastrestart_hours, expected_lastrestarts)
     do n = 1,expected
        if (debug_onroot) then
           timestr = get_timestr(expected_lastrestarts(n),rc=rc)
           call esmf_err(rc, subname, "get expected_lastrestart")
           print '(A)','expected last restart times '//timestr
        endif
        if (lastrestart_times(n) == expected_lastrestarts(n)) then
           nmatches = nmatches + 1
        endif
     enddo

     call assert_equal(nmatches, expected, testname, assertrc, assertmsg)
     call addresult(freqntests, assertrc, trim(assertmsg), '')
  endif
  deallocate(lastrestart_times, expected_lastrestarts, expected_lastrestart_hours)

  ! ------------------
  nt = nt + 1
  write(testname,'(A,I2.2,A)')'test ',nt,' starthour=6, runhours=30, restart freq=3 '
  expected = 5 ! expected time matches
  allocate(lastrestart_times(expected))
  allocate(expected_lastrestarts(expected))
  allocate(expected_lastrestart_hours(expected))

  call run_case(trim(testname), freq=6, start_hour=6, runhours=30, &
       use_filesize              = .true.,                         &
       restart_hours             =[3],                             &
       lastrestart_times         =lastrestart_times,               &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  if (completions == expected) then
     testname = trim(testname)//', restart pairing at file completion'
     nmatches = 0
     expected_lastrestart_hours=[18,24,30,36,36]
     call setup_expected_lastrestart_times(expected_lastrestart_hours, expected_lastrestarts)
     do n = 1,expected
        if (debug_onroot) then
           timestr = get_timestr(expected_lastrestarts(n),rc=rc)
           call esmf_err(rc, subname, "get expected_lastrestart")
           print '(A)','expected last restart times '//timestr
        endif
        if (lastrestart_times(n) == expected_lastrestarts(n)) then
           nmatches = nmatches + 1
        endif
     enddo

     call assert_equal(nmatches, expected, testname, assertrc, assertmsg)
     call addresult(freqntests, assertrc, trim(assertmsg), '')
  endif
  deallocate(lastrestart_times, expected_lastrestarts, expected_lastrestart_hours)

  ! ------------------
  nt = nt + 1
  write(testname,'(A,I2.2,A)')'test ',nt,' starthour=9, runhours=96, restart hours specified '
  expected = 15 ! expected time matches
  allocate(lastrestart_times(expected))
  allocate(expected_lastrestarts(expected))
  allocate(expected_lastrestart_hours(expected))

  call run_case(trim(testname), freq=6, start_hour=9, runhours=93, &
       use_filesize              = .true.,                         &
       restart_hours             =[12,30,51,84],                   &
       lastrestart_times         =lastrestart_times,               &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  if (completions == expected) then
     testname = trim(testname)//', restart pairing at file completion'
     nmatches = 0
     expected_lastrestart_hours=[12,30,30,30,30,51,51,51,51,51,84,84,84,84,84]
     call setup_expected_lastrestart_times(expected_lastrestart_hours, expected_lastrestarts)
     do n = 1,expected
        if (debug_onroot) then
           timestr = get_timestr(expected_lastrestarts(n),rc=rc)
           call esmf_err(rc, subname, "get expected_lastrestart")
           print '(A)','expected last restart times '//timestr
        endif
        if (lastrestart_times(n) == expected_lastrestarts(n)) then
           nmatches = nmatches + 1
        endif
     enddo

     call assert_equal(nmatches, expected, testname, assertrc, assertmsg)
     call addresult(freqntests, assertrc, trim(assertmsg), '')
  endif
  deallocate(lastrestart_times, expected_lastrestarts, expected_lastrestart_hours)

  ! ------------------
  nt = nt + 1
  expected = 4 ! rings at hour=8,9 complete 07,08 files; first finalize completes 09, second complete 10; restart freq=1
  write(testname,'(A,I2.2,A)')'test ',nt,' snapshots, multiple rings complete correctly; 2 at finalize'
  allocate(lastrestart_times(expected))
  allocate(expected_lastrestarts(expected))
  allocate(expected_lastrestart_hours(expected))

  call run_case(trim(testname), freq=1, start_hour=6, runhours=4, &
       use_filesize             =.true.,                          &
       timereduce               ='none',                          &
       restart_hours            =[1],                             &
       lastrestart_times        =lastrestart_times,               &
       completions=completions)

  call assert_equal(completions, expected, testname, assertrc, assertmsg)
  call addresult(freqntests, assertrc, trim(assertmsg), '')

  if (completions == expected) then
     testname = trim(testname)//', restart pairing at file completion'
     nmatches = 0
     expected_lastrestart_hours=[8,9,10,10]
     call setup_expected_lastrestart_times(expected_lastrestart_hours, expected_lastrestarts)
     do n = 1,expected
        if (debug_onroot) then
           timestr = get_timestr(expected_lastrestarts(n),rc=rc)
           call esmf_err(rc, subname, "get expected_lastrestart")
           print '(A)','expected last restart times '//timestr
        endif
        if (lastrestart_times(n) == expected_lastrestarts(n)) then
           nmatches = nmatches + 1
        endif
     enddo

     call assert_equal(nmatches, expected, testname, assertrc, assertmsg)
     call addresult(freqntests, assertrc, trim(assertmsg), '')
  endif
  deallocate(lastrestart_times, expected_lastrestarts, expected_lastrestart_hours)

  ! ------------------
  ! Test results
  ! ------------------
  if (isroot) then
  if (freqntests%nfail > 0) then
     print '(A)', 'FAIL: At least one test failed '
     do n = 1,freqntests%count
        if (.not. freqntests%teststatus(n)) print '(A)', trim(freqntests%testmessage(n)%str)
     enddo
  else
     do n = 1,freqntests%count
        print '(A)', trim(freqntests%testmessage(n)%str)
     enddo
  endif
  print '(3(A,I0))','Total tests = ',freqntests%count,' Passing = ',freqntests%npass,' Failing = ',freqntests%nfail
  endif

  ! cleanup if running locally
  if (isroot) then
    cmdstr = 'rm -f '//trim(outputdir)//'*.nc '//trim(outputdir)//'*.mom6.*'//'  ./PET*'
    call execute_command_line(trim(cmdstr), wait=.true.)
  endif

  call ESMF_Finalize(rc=ierr)
  call esmf_err(ierr, subname, "ESMF_Finalize")

  if (freqntests%nfail > 0) then
     if (isroot) print '(A)','Test failures! '
     stop
  endif

contains
  !> Run a single test case through a simulated modelClock advance cycle
  !!
  !! @param[in]      test               descriptive test name
  !! @param[in]      freq               output frequency
  !! @param[in]      start_hour         clock start time
  !! @param[in]      runhours           clock advance length
  !! @param[in]      timereduce         optional, to specify average or snapshot mode
  !! @param[in]      use_filesize       optional, to specify completion type
  !! @param[in]      nfiles             optional, to specify io-layout file number
  !! @param[in]      restart_hours      optional, restart hours or freq
  !! @param[out]     lastrestart_times  optional, last restart available at file completion
  !! @param[out]     completions        number of file completions found
  subroutine run_case(test, freq, start_hour, runhours, timereduce, use_filesize, nfiles, &
       restart_hours, lastrestart_times, completions)

    character(len=*), intent(in)            :: test
    integer,          intent(in)            :: freq, start_hour, runhours
    character(len=*), intent(in), optional  :: timereduce
    logical,          intent(in), optional  :: use_filesize
    integer,          intent(in), optional  :: nfiles
    integer,          intent(in), optional  :: restart_hours(:)
    type(ESMF_Time),  intent(out), optional :: lastrestart_times(:)
    integer,          intent(out)           :: completions

    character(len=7) :: l_timereduce
    logical          :: l_use_filesize
    integer          :: l_nfiles

    type(ESMF_Clock)             :: modelClock
    type(ESMF_Time)              :: stopTime, lastrestart
    type(ESMF_TimeInterval)      :: elapsedtime
    type(ESMF_Time), allocatable :: restart_times(:)

    type(outputlog_config_type)    :: cf_n
    type(outputlog_state_type)     :: state_n
    type(outputlog_modeltime_type) :: modeltime

    integer :: ierr, rc
    integer :: toffset, count
    integer :: minutes, elapsedhours, n_restarts

    logical :: phantom_file, lstop
    logical :: found_firstcompletion = .false.   ! count only the first time the file completes
    logical :: pending = .false.

    character(len=16)  :: timestr
    character(len=40)  :: importexport
    character(len=20)  :: subname = 'run_case'

    !debug
    integer :: nlen, fsize

    completions = 0

    l_timereduce = 'average'
    if (present(timereduce)) l_timereduce = timereduce
    l_use_filesize = .false.
    if (present(use_filesize)) l_use_filesize = use_filesize
    l_nfiles = 1
    if (present(nfiles)) l_nfiles = nfiles

    if (isroot) then
       cmdstr = 'rm -f '//trim(outputdir)//'*.nc '//trim(outputdir)//'*.mom6.*'//'  ./PET*'
       call execute_command_line(trim(cmdstr), wait=.true.)
    endif
    call MPI_Barrier(comm, ierr)
    if (ierr /= 0) then
       write(0,'(A)') "FATAL ("//trim(subname)//"): MPI_Barrier (post-cleanup) failed"
       stop 99
    endif

    ! mimic outputlog_init setup
    call setup_case(start_hour, runhours, freq, l_nfiles, l_timereduce, debug_onroot, &
         modelClock, cf_n, state_n, rc)
    call esmf_err(rc, subname, "setup_case")

    completions = 0
    if (present(restart_hours)) then
       if (size(restart_hours) == 1) then          ! restart_hours are uniform cadence at given frequency
          n_restarts = runhours/restart_hours(1)
       else
          n_restarts = size(restart_hours)
       endif
       allocate(restart_times(n_restarts))
       call setup_restarttimes(start_hour,n_restarts,restart_hours,restart_times,rc)
       call esmf_err(rc, subname, "setup_restarttimes")

       if (debug_onroot) then
          do n = 1, size(restart_times)
             timestr = get_timestr(restart_times(n), rc=rc)
             call esmf_err(rc, subname, "get restart_times")
             print '(A,i3,A)','Restart time ',n,' defined at '//timestr
          enddo
       endif
    endif
    call ESMF_ClockGet(modelClock, startTime=modeltime%startTime, rc=rc)
    call esmf_err(rc, subname, "get start and currTime")
    lastrestart = modeltime%startTime

    if (debug_onroot) then
       print '(A)','Running test '//test
    endif
    call ESMF_TimeIntervalSet(modeltime%tincrement, m=1, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeIntervalSet(tincrement)")

    do while (.not. ESMF_ClockIsStopTime(modelClock, rc=rc))
       call esmf_err(rc, subname, "advance to stop time")

       call ESMF_ClockGet(modelClock, startTime=modeltime%startTime, currTime=modeltime%currTime, rc=rc)
       call esmf_err(rc, subname, "get start and currTime")
       call ESMF_ClockGetNextTime(modelClock, modeltime%nextTime, rc=rc)
       call esmf_err(rc, subname, "get nextTime")
       importexport = get_importexport(modeltime%currTime, modeltime%nextTime, rc=rc)
       call esmf_err(rc, subname, "get importexport")

       state_n%ringing = .false.
       count = 0
       call ESMF_ClockAdvance(modelClock,ringingalarmcount=count,rc=rc)
       call esmf_err(rc, subname, "ringing alarm count")

       if (count > 0) then
          state_n%ringing = .true.
          call ESMF_AlarmRingerOff(cf_n%alarm, rc=rc)
          call esmf_err(rc, subname, "alarm ringer off")
       endif
       call ESMF_AlarmGet(cf_n%alarm, prevRingTime=state_n%prevring, rc=rc)
       call esmf_err(rc, subname, "get prevRing")

       call ESMF_ClockGet(modelClock, stopTime=stopTime,rc=rc)
       call esmf_err(rc, subname, "get stopTime")
       lstop = (modeltime%nextTime == stopTime)

       if (present(restart_hours)) then
          do n = 1, size(restart_times)
             if (modeltime%nextTime == restart_times(n)) lastrestart = restart_times(n)
          enddo
       endif

       ! ======================================================================
       ! set up file states to mimic FMS
       ! ======================================================================

       ! during continuous polling, file is pending until on first advance on new file window (ringing off)
       if (pending .and. len_trim(state_n%filename)>0) then
          call handlefiles(isroot, state_n%filename, l_use_filesize, 'complete')
          pending = .false.
          if (debug_onroot) then
             call get_file_state(comm, isroot, rootpe, state_n%filename, nlen=nlen, fsize=fsize, rc=rc)
             print '(A,i4,i12,2(A,L))',trim(subname)//' complete file '//state_n%filename//'  '//importexport, &
                 nlen,fsize,' pending ',pending,' ringing ',state_n%ringing
          endif
       endif

       ! during continuous polling, file will be created when ring occurs at tracking window
       if (state_n%ringing) then
          phantom_file = .false.
          found_firstcompletion = .false.
          timestr = get_timestr(modeltime%nextTime - cf_n%filename_fhoffset, rc=rc)
          call esmf_err(rc, subname, "get_timestr(ring filename)")
          if (modeltime%nextTime - cf_n%filename_fhoffset <= modeltime%startTime) phantom_file = .true.

          state_n%filename = trim(outputdir)//trim(cf_n%fnameprefix)//trim(timestr)//'.nc' &
               //trim(cf_n%fnamesuffix)

          if (phantom_file) then
             pending = .false.
             if (debug_onroot) then
                print '(A)',' file '//state_n%filename//' is phantom'
             endif
          else
             call handlefiles(isroot, state_n%filename, l_use_filesize, 'create')
             pending = .true.
             if (debug_onroot) then
               call get_file_state(comm, isroot, rootpe, state_n%filename, nlen=nlen, fsize=fsize, rc=rc)
               print '(A,i4,i12,2(A,L))',trim(subname)//' create file '//state_n%filename//'  '//importexport,  &
                   nlen,fsize,' pending ',pending,' ringing ',state_n%ringing
             endif
          endif
       endif

       ! ======================================================================
       ! end of file preparation for continuous polling
       ! ======================================================================

       call track_freqn(modeltime, cf_n, state_n, comm, isroot, rootpe, outputdir, lastrestart, &
            debug_onroot, .false., rc)
       call esmf_err(rc, subname, "track_freqn (main loop)")
       if (state_n%filecomplete .and. .not.found_firstcompletion) then
          completions = completions + 1
          found_firstcompletion = .true.
          if (present(lastrestart_times)) lastrestart_times(completions) = state_n%time_lastrestart
       endif

       ! first call during finalize; io_infra_end/MOM_infra_end finishes any pending file
       if (modeltime%nextTime == stopTime) then
          if (pending) then
             call handlefiles(isroot, state_n%filename, l_use_filesize, 'complete')
             pending = .false.
             found_firstcompletion = .false.
             if (debug_onroot) then
                call get_file_state(comm, isroot, rootpe, state_n%filename, nlen=nlen, fsize=fsize, rc=rc)
                print '(A,i4,i12,2(A,L))',trim(subname)//' complete file '//state_n%filename//'  '//importexport, &
                    nlen,fsize,' pending ',pending,' ringing ',state_n%ringing
             endif

             state_n%ringing = .false.
             state_n%chkfile_nextAdvance = .true.
             call track_freqn(modeltime, cf_n, state_n, comm, isroot, rootpe, outputdir, lastrestart, &
                  debug_onroot, .false., rc)
             call esmf_err(rc, subname, "track_freqn (plain finalize)")
             if (state_n%filecomplete .and. .not.found_firstcompletion) then
                completions = completions + 1
                found_firstcompletion = .true.
                if (present(lastrestart_times)) lastrestart_times(completions) = state_n%time_lastrestart
             endif
          endif
       endif

       ! ======================================================================
       ! actual feature obtains prevRing from model clock to set file name at lstop
       ! here, must rig the correct filename to stage the file
       ! ======================================================================

       toffset = set_toffset(start_hour, freq)
       elapsedtime = modeltime%nextTime - modeltime%startTime
       call ESMF_TimeIntervalGet(elapsedtime, m=minutes, rc=rc)
       call esmf_err(rc, subname, "get elapsedtime at lstop")
       elapsedhours = toffset + minutes/60

       if (modeltime%nextTime == stopTime) then
          if (mod(elapsedhours,freq) == 0) then
             if (trim(cf_n%timereduce) == 'none') then
                timestr = get_timestr(state_n%prevring, rc=rc)
                call esmf_err(rc, subname, "get_timestr(lstop filename, none)")
             else
                timestr = get_timestr(state_n%prevring-30*cf_n%opt_n*modeltime%tincrement, rc=rc)
                call esmf_err(rc, subname, "get_timestr(lstop filename, average)")
             endif

             state_n%filename = trim(outputdir)//trim(cf_n%fnameprefix)//trim(timestr)//'.nc' &
                  //trim(cf_n%fnamesuffix)
             call handlefiles(isroot, state_n%filename, l_use_filesize, 'create-complete')

             if (debug_onroot) then
                call get_file_state(comm, isroot, rootpe, state_n%filename, nlen=nlen, fsize=fsize, rc=rc)
                print '(A,i4,i12,2(A,L))',trim(subname)//' create-complete file '//state_n%filename//'  '//importexport, &
                    nlen,fsize,' pending ',pending,' ringing ',state_n%ringing
             endif
          endif
          ! ======================================================================
          ! end of file preparation for lstop
          ! ======================================================================

          state_n%ringing = .false.
          state_n%chkfile_nextAdvance = .false.
          found_firstcompletion = .false.
          call track_freqn(modeltime, cf_n, state_n, comm, isroot, rootpe, outputdir, lastrestart, &
               debug_onroot, .true., rc)
          call esmf_err(rc, subname, "track_freqn (lstop)")

          if (state_n%filecomplete .and. .not.found_firstcompletion) then
             completions = completions + 1
             found_firstcompletion = .true.
             if (present(lastrestart_times)) lastrestart_times(completions) = state_n%time_lastrestart
          endif
       endif

    enddo ! do while
    if (debug_onroot) then
       print '(A,/)','Completed test '//test
    endif
  end subroutine run_case

end program test_outputlog_freqn
