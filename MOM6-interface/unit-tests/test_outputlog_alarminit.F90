!> @file test_outputlog_alarminit.F90
!> @brief Test code for outputlog Alarm Initialization
!!
!! Probes the use of AlarmInit from MOM6 NUOPC cap to intialize alarms for
!! use by outputlog feature. Tests to ensure the outputlog alarms will trigger
!! for start hours which are not multiples of 6 (eg IAU use cases)
!!
!> @date 08-01-2026

!> Main program for testing outputlog alarm initialization
program test_outputlog_alarminit

  use test_utils

  use ESMF,  only : ESMF_Initialize, ESMF_Finalize, ESMF_SUCCESS, ESMF_FAILURE
  use ESMF,  only : ESMF_CALKIND_GREGORIAN, ESMF_Calendar, ESMF_CalendarCreate, ESMF_Alarm
  use ESMF,  only : ESMF_Clock, ESMF_ClockCreate, ESMF_ClockGet, ESMF_ClockAdvance
  use ESMF,  only : ESMF_Time, ESMF_TimeSet, ESMF_TimeGet, ESMF_TimeInterval, ESMF_TimeIntervalSet
  use ESMF,  only : operator(==), operator(/=), operator(+), operator(-), operator(*)

  use MOM_cap_time,          only : AlarmInit
  use mom_outputlog_methods, only : set_toffset

  implicit none

  integer, parameter :: base_yy = 2021, base_mm = 3, base_dd = 22
  integer, parameter :: maxtests = 25

  character(len=128) :: testname
  character(len=256) :: errmsg
  character(len=256) :: assertmsg
  character(len=20)  :: subname = 'test_alarminit'

  type(testsummary)  :: alarmtests

  logical :: is_passing, assertrc
  integer :: teststart, testfreq
  integer :: rc,nt,n,ierr
  integer :: toffset
  type(ESMF_Time) :: ringTime, startTime, expectedTime
  ! debug printing
  logical :: verbose = .false.

  ! initialize test tracker
  call alarmtests%init(maxtests)

  call ESMF_Initialize(defaultCalKind=ESMF_CALKIND_GREGORIAN, rc=rc)
  call esmf_err(rc, subname, "ESMF_Initialize")

  nt = 0
  ! ===========================================================================
  ! test set_toffset directly
  ! ===========================================================================

  nt = nt + 1
  teststart = 21; testfreq = 6
  write(testname,'(3(A,I2.2))')'test ',nt,' set_toffset: start_hour ',teststart,' freq ',testfreq
  is_passing = (set_toffset(teststart, testfreq) == 3)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), 'expected toffset=3')

  nt = nt + 1
  teststart = 15; testfreq = 6
  write(testname,'(3(A,I2.2))')'test ',nt,' set_toffset: start_hour ',teststart,' freq ',testfreq
  is_passing = (set_toffset(teststart, testfreq) == 3)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), 'expected toffset=3')

  nt = nt + 1
  teststart = 11; testfreq = 6
  write(testname,'(3(A,I2.2))')'test ',nt,' set_toffset: start_hour ',teststart,' freq ',testfreq
  is_passing = (set_toffset(teststart, testfreq) == 1)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), 'expected toffset=1')

  nt = nt + 1
  teststart = 10; testfreq = 3
  write(testname,'(3(A,I2.2))')'test ',nt,' set_toffset: start_hour ',teststart,' freq ',testfreq
  is_passing = (set_toffset(teststart, testfreq) == 2)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), 'expected toffset=2')

  ! already-aligned cases: mod(hour,freq)==0, toffset must stay 0
  nt = nt + 1
  teststart = 9; testfreq = 3
  write(testname,'(3(A,I2.2))')'test ',nt,' set_toffset: start_hour ',teststart,' freq ',testfreq
  is_passing = (set_toffset(teststart, testfreq) == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), 'expected toffset=0 (already aligned)')

  nt = nt + 1
  teststart = 6; testfreq = 6
  write(testname,'(3(A,I2.2))')'test ',nt,' set_toffset: start_hour ',teststart,' freq ',testfreq
  is_passing = (set_toffset(teststart, testfreq) == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), 'expected toffset=0 (already aligned)')

  ! freqs = 1 and 24: explicitly toffset zero regardless of start hour
  nt = nt + 1
  teststart = 5; testfreq = 1
  write(testname,'(3(A,I2.2))')'test ',nt,' set_toffset: start_hour ',teststart,' freq ',testfreq
  is_passing = (set_toffset(teststart, testfreq) == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), 'expected toffset=0 (freq=1 always excluded)')

  nt = nt + 1
  teststart = 5; testfreq = 24
  write(testname,'(3(A,I2.2))')'test ',nt,' set_toffset: start_hour ',teststart,' freq ',testfreq
  is_passing = (set_toffset(teststart, testfreq) == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), 'expected toffset=0 (freq=24 always excluded)')

  ! ===========================================================================
  ! test capture of ringtime via test by setting a small dt so that the alarm
  ! never rings
  ! ===========================================================================

  nt = nt + 1
  teststart = 0; testfreq = 6
  write(testname,'(3(A,I2.2))')'test ',nt,' start_hour ',teststart,' freq ',testfreq
  call run_case(testfreq, teststart, ierr, errmsg, dt=60)   ! 200*60s ~ 3.3h, well short of the 6h needed

  is_passing = (ierr /= 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ===========================================================================
  ! test alignment for start hours, freqs=6,24; ring hour must land on a
  ! multiple of freq
  ! ===========================================================================

  testfreq = 6
  do n = 1,8
     nt = nt + 1
     teststart = (n-1)*3
     write(testname,'(3(A,I2.2))')'test ',nt,' start_hour ',teststart,' freq ',testfreq
     call run_case(testfreq, teststart, ierr, errmsg)

     is_passing = (ierr == 0)
     call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
     call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))
  enddo

  ! ------------------
  testfreq = 24
  do n = 1,8
     nt = nt + 1
     teststart = (n-1)*3
     write(testname,'(3(A,I2.2))')'test ',nt,' start_hour ',teststart,' freq ',testfreq
     call run_case(testfreq, teststart, ierr, errmsg)

     is_passing = (ierr == 0)
     call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
     call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))
  enddo

  ! ===========================================================================
  ! test alignment for freq=1,3 under the generalized per-frequency toffset:
  ! ring hour must land on a multiple of freq
  ! ===========================================================================

  nt = nt + 1
  teststart = 0; testfreq = 1
  write(testname,'(3(A,I2.2))')'test ',nt,' start_hour ',teststart,' freq ',testfreq
  call run_case(testfreq, teststart, ierr, errmsg)

  is_passing = (ierr == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  nt = nt + 1
  teststart = 9; testfreq = 1
  write(testname,'(3(A,I2.2))')'test ',nt,' start_hour ',teststart,' freq ',testfreq
  call run_case(testfreq, teststart, ierr, errmsg)

  is_passing = (ierr == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  nt = nt + 1
  teststart = 0; testfreq = 3
  write(testname,'(3(A,I2.2))')'test ',nt,' start_hour ',teststart,' freq ',testfreq
  call run_case(testfreq, teststart, ierr, errmsg)

  is_passing = (ierr == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  nt = nt + 1
  teststart = 9; testfreq = 3
  write(testname,'(3(A,I2.2))')'test ',nt,' start_hour ',teststart,' freq ',testfreq
  call run_case(testfreq, teststart, ierr, errmsg)

  is_passing = (ierr == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  ! start=21, freq=3 crosses midnight (21+3=24 -> hour 0, next day)
  nt = nt + 1
  teststart = 21; testfreq = 3
  write(testname,'(3(A,I2.2))')'test ',nt,' start_hour ',teststart,' freq ',testfreq
  call run_case(testfreq, teststart, ierr, errmsg)

  is_passing = (ierr == 0)
  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ===========================================================================
  ! Literal-value tests: hand-verified expected ring times, independent of
  ! any re-derived formula.
  ! ===========================================================================

  ! IAU case: start 21 is one interval before the nominal grid crosses midnight.
  ! The averaging window centered here starts the following day.
  nt = nt + 1
  teststart = 21; testfreq = 6
  write(testname,'(3(A,I2.2))')'test ',nt,' literal check: start_hour ',teststart,' freq ',testfreq

  call find_ring_time(testfreq, teststart, ringTime, startTime, toffset, ierr, errmsg)
  if (ierr == 0) then
     call ESMF_TimeSet(expectedTime, yy=base_yy, mm=base_mm, dd=base_dd+1, h=0, rc=rc)
     call esmf_err(rc, subname, "ESMF_TimeSet(expectedTime)")
     is_passing = (ringTime == expectedTime)
  else
     is_passing = .false.
  endif

  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  ! IAU case, same-day: 3h before the nominal 18:00.
  nt = nt + 1
  teststart = 15; testfreq = 6
  write(testname,'(3(A,I2.2))')'test ',nt,' literal check: start_hour ',teststart,' freq ',testfreq

  call find_ring_time(testfreq, teststart, ringTime, startTime, toffset, ierr, errmsg)
  if (ierr == 0) then
     call ESMF_TimeSet(expectedTime, yy=base_yy, mm=base_mm, dd=base_dd, h=18, rc=rc)
     call esmf_err(rc, subname, "ESMF_TimeSet(expectedTime)")
     is_passing = (ringTime == expectedTime)
  else
     is_passing = .false.
  endif

  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  ! Arbitrary, non-IAU, non-multiple-of-3 start hour
  nt = nt + 1
  teststart = 11; testfreq = 6
  write(testname,'(3(A,I2.2))')'test ',nt,' literal check: start_hour ',teststart,' freq ',testfreq

  call find_ring_time(testfreq, teststart, ringTime, startTime, toffset, ierr, errmsg)
  if (ierr == 0) then
     call ESMF_TimeSet(expectedTime, yy=base_yy, mm=base_mm, dd=base_dd, h=12, rc=rc)
     call esmf_err(rc, subname, "ESMF_TimeSet(expectedTime)")
     is_passing = (ringTime == expectedTime)
  else
     is_passing = .false.
  endif

  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  ! Non-IAU restart scenario: start=10 is not a multiple of 3, so the
  ! generalized per-frequency toffset must apply here too. 12:15 is the
  ! first usable 3-hourly average, so the ring must land at hour 12.
  nt = nt + 1
  teststart = 10; testfreq = 3
  write(testname,'(3(A,I2.2))')'test ',nt,' literal check: start_hour ',teststart,' freq ',testfreq

  call find_ring_time(testfreq, teststart, ringTime, startTime, toffset, ierr, errmsg)
  if (ierr == 0) then
     call ESMF_TimeSet(expectedTime, yy=base_yy, mm=base_mm, dd=base_dd, h=12, rc=rc)
     call esmf_err(rc, subname, "ESMF_TimeSet(expectedTime)")
     is_passing = (ringTime == expectedTime)
  else
     is_passing = .false.
  endif

  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  ! freq=24, rRing must land at exactly start+24h, same as if toffset were
  ! always 0, regardless of start_hour.
  nt = nt + 1
  teststart = 10; testfreq = 24
  write(testname,'(3(A,I2.2))')'test ',nt,' literal check: start_hour ',teststart,' freq ',testfreq

  call find_ring_time(testfreq, teststart, ringTime, startTime, toffset, ierr, errmsg)
  if (ierr == 0) then
     call ESMF_TimeSet(expectedTime, yy=base_yy, mm=base_mm, dd=base_dd+1, h=teststart, rc=rc)
     call esmf_err(rc, subname, "ESMF_TimeSet(expectedTime)")
     is_passing = (ringTime == expectedTime)
  else
     is_passing = .false.
  endif

  call assert_equal(is_passing, .true., testname, assertrc, assertmsg)
  call addresult(alarmtests, assertrc, trim(assertmsg), trim(errmsg))

  ! ------------------
  ! Test results
  ! ------------------

  if (alarmtests%nfail > 0) then
     print '(A)', 'FAIL: At least one test failed '
     do n = 1,alarmtests%count
        if (.not. alarmtests%teststatus(n)) print '(A)', trim(alarmtests%testmessage(n)%str)//'  [' &
             //trim(alarmtests%errmessage(n)%str)//']'
     enddo
     stop 1
  else
     do n = 1,alarmtests%count
        if (verbose .and. len_trim(alarmtests%errmessage(n)%str) > 0) then
           print '(A)', trim(alarmtests%testmessage(n)%str)//'  ['//trim(alarmtests%errmessage(n)%str)//']'
        else
           print '(A)', trim(alarmtests%testmessage(n)%str)
        endif
     enddo
  endif
  print '(3(A,I0))','Total tests = ',alarmtests%count,' Passing = ',alarmtests%npass,' Failing = ',alarmtests%nfail

  call ESMF_Finalize(rc=rc)
  call esmf_err(rc, subname, "ESMF_Finalize")

contains
  !> Steps the clock forward (small dt, checked every step) until the given alarm actually rings,
  !! or gives up after max_steps.
  !!
  !! @param[in]      freq          output frequency (hours)
  !! @param[in]      start_hour    start_hour for clock
  !! @param[out]     ringTime      ESMF_Time at ring
  !! @param[out]     startTime     ESMF_Time at start
  !! @param[out]     toffset       required hour offset to set ring at valid interval
  !! @param[out]     ierr          error return code
  !! @param[out]     errmsg        returned error message
  !! @param[in]      dt            optional timestep used to test no-ring
  subroutine find_ring_time(freq, start_hour, ringTime, startTime, toffset, ierr, errmsg, dt)

    integer,           intent(in)  :: freq, start_hour
    type(ESMF_Time),   intent(out) :: ringTime
    type(ESMF_Time),   intent(out) :: startTime
    integer,           intent(out) :: toffset
    integer,           intent(out) :: ierr
    character(len=*),  intent(out) :: errmsg
    integer, optional, intent(in)  :: dt

    type(ESMF_Clock)        :: clock
    type(ESMF_Calendar)     :: cal
    type(ESMF_Time)         :: refTime
    type(ESMF_TimeInterval) :: timeStep, tincrement, alarmoffset
    type(ESMF_Alarm)        :: alarm

    integer :: rc, use_dt
    integer :: rcnt, step
    integer, parameter :: max_steps = 200   ! fixed everywhere -- covers 100h at the default dt=1800s;
                                            ! only dt varies per-case to change effective coverage
    logical :: rang

    ierr = 0
    errmsg = ''

    use_dt = 1800
    if (present(dt)) use_dt = dt

    cal = ESMF_CalendarCreate(ESMF_CALKIND_GREGORIAN, rc=rc)
    call esmf_err(rc, subname,  "ESMF_CalendarCreate")
    call ESMF_TimeSet(startTime, yy=base_yy, mm=base_mm, dd=base_dd, h=start_hour, calendar=cal, rc=rc)
    call esmf_err(rc, subname,  "ESMF_TimeSet(startTime)")
    call ESMF_TimeIntervalSet(timeStep, s=use_dt, rc=rc)
    call esmf_err(rc, subname,  "ESMF_TimeIntervalSet(timeStep)")

    clock = ESMF_ClockCreate(timeStep=timeStep, startTime=startTime, rc=rc)
    call esmf_err(rc, subname,  "ESMF_ClockCreate")
    call ESMF_TimeIntervalSet(tincrement, m=1, rc=rc)
    call esmf_err(rc, subname,  "ESMF_TimeIntervalSet(tincrement)")

    toffset = set_toffset(start_hour, freq)

    alarmoffset = toffset*60*tincrement
    refTime = startTime + alarmoffset

    call AlarmInit(clock,      &
         alarm     = alarm,    &
         option    = 'nhours', &
         opt_n     = freq,     &
         opt_ymd   = -999,     &
         RefTime   = refTime,  &
         alarmname = 'test_alarm', rc=rc)
    call esmf_err(rc, subname,  "AlarmInit")

    ! Find first ring
    rang = .false.
    do step = 1, max_steps
       call ESMF_ClockAdvance(clock, ringingAlarmCount=rcnt, rc=rc)
       call esmf_err(rc, subname,  "ESMF_ClockAdvance")
       if (rcnt > 0) then
          rang = .true.
          exit
       endif
    end do
    if (.not. rang) then
       ierr = 1
       errmsg = 'ERROR: alarm never rang within '//itoa(max_steps)//' steps (dt='//itoa(use_dt)//'s)'
       return
    end if

    if (verbose) then
       print '(A,I0,A)', "  alarm first rang after ", step, " steps of stepping the clock forward"
       print '(A,I0,A,/)', "  RefTime passed to AlarmInit is shifted by that same ", toffset, "h"
    endif

    call ESMF_ClockGet(clock, currTime=ringTime, rc=rc)
    call esmf_err(rc, subname,  "ESMF_ClockGet(currTime at ringTime)")
  end subroutine find_ring_time
  !> Runs one freq/start_hour case end to end. Returns ierr, where ierr==0 means the
  !! alarm rang in time and both the primary (structural) and secondary (regression)
  !! checks passed
  !!
  !! @param[in]      freq          output frequency (hours)
  !! @param[in]      start_hour    start_hour for clock
  !! @param[out]     ierr          error return code
  !! @param[out]     errmsg        returned error message
  !! @param[in]      dt            optional timestep used to test no-ring
  subroutine run_case(freq, start_hour, ierr, errmsg, dt)

    integer,           intent(in)  :: freq, start_hour
    integer,           intent(out) :: ierr
    character(len=*),  intent(out) :: errmsg
    integer, optional, intent(in)  :: dt

    type(ESMF_Time)         :: startTime, ringTime, regressionTime, expectedTime
    type(ESMF_TimeInterval) :: regressionInterval, freqInterval

    integer :: rc
    integer :: toffset, ring_day, ring_hour
    logical :: primary_ok, secondary_ok

    call find_ring_time(freq, start_hour, ringTime, startTime, toffset, ierr, errmsg, dt)
    if (ierr /= 0) return   ! never rang

    call ESMF_TimeGet(ringTime, dd=ring_day, h=ring_hour, rc=rc)
    call esmf_err(rc, subname,  "ESMF_TimeGet(ringTime)")

    ! --- primary (structural) check. verify that the ringhour is a multiple of freq
    if (freq == 24) then
       call ESMF_TimeIntervalSet(freqInterval, h=freq, rc=rc)
       call esmf_err(rc, subname, "ESMF_TimeIntervalSet(freqInterval)")
       expectedTime = startTime + freqInterval
       primary_ok = (ringTime == expectedTime)
       if (.not. primary_ok) then
          errmsg = trim(errmsg)//'PRIMARY FAIL: ring did not occur at exactly start+24h. '
       endif
    else
       primary_ok = (mod(ring_hour, freq) == 0)
       if (.not. primary_ok) then
          errmsg = trim(errmsg)//'PRIMARY FAIL: ring hour '//itoa(ring_hour)//' is not a multiple of freq='// &
               itoa(freq)//'. '
       endif
    endif

    ! --- secondary (regression) check. verify manual construction of timeinterval based
    ! on toffset rings at the right time
    call ESMF_TimeIntervalSet(regressionInterval, h=predicted_ring_offset(freq, toffset), rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeIntervalSet(regressionInterval)")
    regressionTime = startTime + regressionInterval
    secondary_ok = (ringTime == regressionTime)
    if (.not. secondary_ok) then
       errmsg = trim(errmsg)//'SECONDARY FAIL: ring time did not match the hand-derived regression value. '
    endif

    if (.not. (primary_ok .and. secondary_ok)) ierr = 1
  end subroutine run_case
  !> Independent re-derivation of AlarmInit's rewind-then-advance loop, for
  !! the regression check only. Returns the elapsed-hours offset from startTime
  !! to construct ESMF_TimeInterval to compare full ESMF_Time objects
  !!
  !! @param[in]    freq     alarm freq (hour)
  !! @param[in]    toffset  establised time offset (hour)
  function predicted_ring_offset(freq, toffset) result(val)
    integer, intent(in) :: freq, toffset
    integer :: val

    val = toffset - freq
    do while (val <= 0)
       val = val + freq
    end do
  end function predicted_ring_offset
end program test_outputlog_alarminit
