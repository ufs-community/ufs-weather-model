!> @file test_helpers.F90
!> @brief Shared setup for outputlog unit tests.
!!
!> @date 08-12-2026
module test_helpers

  use ESMF
  use mom_outputlog_methods, only : outputlog_config_type, outputlog_state_type, outputlog_modeltime_type
  use mom_outputlog_methods, only : setup_freq_config, get_timestr, set_toffset
  use MOM_cap_time,          only : AlarmInit
  use nc_fixture_mod,        only : create_schema, write_record, write_padding, write_bulk_data
  use test_utils,            only : esmf_err

  implicit none
  private

  public :: base_yy, base_mm, base_dd
  public :: setup_case, handlefiles, setup_restarttimes, setup_expected_lastrestart_times

  integer, parameter :: base_yy = 2021   !< a standard start year
  integer, parameter :: base_mm = 3      !< a standard start month
  integer, parameter :: base_dd = 22     !< a standard start day

contains
  !> Build a real ESMF_Clock/alarm plus cf_n/state_n for one test case
  !!
  !! @param[in]     start_hour     model start hour
  !! @param[in]     runhours       total run length in hours
  !! @param[in]     freq           output frequency in hours
  !! @param[in]     l_nfiles       number of IO-layout files (1 = single file)
  !! @param[in]     l_timereduce   'average' or 'none'
  !! @param[in]     debug_onroot   enable verbose setup printing
  !! @param[out]    modelClock     the constructed ESMF_Clock
  !! @param[out]    cf_n           this frequency's config
  !! @param[out]    state_n        this frequency's state
  !! @param[out]    rc             return code
  subroutine setup_case(start_hour, runhours, freq, l_nfiles, l_timereduce, debug_onroot, &
       modelClock, cf_n, state_n, rc)

    integer,                     intent(in)  :: start_hour, runhours, freq, l_nfiles
    character(len=*),            intent(in)  :: l_timereduce
    logical,                     intent(in)  :: debug_onroot
    type(ESMF_Clock),            intent(out) :: modelClock
    type(outputlog_config_type), intent(out) :: cf_n
    type(outputlog_state_type),  intent(out) :: state_n
    integer,                     intent(out) :: rc

    type(outputlog_modeltime_type) :: modeltime
    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep, alarmoffset

    integer :: toffset, hour
    character(len=16)  :: startstr, stopstr
    character(len=120) :: subname = 'setup_case'

    rc = ESMF_SUCCESS

    call ESMF_TimeSet(startTime, yy=base_yy, mm=base_mm, dd=base_dd, h=start_hour, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeSet(startTime)")
    call ESMF_TimeSet(stopTime,  yy=base_yy, mm=base_mm, dd=base_dd, h=start_hour+runhours, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeSet(stopTime)")

    call ESMF_TimeIntervalSet(timeStep, s=1800, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeIntervalSet(timeStep)")
    call ESMF_TimeIntervalSet(modeltime%tincrement, m=1, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeIntervalSet(tincrement)")
    modelClock  = ESMF_ClockCreate(name="Model",timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
    call esmf_err(rc, subname, "ESMF_ClockCreate")

    call ESMF_ClockGet(modelclock, currTime=modeltime%currTime, startTime=startTime, stopTime=stopTime, rc=rc)
    call esmf_err(rc, subname, "ESMF_ClockGet start,stop time")
    startstr = get_timestr(startTime, rc=rc)
    call esmf_err(rc, subname, "get_timestr(startTime)")
    stopstr = get_timestr(stopTime, rc=rc)
    call esmf_err(rc, subname, "get_timestr(stopTime)")
    if (debug_onroot) then
       print '(/,A)','Clock will run from '//startstr//' to '//stopstr
    endif

    ! initialize as in production
    cf_n%opt_n       = freq
    cf_n%requested   = .true.
    cf_n%timereduce  = l_timereduce
    cf_n%fnameprefix = 'ocn_'

    call setup_freq_config(freq, l_nfiles, modeltime, cf_n, state_n, rc)
    call esmf_err(rc, subname, "setup_freq_config")

    call ESMF_TimeGet(modeltime%currTime, h=hour, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeGet(hour)")
    toffset = set_toffset(hour, freq)
    alarmoffset = toffset*60*modeltime%tincrement

    call AlarmInit(modelclock,              &
         alarm     = cf_n%alarm,            &
         option    = 'nhours',              &
         opt_n     = freq,                  &
         opt_ymd   = -999,                  &
         RefTime   = modeltime%currTime+alarmoffset,  &
         alarmname = cf_n%alarm_name, rc=rc)
    call esmf_err(rc, subname, "call AlarmInit")

  end subroutine setup_case
  !> Create/complete a netCDF fixture file, matching the real DATM/ATM
  !! completion contract.
  !!
  !! @param[in]  isroot        .true. on the root PE
  !! @param[in]  fname          file name
  !! @param[in]  use_filesize  .true. for ATM-style (fsize-based) completion,
  !!                           .false. for DATM-style (nlen-based)
  !! @param[in]  mode          'create', 'complete', or 'create-complete'
  subroutine handlefiles(isroot, fname, use_filesize, mode)

    logical,          intent(in)  :: isroot
    character(len=*), intent(in)  :: fname
    logical,          intent(in)  :: use_filesize
    character(len=*), intent(in)  :: mode

    select case (mode)
    case ('create')
       if (isroot) then
          call create_schema(fname)
          if (use_filesize) then
             call write_record(fname)
          else
             call write_padding(fname)
          endif
       endif

    case ('complete')
       if (isroot) then
          if (use_filesize) then
             call write_bulk_data(fname)   ! fsize grows past createsize
          else
             call write_record(fname)      ! nlen 0->1
          endif
       endif

    case('create-complete')
       if (isroot) then
          call create_schema(fname)
          if (use_filesize) then
             call write_record(fname)
             call write_bulk_data(fname)   ! fsize grows past createsize
          else
             call write_padding(fname)
             call write_record(fname)      ! nlen 0->1
          endif
       endif

    case default
       if (isroot) then
          print '(A)',' ERROR: unknown case '
          stop
       endif
    end select

  end subroutine handlefiles
  !> Build a list of restart times from a list of restart hours
  !!
  !! @param[in]     start_hour     clock start hour
  !! @param[in]     n_restarts     number of restart times specified
  !! @param[in]     restart_hours  restart frequency, either cadence or list of hours
  !! @param[out]    restart_times  restart times
  !! @param[out]    rc             return code
  subroutine setup_restarttimes(start_hour, n_restarts, restart_hours, restart_times, rc)

    integer,         intent(in)  :: start_hour
    integer,         intent(in)  :: n_restarts
    integer,         intent(in)  :: restart_hours(:)
    type(ESMF_Time), intent(out) :: restart_times(:)
    integer,         intent(out) :: rc

    integer :: n

    type(ESMF_Time)         :: startTime
    type(ESMF_TimeInterval) :: tincrement
    character(len=120)      :: subname = 'setup_restarttimes'

    rc = ESMF_SUCCESS

    call ESMF_TimeSet(startTime, yy=base_yy, mm=base_mm, dd=base_dd, h=start_hour, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeSet(startTime)")
    call ESMF_TimeIntervalSet(tincrement, m=1, rc=rc)
    call esmf_err(rc, subname, "ESMF_TimeIntervalSet(tincrement)")

    if (size(restart_hours) == 1) then
       restart_times(1) = startTime + restart_hours(1)*60*tincrement
       do n = 2, n_restarts
          restart_times(n) = restart_times(n-1) + restart_hours(1)*60*tincrement
       enddo
    else
       do n = 1, size(restart_hours)
          call ESMF_TimeSet(restart_times(n), yy=base_yy, mm=base_mm, dd=base_dd, h=restart_hours(n), rc=rc)
          call esmf_err(rc, subname, "ESMF_TimeSet(restart_times)")
       enddo
    endif
  end subroutine setup_restarttimes
  !> Build a list of expected lastrestart times
  !!
  !! @param[in]    expected_lastrestart_hours  elapsed hour relative to starting day when restart is written
  !! @param[out]   expected_lastrestarts       expected last restart times
  subroutine setup_expected_lastrestart_times(expected_lastrestart_hours,expected_lastrestarts)

    integer,         intent(in)  :: expected_lastrestart_hours(:)
    type(ESMF_Time), intent(out) :: expected_lastrestarts(:)

    integer :: n, rc
    character(len=120) :: subname = 'setup_lastrestart_times'

    do n = 1,size(expected_lastrestart_hours)
       call ESMF_TimeSet(expected_lastrestarts(n), yy=base_yy, mm=base_mm, dd=base_dd, h=expected_lastrestart_hours(n), rc=rc)
       call esmf_err(rc, subname, "get expected_lastrestart")
    enddo
  end subroutine setup_expected_lastrestart_times

end module test_helpers
