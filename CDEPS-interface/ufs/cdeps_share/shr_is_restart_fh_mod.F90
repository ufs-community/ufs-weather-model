module shr_is_restart_fh_mod

  ! Common methods for components to check if it's time to write forecast hour-based restarts

  !use dshr_methods_mod , only : chkerr
  use ESMF, only : ESMF_ConfigCreate, ESMF_ConfigDestroy, ESMF_ConfigLoadFile, &
                     ESMF_ConfigGetLen, ESMF_ConfigGetAttribute, ESMF_TimePrint, &
                     ESMF_LOGMSG_INFO, ESMF_LogWrite, ESMF_TimeInterval, &
                     ESMF_Time, ESMF_KIND_R8, ESMF_Config, ESMF_Clock, &
                     ESMF_TimeIntervalSet, ESMF_TimePrint, operator(+), operator(==), &
                     ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU

  implicit none
  private

  type :: is_restart_fh_type
    logical :: write_restartfh = .false.
    type(ESMF_Time), allocatable :: restartFhTimes(:)
  end type is_restart_fh_type

  public :: init_is_restart_fh, is_restart_fh, finalize_restart_fh, is_restart_fh_type
  public :: log_restart_fh

contains

  !-----------------------------------------------------------------------
  subroutine init_is_restart_fh(currentTime, dtime, lLog, restartfh_info, key)
    !
    ! !DESCRIPTION:
    ! Process restart_fh attribute from model_configure in UFS
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ESMF_Time), intent(in) :: currentTime
    integer, intent(in)         :: dtime ! time step (s)
    logical, intent(in)         :: lLog ! If true, this task logs restart_fh info
    type(is_restart_fh_type), intent(out) :: restartfh_info !restart_fh info for each task
    character(len=*),intent(in),optional :: key
    !
    ! !LOCAL VARIABLES:
    character(len=256)           :: timestr
    integer                      :: n, nfh, fh_s, rc, nhours_fcst, ntimeout, fhi, ifh
    logical                      :: isPresent
    real(kind=ESMF_KIND_R8), allocatable :: restart_fh(:)
    type(ESMF_TimeInterval)      :: fhInterval
    type(ESMF_Config)            :: CF_mc
    character(len=64)            :: key_name
    !-----------------------------------------------------------------------

    ! Default is to create restart times, but other "keys" could use this routine.
    if (present(key)) then
       key_name = trim(key)//':'
    else
       key_name = 'restart_fh:'
    endif
    
    ! set up Times to write non-interval restarts
    inquire(FILE='model_configure', EXIST=isPresent)
    if (isPresent) then !model_configure exists. this is ufs run
      CF_mc = ESMF_ConfigCreate(rc=rc)
      call ESMF_ConfigLoadFile(config=CF_mc,filename='model_configure' ,rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
      nfh = ESMF_ConfigGetLen(config=CF_mc, label=trim(key_name),rc=rc)

      if (nfh .gt. 0) then
        allocate(restart_fh(1:nfh))
        allocate(restartfh_info%restartFhTimes(1:nfh)) !not deallocated here
        call ESMF_ConfigGetAttribute(CF_mc,valueList=restart_fh,label=trim(key_name), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return

        ! Support case where restart_fh: "fh_interval -1"
        if (nfh == 2 .and. int(restart_fh(2)) == -1) then
           call ESMF_ConfigGetAttribute(CF_mc,value=nhours_fcst,label='nhours_fcst:', rc=rc)
           fhi = int(restart_fh(1))
           ntimeout = nhours_fcst/fhi
           nfh = ntimeout
           deallocate(restart_fh)
           deallocate(restartfh_info%restartFhTimes)
           allocate(restart_fh(1:nfh))
           allocate(restartfh_info%restartFhTimes(1:nfh))
           do ifh = 1,nfh
              restart_fh(ifh) = ifh*fhi
           end do
        endif
        
        ! create a list of times at each restart_fh
        do n = 1,nfh
          fh_s = NINT(3600*restart_fh(n))
          call ESMF_TimeIntervalSet(fhInterval, s=fh_s, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          restartfh_info%restartFhTimes(n) = currentTime + fhInterval
          call ESMF_TimePrint(restartfh_info%restartFhTimes(n), options="string", &
                              preString="restart_fh at ", unit=timestr, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
          if (lLog) then
            if (mod(fh_s,dtime) /= 0) then
              call ESMF_LogWrite('restart time NOT to be written for '//trim(timestr), ESMF_LOGMSG_INFO)
            else
              call ESMF_LogWrite('restart time to be written for '//trim(timestr), ESMF_LOGMSG_INFO)
            end if
          end if
        end do
        deallocate(restart_fh)
      end if !nfh>0
      call ESMF_ConfigDestroy(CF_mc, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    end if !model_configure

  end subroutine init_is_restart_fh

  subroutine is_restart_fh(clock, restartfh_info, lWrite)
    !
    ! !DESCRIPTION:
    ! True/false if time to write restart
    !
    ! !USES:
    use ESMF, only : ESMF_ClockGetNextTime

    !
    ! !ARGUMENTS:
    type(ESMF_Clock), intent(in) :: clock
    type(is_restart_fh_type), intent(inout) :: restartfh_info
    logical, intent(out)         :: lWrite ! time to write?
    !
    ! !LOCAL VARIABLES:
    integer                    :: nfh, rc
    type(ESMF_Time)            :: nextTime
    !-----------------------------------------------------------------------

    restartfh_info%write_restartfh = .false.
    if (allocated(restartfh_info%restartFhTimes)) then
      ! check if next time is == to any restartfhtime
      do nfh = 1,size(restartfh_info%restartFhTimes)
        call ESMF_ClockGetNextTime(clock, nextTime=nexttime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
        if (nextTime == restartfh_info%restartFhTimes(nfh)) restartfh_info%write_restartfh = .true.
      end do
    end if

    lWrite = restartfh_info%write_restartfh

  end subroutine is_restart_fh
  !===============================================================================
  !> Write a log file
  !!
  !> @details Write a log file for a named component when a restart file is written
  !!
  !! @param[in]   nextTime       the ESMF time at the end of a ModelAdvance
  !! @param[in]   startTime      the ESMF time at the Model Start
  !! @param[in]   complog        the named component
  !! @param[in]   prefixtime     optional, if true log filename has time prefix
  !! @param[in]   lastrestart    optional, if present, write the time of the last restart
  !! @param[in]   lastoutput     optional, if present, write the filename written at this FH
  !! @param[out]  rc return code
  !!
  !> @authorDenise.Worthen@noaa.gov
  !> @date 04-14-2025
  subroutine log_restart_fh(myTime, startTime, complog, prefixtime, lastrestart, lastoutput, rc)

    use ESMF,              only : ESMF_SUCCESS, ESMF_MAXSTR, ESMF_Time, ESMF_TimeInterval
    use ESMF,              only : ESMF_TimeGet, ESMF_TimeIntervalGet
    use ESMF,              only : operator(==), operator(-)

    type(ESMF_Time),  intent(in)           :: myTime, startTime
    character(len=*), intent(in)           :: complog
    logical,          intent(in), optional :: prefixtime
    type(ESMF_Time),  intent(in), optional :: lastrestart
    character(len=*), intent(in), optional :: lastoutput
    integer,         intent(out)           :: rc

    ! local variables
    type(ESMF_TimeInterval)     :: elapsedTime
    real(ESMF_KIND_R8)          :: fhour
    character(ESMF_MAXSTR)      :: filename
    character(ESMF_MAXSTR)      :: nexttimestring
    character(ESMF_MAXSTR)      :: timestring
    integer                     :: fh_logunit
    integer                     :: yr,mon,day,hour,minute,sec ! time units
    logical                     :: lprefix
    character(ESMF_MAXSTR)      :: lastout
    character(len=*), parameter :: subname='(log_restart_fh)'
    !-----------------------------------------------------------------------

    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    rc = ESMF_SUCCESS

    lprefix = .false.
    if (present(prefixtime)) then
       lprefix = prefixtime
    end if
    lastout = ''
    if (present(lastoutput)) then
       lastout = trim(lastoutput)
    end if

    elapsedTime = myTime - startTime
    call ESMF_TimeIntervalGet(elapsedTime, h_r8=fhour,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    call ESMF_TimeGet(myTime, yy=yr, mm=mon, dd=day, h=hour, m=minute, s=sec, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    write(nexttimestring,'(6i8)')yr,mon,day,hour,minute,sec
    if (lprefix) then
       write(filename,'(i4.4,2(i2.2),A,3(i2.2),A)') yr, mon, day,'.', hour, minute, sec,'.'//trim(complog)
    else
       write(filename,'(a,i4.4)')'log.'//trim(complog)//'.f',int(fhour)
    end if
    if (present(lastrestart)) then
       call ESMF_TimeGet(lastrestart, yy=yr, mm=mon, dd=day, h=hour, m=minute, s=sec, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
       write(timestring,'(6i8)')yr,mon,day,hour,minute,sec
    end if

    open(newunit=fh_logunit,file=trim(filename))
    write(fh_logunit,'(a)')'completed: '//trim(complog)
    write(fh_logunit,'(a,f10.3)')'forecast hour:',fhour
    write(fh_logunit,'(a)')'valid time: '//trim(nexttimestring)
    if (len_trim(lastout) > 0) then
       write(fh_logunit,'(a)')'last output: '//trim(lastout)
    end if
    if (present(lastrestart)) then
       write(fh_logunit,'(a)')'last restart: '//trim(timestring)
    end if
    close(fh_logunit)

  end subroutine log_restart_fh

  subroutine finalize_restart_fh(restartfh_info)
    !
    ! !DESCRIPTION:
    ! Clean-up...release allocated memory
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(is_restart_fh_type), intent(inout) :: restartfh_info
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------

    if (allocated(restartfh_info%restartFhTimes)) deallocate(restartfh_info%restartFhTimes)

  end subroutine finalize_restart_fh

end module shr_is_restart_fh_mod
