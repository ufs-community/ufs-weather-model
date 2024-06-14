#define ESMF_ERR_ABORT(rc) if (ESMF_LogFoundError(rc, msg="Aborting UFS", line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
!-----------------------------------------------------------------------
!
      PROGRAM UFS
!
!-----------------------------------------------------------------------
!***  Main Program for UFS.
!***  Define ESMF data types and procedures.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2007-       Black   - Modified from Wei-yu's version
!   2007-09     Black   - Create the Clock here.
!   2009-08     Colon   - Unified NEM-NMM & NEMS-GFS
!   2009-06-29  Black   - Modified for addition of NMM nesting;
!                         added new ATM Driver Component.
!   2009-09     Lu      - Add w3tage calls for resource statistics
!   2009-08     W. Yang - Ensemble GEFS Concurrency Code.
!   2010-03     Jovic/Black - Revised to create NEMS gridded component
!                             for new structure.
!   2010-04     Yang    - Add GEFS and GFS for the revised NEMS.
!   2010-11     Yang    - Add the "Generic Core" to NEMS
!   2011-02     Yang    - Updated to use both the ESMF 4.0.0rp2 library,
!                         ESMF 5 series library and the the
!                         ESMF 3.1.0rp2 library.
!   2011-05     Theurich & Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-10     Yang    - Modified for using the ESMF 5.2.0r library.
!   2013-07     Theurich - Macro based ESMF error handling
!   2016-11     Trahan  - Resource usage reporting
!
!-----------------------------------------------------------------------
!
      USE MPI
      USE ESMF
!
!-----------------------------------------------------------------------
!***  USE the UFSDriver module.  Although it
!***  contains the calls to Register and the top level Initialize,
!***  Run, and Finalize, only the Register routine is public.
!-----------------------------------------------------------------------
!
      USE UFSDriver, only : UFSDriver_SS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  Local Variables.
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE                                                   &  !<-- The MPI task ID
                ,NSECONDS_FCST                                          &  !<-- Length of forecast in seconds
                ,YY,MM,DD                                               &  !<-- Time variables for date
                ,HH,MNS,SEC                                                !<-- Time variables for time of day

      REAL(ESMF_KIND_R8) :: fhrot                                          !< forecast hour at restart time
!
      REAL :: NHOURS_FCST                                                  !<-- Length of forecast in hours

      TYPE(ESMF_TimeInterval) :: RUNDURATION                            &  !<-- The ESMF time. The total forecast hours.
                                ,TIMESTEP                               &  !<-- The ESMF timestep length (we only need a dummy here)
                                ,restartOffset
!
      TYPE(ESMF_Time) :: CURRTIME                                       &  !<-- The ESMF current time.
                        ,STARTTIME                                         !<-- The ESMF start time.
!
      TYPE(ESMF_VM) :: VM                                                  !<-- The ESMF virtual machine,
                                                                           !    which contains and manages
                                                                           !    the computer CPU resource
                                                                           !    for the ESMF grid components.
!
      TYPE(ESMF_GridComp) :: UFSDriverComp                                 !<-- The UFS Driver gridded component.
!
      TYPE(ESMF_Clock) :: CLOCK_MAIN                                       !<-- The ESMF time management clock
!
      TYPE(ESMF_Config) :: CF_MAIN                                         !<-- The Configure object
!
      CHARACTER(ESMF_MAXSTR) :: MESSAGE_CHECK
!
      CHARACTER(len=MPI_MAX_LIBRARY_VERSION_STRING) :: library_version
      INTEGER :: resultlen
!
      INTEGER :: RC, RC_USER                                               !<-- The running error signal
!
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Initialize the ESMF framework.
!-----------------------------------------------------------------------
!
      CALL ESMF_Initialize(configFileName="ufs.configure"               & !<-- top level configuration
                          ,defaultCalKind =ESMF_CALKIND_GREGORIAN       & !<-- Set up the default calendar.
                          ,VM             =VM                           & !<-- The ESMF Virtual Machine
                          ,rc             =RC)
      ESMF_ERR_ABORT(RC)

#if 0
!***  Useful when debugging ESMF-managed threading and resource control
      call ESMF_VMLog(vm, rc=rc)
      ESMF_ERR_ABORT(RC)
#endif
!
!-----------------------------------------------------------------------
!***  Extract the MPI task ID.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Obtain the local task ID"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The ESMF Virtual Machine
                     ,localpet=MYPE                                     &  !<-- The local MPI task ID
                     ,rc      =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Print subversion version and other status information.
!-----------------------------------------------------------------------
!
      if (mype == 0) then
        call w3tagb('ufs-weather-model',0,0,0,'np23')
        call MPI_Get_library_version(library_version, resultlen, rc)
        write(*,'(A,A)')      'MPI Library = ', library_version(1:resultlen)
        write(*,'(A,I0,A,I0)')'MPI Version = ', mpi_version,'.',mpi_subversion
      endif
!
!-----------------------------------------------------------------------
!***  Set up the default log.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Up ESMF Log"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_LogSet(flush      =.false.                            &
                      ,trace      =.false.                            &
                      ,rc         =RC)
      ESMF_ERR_ABORT(RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create and load the Configure object which will hold the contents
!***  of the Main configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Main Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CF_MAIN=ESMF_ConfigCreate(rc=RC)
!
      CALL ESMF_ConfigLoadFile(config  =CF_MAIN                         &  !<-- The Configure object
                              ,filename='model_configure'               &  !<-- The name of the configure file
                              ,rc      =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the UFS Driver gridded component which will create and
!***  control the ATM (atmoshpere), OCN (ocean), ICE (sea ice), etc.
!***  gridded components.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the UFS Driver Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      UFSDriverComp=ESMF_GridCompCreate(name   ='UFS Driver Grid Comp'  &
                                         ,rc     = RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the UFS Driver gridded component's Initialize, Run and
!***  Finalize routines.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register UFS Driver Gridded Component Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetServices(UFSDriverComp                       &
                                   ,UFSDriver_SS                        &
                                   ,rc=RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the main ESMF Clock.
!***  The Clock is needed for all calls to Init, Run, and Finalize.
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set the start time in the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Year from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =YY                            &
                                  ,label ='start_year:'                 &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Month from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =MM                            &
                                  ,label ='start_month:'                &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Day from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =DD                            &
                                  ,label ='start_day:'                  &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Hour from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =HH                            &
                                  ,label ='start_hour:'                 &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Minute from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =MNS                           &
                                  ,label ='start_minute:'               &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Second from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =SEC                           &
                                  ,label ='start_second:'               &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Start Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeSet(time=STARTTIME                                  &  !<-- The start time of the forecast (ESMF)
                       ,yy  =YY                                         &  !<-- Year from config file
                       ,mm  =MM                                         &  !<-- Month from config file
                       ,dd  =DD                                         &  !<-- Day from config file
                       ,h   =HH                                         &  !<-- Hour from config file
                       ,m   =MNS                                        &  !<-- Minute from config file
                       ,s   =SEC                                        &  !<-- Second from config file
                       ,rc  =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set the run duration in the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Forecast Length from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF_MAIN                       &
                                  ,value =NHOURS_FCST                   &
                                  ,label ='nhours_fcst:'                &
                                  ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NSECONDS_FCST=nint(NHOURS_FCST*3600.)                             !<-- The forecast length (sec) (REAL)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Length"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION                &  !<-- The forecast length (s) (ESMF)
                               ,s           =NSECONDS_FCST              &
                               ,rc          =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now the Main Clock can be created.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

! NUOPC requires a correct timeStep to function. Here TIMESTEP was just created
! pro forma, but now is replaced with RUNDURATION to become meaningful.

      TIMESTEP = RUNDURATION

      CLOCK_MAIN=ESMF_ClockCreate(name       ='CLOCK_MAIN'              &  !<-- The top-level ESMF Clock
                                 ,timeStep   =TIMESTEP                  &  !<-- Timestep needed by the Clock (ESMF)
                                 ,startTime  =STARTTIME                 &  !<-- The integration start time (ESMF)
                                 ,runDuration=RUNDURATION               &  !<-- The integration duration (ESMF)
                                 ,rc         =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Adjust the currTime of the main clock: CLOCK_MAIN
!***  if the fhrot is > 0
!***  This will correctly set the UFS Driver clocks in case of
!***  Restart-From-History.
!-----------------------------------------------------------------------
      CALL ESMF_ConfigGetAttribute(config   = CF_MAIN          &
                                   ,value   = fhrot            &
                                   ,label   = 'fhrot:'         &
                                   ,default = 0.0_ESMF_KIND_R8 &
                                   ,rc      = RC)
      ESMF_ERR_ABORT(RC)

      if (fhrot > 0) then
        CALL ESMF_TimeIntervalSet(restartOffset, h_r8=fhrot, rc=RC)
        ESMF_ERR_ABORT(RC)
        CURRTIME = STARTTIME + restartOffset
        call ESMF_ClockSet(CLOCK_MAIN, currTime=CURRTIME, &
                                       timeStep=(TIMESTEP-restartOffset), &
                                       rc=RC)
        ESMF_ERR_ABORT(RC)
      endif
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the INITIALIZE step for the UFS Driver component.
!***  The Initialize routine that is called here as well as the
!***  Run and Finalize routines invoked below are those specified
!***  in the Register routine called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the UFS Driver Component Initialize Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompInitialize(gridcomp   =UFSDriverComp            &
                                  ,clock      =CLOCK_MAIN               &
                                  ,userRc     =RC_USER                  &
                                  ,rc         =RC)
      ESMF_ERR_ABORT(RC)
      ESMF_ERR_ABORT(RC_USER)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the RUN step for the UFS Driver component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the UFS Driver Component Run Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompRun(gridcomp   =UFSDriverComp                   &
                           ,clock      =CLOCK_MAIN                      &
                           ,userRc     =RC_USER                         &
                           ,rc         =RC)
      ESMF_ERR_ABORT(RC)
      ESMF_ERR_ABORT(RC_USER)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Execute the FINALIZE step for the UFS Driver component.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the UFS Driver Component Finalize Step"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompFinalize(gridcomp   =UFSDriverComp              &
                                ,clock      =CLOCK_MAIN                 &
                                ,userRc     =RC_USER                    &
                                ,rc         =RC)
      ESMF_ERR_ABORT(RC)
      ESMF_ERR_ABORT(RC_USER)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_MAIN                           &
                            ,rc   =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Destroy the Main Configure object.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Main Configure Object"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigDestroy(config=CF_MAIN                            &
                             ,rc    =RC)
      ESMF_ERR_ABORT(RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK = "Destroy ESMF Grid Comp and Cpl Comp"
!      CALL ESMF_LogWrite(MESSAGE_CHECK, ESMF_LOGMSG_INFO, rc = RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompDestroy(gridcomp=UFSDriverComp                  &
                               ,rc      =RC)
      ESMF_ERR_ABORT(RC)

!-----------------------------------------------------------------------
!***  Shut down the ESMF system.
!-----------------------------------------------------------------------
!
      CALL ESMF_Finalize()
!
!-----------------------------------------------------------------------
!
      if (mype==0) call w3tage('ufs      ')
!
!-----------------------------------------------------------------------
!
      END PROGRAM UFS
!
!-----------------------------------------------------------------------
