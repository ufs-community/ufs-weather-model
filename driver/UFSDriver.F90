!-----------------------------------------------------------------------
!
      MODULE UFSDriver
!
!-----------------------------------------------------------------------
!***  This module contains codes directly related to the UFS Driver component.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  2010-03-24  Black - Created UFS Driver component module.
!  2010-04     Yang  - Added Ensemble capability.
!  2011-05-11  Theurich & Yang - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!  2011-10-04  Yang - Modified for using the ESMF 5.2.0r library.
!  2012-02     Tripp - Added ESMF superstructure to support an OCN model
!  2013-06     Theurich - Reworked OCN dependency to be NUOPC based
!  2013-07     Theurich - Macro based ESMF error handling
!-----------------------------------------------------------------------
!
!***  The UFS Driver component lies in the hierarchy seen here:
!
!          Main program
!               |
!          UFS Driver component
!              /|\
!             / | \
!          ATM/OCN/ICE/WAV/LND/IPM/HYD .. components
!          |    |   |
!          |    |   (CICE, etc.)
!          |    |
!          |    (MOM6, HYCOM, etc.)
!          |
!          (FV3, etc.)
!
!-----------------------------------------------------------------------
!
      use ESMF

      use NUOPC
      use NUOPC_Driver, &
        Driver_routine_SS             => SetServices, &
        Driver_label_SetModelServices => label_SetModelServices, &
        Driver_label_SetRunSequence   => label_SetRunSequence, &
        Driver_label_SetRunClock      => label_SetRunClock
      use NUOPC_Connector, only: conSS => SetServices
      use NUOPC_Model, only: SetVM

  ! - Handle build time ATM options:
#ifdef FRONT_FV3
      use FRONT_FV3,        only: FV3_SS   => SetServices
#endif
#ifdef FRONT_CDEPS_DATM
      use FRONT_CDEPS_DATM, only: DATM_SS  => SetServices
#endif
  ! - Handle build time OCN options:
#ifdef FRONT_HYCOM
      use FRONT_HYCOM,      only: HYCOM_SS  => SetServices
#endif
#ifdef FRONT_MOM6
      use FRONT_MOM6,       only: MOM6_SS   => SetServices, &
                                  MOM6_SV   => SetVM
#endif
#ifdef FRONT_CDEPS_DOCN
      use FRONT_CDEPS_DOCN, only: DOCN_SS  => SetServices
#endif
  ! - Handle build time ICE options:
#ifdef FRONT_CICE6
      use FRONT_CICE6,      only: CICE6_SS => SetServices, &
                                  CICE6_SV => SetVM
#endif
  ! - Handle build time WAV options:
#ifdef FRONT_WW3
      use FRONT_WW3,        only: WW3_SS  => SetServices, &
                                  WW3_SV  => SetVM
#endif
  ! - Handle build time LND options:
#ifdef FRONT_NOAH
      use FRONT_NOAH,       only: NOAH_SS  => SetServices
#endif
#ifdef FRONT_NOAHMP
      use FRONT_NOAHMP,     only: NOAHMP_SS  => SetServices
#endif
#ifdef FRONT_LIS
      use FRONT_LIS,        only: LIS_SS   => SetServices
#endif
  ! - Handle build time IPM options:
#ifdef FRONT_IPE
      use FRONT_IPE,        only: IPE_SS   => SetServices
#endif
  ! - Handle build time AQM options:
#ifdef FRONT_AQM
      use FRONT_AQM,        only: AQM_SS  => SetServices
#endif
  ! - Handle build time GOCART options:
#ifdef FRONT_GOCART
      use FRONT_GOCART,     only: GOCART_SS  => SetServices, &
                                  GOCART_SV  => SetVM
#endif
  ! - Mediator
#ifdef FRONT_CMEPS
      use MED,              only: MED_SS     => SetServices, &
                                  MED_SV     => SetVM
#endif
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      private
!
      public :: UFSDriver_SS
!
!-----------------------------------------------------------------------
!
      logical, private :: flag_verbose_diagnostics = .false.
      logical, private :: printattr = .false.

      character(len=*),parameter :: u_FILE_u = &
           __FILE__

      contains

      logical function ChkErr(rc, line, file)
        integer, intent(in) :: rc            !< return code to check
        integer, intent(in) :: line          !< Integer source line number
        character(len=*), intent(in) :: file !< User-provided source file name
        integer :: lrc
        ChkErr = .false.
        lrc = rc
        if (ESMF_LogFoundError(rcToCheck=lrc, msg=ESMF_LOGERR_PASSTHRU, line=line, file=file)) then
           ChkErr = .true.
        endif
      end function ChkErr

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE UFSDriver_SS(driver,RC)
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      type(ESMF_GridComp) :: driver
!
      integer,intent(out) :: rc
!
!---------------------
!***  Local Variables
!---------------------
!
      type(ESMF_Config)             :: config
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC = ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      ! Derive from NUOPC_Driver
      call NUOPC_CompDerive(driver, Driver_routine_SS, rc=RC)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! specializations:

      call NUOPC_CompSpecialize(driver, &
        specLabel=Driver_label_SetModelServices, specRoutine=SetModelServices, &
        rc=RC)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      call NUOPC_CompSpecialize(driver, &
        specLabel=Driver_label_SetRunSequence, specRoutine=SetRunSequence, &
        rc=RC)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

#ifndef JEDI_DRIVER
      ! The UFS Driver component is currently the top-level driver and
      ! does not need to coordinate Clocks with its parent.
      call ESMF_MethodRemove(driver, Driver_label_SetRunClock, rc=RC)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call NUOPC_CompSpecialize(driver, &
        specLabel=Driver_label_SetRunClock, specRoutine=NUOPC_NoOp, rc=RC)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
#endif

      ! register an internal initialization method
      call NUOPC_CompSetInternalEntryPoint(driver, ESMF_METHOD_INITIALIZE, &
        phaseLabelList=(/"IPDv04p2"/), userRoutine=ModifyCplLists, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! create, open, and set the config
      config = ESMF_ConfigCreate(rc=RC)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_ConfigLoadFile(config, "nems.configure", rc=RC)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return
      call ESMF_GridCompSet(driver, config=config, rc=RC)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! Load the required entries from the fd_nems.yaml file
      call NUOPC_FieldDictionarySetup("fd_nems.yaml", rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

!-----------------------------------------------------------------------
!
      END SUBROUTINE UFSDriver_SS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!

      subroutine SetModelServices(driver, rc)

        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        ! local variables
        integer                         :: localrc, stat, i, j, petCount
        character(ESMF_MAXSTR)          :: name
        type(ESMF_GridComp)             :: comp
        type(ESMF_Config)               :: config
        type(ESMF_Info)                 :: info
        character(len=32), allocatable  :: compLabels(:)
        integer, allocatable            :: petList(:)
        character(len=10)               :: value
        character(len=20)               :: model, prefix
        character(len=160)              :: msg
        integer                         :: petListBounds(2)
        integer                         :: ompNumThreads
        integer                         :: componentCount
        type(NUOPC_FreeFormat)          :: attrFF, fdFF
        logical                         :: found_comp
        logical                         :: isPresent

        rc = ESMF_SUCCESS

        ! query the Component for info
        call ESMF_GridCompGet(driver, name=name, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        ! get petCount and config
        call ESMF_GridCompGet(driver, petCount=petCount, config=config, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        ! read and ingest free format driver attributes
        attrFF = NUOPC_FreeFormatCreate(config, label="EARTH_attributes::", &
          relaxedflag=.true., rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_CompAttributeIngest(driver, attrFF, addFlag=.true., rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        ! dump the current field dictionary into the Log file
        call ESMF_AttributeGet(driver, name="DumpFieldDictionary", &
          value=value, defaultValue="false", &
          convention="NUOPC", purpose="Instance", rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        if (trim(value)=="true") then
          call ESMF_LogWrite( &
            "===>===>===>===> Begin Dumping Field Dictionary <===<===<===<===",&
            ESMF_LOGMSG_INFO)
          call NUOPC_FieldDictionaryEgest(fdFF, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_FreeFormatLog(fdFF, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite( &
            "===>===>===>===> Done Dumping Field Dictionary <===<===<===<===", &
            ESMF_LOGMSG_INFO)
        endif

        ! determine the generic component labels
        componentCount = ESMF_ConfigGetLen(config, &
          label="EARTH_component_list:", rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        allocate(compLabels(componentCount), stat=stat)
        if (ESMF_LogFoundAllocError(statusToCheck=stat, &
          msg="Allocation of compLabels failed.", &
          line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
          return  ! bail out
        call ESMF_ConfigGetAttribute(config, valueList=compLabels, &
          label="EARTH_component_list:", count=componentCount, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        call ReadAttributes(driver, config, "DRIVER_attributes::",  relaxedflag=.true., &
          formatprint=printattr, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        call ReadAttributes(driver, config, "ALLCOMP_attributes::", relaxedflag=.true., &
          formatprint=printattr, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        ! get starttype and set read_restart attribute in driver config list
        call InitRestart(driver, config, rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        ! determine information for each component and add to the driver
        do i=1, componentCount
          ! construct component prefix
          prefix=trim(compLabels(i))
          ! read in petList bounds
          call ESMF_ConfigGetAttribute(config, petListBounds, &
            label=trim(prefix)//"_petlist_bounds:", default=-1, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! handle the default situation
          if (petListBounds(1)==-1 .or. petListBounds(2)==-1) then
            petListBounds(1) = 0
            petListBounds(2) = petCount - 1
          endif
          ! read in model instance name
          call ESMF_ConfigGetAttribute(config, model, &
            label=trim(prefix)//"_model:", default="none", rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! check that there was a model instance specified
          if (trim(model) == "none") then
            ! Error condition: no model was specified
            write (msg, *) "No model was specified for component: ",trim(prefix)
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return  ! bail out
          endif
          ! set petList for this component
          allocate(petList(petListBounds(2)-petListBounds(1)+1))
          do j=petListBounds(1), petListBounds(2)
            petList(j-petListBounds(1)+1) = j ! PETs are 0 based
          enddo

! *** read in number of OpenMP threads for this component
          call ESMF_ConfigGetAttribute(config, ompNumThreads, &
            label=trim(prefix)//"_omp_num_threads:", default=-1, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

! *** create info object
          info = ESMF_InfoCreate(rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

! *** set up PePerPet NUOPC hints
          if (ompNumThreads /= -1) then
            call ESMF_InfoSet(info, key="/NUOPC/Hint/PePerPet/MaxCount", &
              value=ompNumThreads, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif

          found_comp = .false.
#ifdef FRONT_FV3
          if (trim(model) == "fv3") then
            call NUOPC_DriverAddComp(driver, trim(prefix), FV3_SS, &
              info=info, petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#if defined FRONT_CDEPS_DATM
          if (trim(model) == "datm" ) then
            !TODO: Remove bail code and pass info and SetVM to DriverAddComp
            !TODO: once component supports threading.
            if (ompNumThreads > 1) then
              write (msg, *) "ESMF-aware threading NOT implemented for model: "//&
                trim(model)
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
                file=__FILE__, rcToReturn=rc)
              return  ! bail out
            endif
            call NUOPC_DriverAddComp(driver, trim(prefix), DATM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_HYCOM
          if (trim(model) == "hycom") then
            call NUOPC_DriverAddComp(driver, trim(prefix), HYCOM_SS, &
              SetVM, info=info, petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_MOM6
          if (trim(model) == "mom6") then
            call NUOPC_DriverAddComp(driver, trim(prefix), MOM6_SS, &
               MOM6_SV, info=info, petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_CDEPS_DOCN
          if (trim(model) == "docn") then
            !TODO: Remove bail code and pass info and SetVM to DriverAddComp
            !TODO: once component supports threading.
            if (ompNumThreads > 1) then
              write (msg, *) "ESMF-aware threading NOT implemented for model: "//&
                trim(model)
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
                file=__FILE__, rcToReturn=rc)
              return  ! bail out
            endif
            call NUOPC_DriverAddComp(driver, trim(prefix), DOCN_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_CICE6
          if (trim(model) == "cice6") then
            call NUOPC_DriverAddComp(driver, trim(prefix), CICE6_SS, &
              CICE6_SV, info=info, petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_WW3
          if (trim(model) == "ww3") then
            call NUOPC_DriverAddComp(driver, trim(prefix), WW3_SS, &
              WW3_SV, info=info, petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_NOAH
          if (trim(model) == "noah") then
            !TODO: Remove bail code and pass info and SetVM to DriverAddComp
            !TODO: once component supports threading.
            if (ompNumThreads > 1) then
              write (msg, *) "ESMF-aware threading NOT implemented for model: "//&
                trim(model)
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
                file=__FILE__, rcToReturn=rc)
              return  ! bail out
            endif
            call NUOPC_DriverAddComp(driver, trim(prefix), NOAH_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_NOAHMP
          if (trim(model) == "noahmp") then
            call NUOPC_DriverAddComp(driver, trim(prefix), NOAHMP_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_LIS
          if (trim(model) == "lis") then
            !TODO: Remove bail code and pass info and SetVM to DriverAddComp
            !TODO: once component supports threading.
            if (ompNumThreads > 1) then
              write (msg, *) "ESMF-aware threading NOT implemented for model: "//&
                trim(model)
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
                file=__FILE__, rcToReturn=rc)
              return  ! bail out
            endif
            call NUOPC_DriverAddComp(driver, trim(prefix), LIS_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_IPE
          if (trim(model) == "ipe") then
            !TODO: Remove bail code and pass info and SetVM to DriverAddComp
            !TODO: once component supports threading.
            if (ompNumThreads > 1) then
              write (msg, *) "ESMF-aware threading NOT implemented for model: "//&
                trim(model)
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
                file=__FILE__, rcToReturn=rc)
              return  ! bail out
            endif
            call NUOPC_DriverAddComp(driver, trim(prefix), IPE_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_AQM
          if (trim(model) == "aqm") then
            !TODO: Remove bail code and pass info and SetVM to DriverAddComp
            !TODO: once component supports threading.
            if (ompNumThreads > 1) then
              write (msg, *) "ESMF-aware threading NOT implemented for model: "//&
                trim(model)
              call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
                file=__FILE__, rcToReturn=rc)
              return  ! bail out
            endif
            call NUOPC_DriverAddComp(driver, trim(prefix), AQM_SS, &
              petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_GOCART
          if (trim(model) == "gocart") then
            call NUOPC_DriverAddComp(driver, trim(prefix), GOCART_SS, &
              GOCART_SV, info=info, petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
#ifdef FRONT_CMEPS
          if (trim(model) == "cmeps") then
            call NUOPC_DriverAddComp(driver, trim(prefix), MED_SS, &
               MED_SV, info=info, petList=petList, comp=comp, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
            found_comp = .true.
          end if
#endif
          if (.not. found_comp) then
            write(msg,*) 'No component ',trim(model),' found'
            call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msg, line=__LINE__, &
              file=__FILE__, rcToReturn=rc)
            return
          endif

          call AddAttributes(comp, driver, config, trim(prefix), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! clean-up
          deallocate(petList)
          call ESMF_InfoDestroy(info, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
        enddo

        ! clean-up
        deallocate(compLabels)
      end subroutine

  !-----------------------------------------------------------------------------

  subroutine SetRunSequence(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    character(ESMF_MAXSTR)          :: name
    type(ESMF_Config)               :: config
    type(NUOPC_FreeFormat)          :: runSeqFF

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! read free format run sequence from config
    call ESMF_GridCompGet(driver, config=config, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    runSeqFF = NUOPC_FreeFormatCreate(config, label="runSeq::", rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ingest FreeFormat run sequence
    call NUOPC_DriverIngestRunSequence(driver, runSeqFF, &
      autoAddConnectors=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Diagnostic output
    if(flag_verbose_diagnostics) then
       call NUOPC_DriverPrint(driver, orderflag=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

  end subroutine

  !-----------------------------------------------------------------------------

  recursive subroutine ModifyCplLists(driver, importState, exportState, clock, &
    rc)
    type(ESMF_GridComp)  :: driver
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    character(len=160)              :: name, msg
    type(ESMF_CplComp), pointer     :: connectorList(:)
    integer                         :: i, j, cplListSize
    character(len=160), allocatable :: cplList(:)
    character(len=160)              :: value

    rc = ESMF_SUCCESS

    call ESMF_LogWrite("Driver is in ModifyCplLists()", ESMF_LOGMSG_INFO)

    nullify(connectorList)
    call NUOPC_DriverGetComp(driver, compList=connectorList, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    write (msg,*) "Found ", size(connectorList), " Connectors."// &
      " Modifying CplList Attribute...."
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO)

    do i=1, size(connectorList)
      ! query Connector i for its name
      call ESMF_CplCompGet(connectorList(i), name=name, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      ! access CplList for Connector i
      call NUOPC_CompAttributeGet(connectorList(i), name="CplList", &
        itemCount=cplListSize, rc=rc)
      if (ChkErr(rc,__LINE__,u_FILE_u)) return

      if (cplListSize>0) then
        allocate(cplList(cplListSize))
        call NUOPC_CompAttributeGet(connectorList(i), name="CplList", &
          valueList=cplList, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return

        ! go through all of the entries in the cplList and add options
        do j=1, cplListSize
          cplList(j) = trim(cplList(j))//":DumpWeights=true"
          cplList(j) = trim(cplList(j))//":SrcTermProcessing=1:TermOrder=SrcSeq"
          ! add connection options read in from configuration file
          call ESMF_AttributeGet(connectorList(i), name="ConnectionOptions", &
            value=value, defaultValue="", rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          cplList(j) = trim(cplList(j))//trim(value)
        enddo
        ! store the modified cplList in CplList attribute of connector i
        call NUOPC_CompAttributeSet(connectorList(i), &
          name="CplList", valueList=cplList, rc=rc)
        if (ChkErr(rc,__LINE__,u_FILE_u)) return
        deallocate(cplList)
      endif
    enddo

    deallocate(connectorList)
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine ReadAttributes(gcomp, config, label, relaxedflag, formatprint, rc)

    use ESMF  , only : ESMF_GridComp, ESMF_Config, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use NUOPC , only : NUOPC_FreeFormatCreate, NUOPC_CompAttributeIngest
    use NUOPC , only : NUOPC_FreeFormatDestroy, NUOPC_FreeFormat

    ! input/output arguments
    type(ESMF_GridComp) , intent(inout)        :: gcomp
    type(ESMF_Config)   , intent(in)           :: config
    character(len=*)    , intent(in)           :: label
    logical             , intent(in), optional :: relaxedflag
    logical             , intent(in), optional :: formatprint
    integer             , intent(inout)        :: rc

    ! local variables
    type(NUOPC_FreeFormat)  :: attrFF
    character(len=*), parameter :: subname = "(UFSDriver.F90:ReadAttributes)"
    !-------------------------------------------

    rc = ESMF_SUCCESS

    if (present(relaxedflag)) then
       attrFF = NUOPC_FreeFormatCreate(config, label=trim(label), relaxedflag=.true., rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       attrFF = NUOPC_FreeFormatCreate(config, label=trim(label), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    call NUOPC_CompAttributeIngest(gcomp, attrFF, addFlag=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present (formatprint)) then
       call ESMF_LogWrite('ReadAttributes '//trim(label)//' start:', ESMF_LOGMSG_INFO)
       call NUOPC_FreeFormatLog(attrFF, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('ReadAttributes '//trim(label)//' end:', ESMF_LOGMSG_INFO)
    end if

    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine ReadAttributes

  subroutine InitRestart(driver, config, rc)

    !-----------------------------------------------------
    ! Determine if will restart and read pointer file if appropriate
    !-----------------------------------------------------

    use ESMF         , only : ESMF_GridComp, ESMF_VM, ESMF_GridCompGet, ESMF_VMGet, ESMF_SUCCESS
    use ESMF         , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC        , only : NUOPC_CompAttributeGet, NUOPC_CompAttributeSet, NUOPC_CompAttributeAdd

    ! input/output variables
    type(ESMF_GridComp)    , intent(inout) :: driver
    type(ESMF_Config)      , intent(inout) :: config
    integer                , intent(out)   :: rc

    ! local variables
    logical           :: read_restart   ! read the restart file, based on start_type
    character(len=ESMF_MAXSTR) :: cvalue         ! temporary
    character(len=ESMF_MAXSTR) :: attribute      !
    character(len=*) , parameter :: subname = "(UFSDriver.F90:InitRestart)"
    !-------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    !-----------------------------------------------------
    ! Carry out restart if appropriate
    !-----------------------------------------------------

    read_restart = IsRestart(driver, config, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    attribute = 'read_restart'
    ! Add read_restart to driver attributes
    call NUOPC_CompAttributeAdd(driver, attrList=(/trim(attribute)/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    write(cvalue,*) read_restart
    call NUOPC_CompAttributeSet(driver, name='read_restart', value=trim(cvalue), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite('Set attribute read_restart in driver ', ESMF_LOGMSG_INFO)

  end subroutine InitRestart

  function IsRestart(gcomp, config, rc)

    use ESMF         , only : ESMF_GridComp, ESMF_SUCCESS
    use ESMF         , only : ESMF_LogSetError, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_RC_NOT_VALID
    use NUOPC        , only : NUOPC_CompAttributeGet

    ! input/output variables
    logical                                :: IsRestart
    type(ESMF_GridComp)    , intent(inout) :: gcomp
    type(ESMF_Config)      , intent(inout) :: config
    integer                , intent(out)   :: rc

    ! locals
    logical                               :: isPresent, isSet
    character(len=ESMF_MAXSTR)            :: start_type     ! Type of startup
    character(len=ESMF_MAXSTR)            :: msgstr
    character(len=*) , parameter :: start_type_start = "startup"
    character(len=*) , parameter :: start_type_cont  = "continue"
    character(len=*) , parameter :: start_type_brnch = "branch"
    character(len=*) , parameter  :: subname = "(UFSDriver.F90:IsRestart)"
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! First Determine if restart is read
    call NUOPC_CompAttributeGet(gcomp, name='start_type', value=start_type, &
       isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if( isPresent .and. isSet) then
       if ((trim(start_type) /= start_type_start) .and.  &
           (trim(start_type) /= start_type_cont ) .and.  &
           (trim(start_type) /= start_type_brnch)) then
          write (msgstr, *) subname//': start_type invalid = '//trim(start_type)
          call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg=msgstr, line=__LINE__, file=__FILE__, rcToReturn=rc)
          return
       end if

       !TODO: this is hard-wired to CIME start/continue types in terms of gcomp
       IsRestart = .false.
       if (trim(start_type) == trim(start_type_cont) .or. trim(start_type) == trim(start_type_brnch)) then
          IsRestart = .true.
       end if
    else
       IsRestart = .false.
       call ESMF_LogWrite('No start_type attribute found, setting read_restart false ', ESMF_LOGMSG_INFO)
    endif

  end function IsRestart

  subroutine AddAttributes(gcomp, driver, config, compname, rc)

    ! Add specific set of attributes to components from driver attributes

    use ESMF  , only : ESMF_GridComp, ESMF_Config, ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF  , only : ESMF_LogFoundAllocError, ESMF_ConfigGetLen, ESMF_ConfigGetAttribute
    use NUOPC , only : NUOPC_CompAttributeAdd, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet

    ! input/output parameters
    type(ESMF_GridComp) , intent(inout) :: gcomp
    type(ESMF_GridComp) , intent(in)    :: driver
    type(ESMF_Config)   , intent(inout) :: config
    character(len=*)    , intent(in)    :: compname
    integer             , intent(inout) :: rc

    ! local variables
    integer                        :: n
    integer                        :: stat
    integer                        :: inst_index
    logical                        :: is_present
    character(len=ESMF_MAXSTR)     :: cvalue
    character(len=32), allocatable :: compLabels(:)
    character(len=32), allocatable :: attrList(:)
    integer                        :: componentCount
    character(len=*), parameter    :: subname = "(UFSDriver.F90:AddAttributes)"
    logical                        :: lvalue = .false.
    !-------------------------------------------

    rc = ESMF_Success
    call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)

    !------
    ! Add restart flag to gcomp attributes
    !------
    allocate(attrList(1))
    attrList =  (/"read_restart"/)

    call NUOPC_CompAttributeAdd(gcomp, attrList=attrList, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    do n = 1,size(attrList)
       if (trim(attrList(n)) == "read_restart") then
          call NUOPC_CompAttributeGet(driver, name="read_restart", value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          read(cvalue,*) lvalue

          if (.not. lvalue) then
            call NUOPC_CompAttributeGet(driver, name=trim(attrList(n)), value=cvalue, rc=rc)
            if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call NUOPC_CompAttributeSet(gcomp, name=trim(attrList(n)), value=trim(cvalue), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          print *, trim(attrList(n))
          call NUOPC_CompAttributeGet(driver, name=trim(attrList(n)), value=cvalue, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call NUOPC_CompAttributeSet(gcomp, name=trim(attrList(n)), value=trim(cvalue), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call ESMF_LogWrite('Added attribute '//trim(attrList(n))//' to '//trim(compname), ESMF_LOGMSG_INFO)
    enddo
    deallocate(attrList)

    !------
    ! Add component specific attributes
    !------

    call ReadAttributes(gcomp, config, trim(compname)//"_attributes::", relaxedflag=.true., &
       formatprint=printattr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(gcomp, config, trim(compname)//"_modelio::", relaxedflag=.true., &
       formatprint=printattr, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call ReadAttributes(gcomp, config, "ALLCOMP_attributes::", relaxedflag=.true., &
       formatprint=printattr, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine AddAttributes
!
!-----------------------------------------------------------------------
!
      END MODULE UFSDriver
!
!-----------------------------------------------------------------------
