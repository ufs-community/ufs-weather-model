! model test after CTSM

module lnd_comp_nuopc

  !-----------------------------------------------------------------------------
  ! LND Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC                  , only : NUOPC_CompDerive, NUOPC_CompSetEntryPoint, NUOPC_CompSpecialize
  use NUOPC                  , only : NUOPC_CompFilterPhaseMap, NUOPC_CompAttributeGet, NUOPC_CompAttributeSet
  use NUOPC_Model            , only : model_routine_SS           => SetServices
  use NUOPC_Model            , only : SetVM
  use NUOPC_Model            , only : model_label_Advance        => label_Advance
  use NUOPC_Model            , only : model_label_DataInitialize => label_DataInitialize
  use NUOPC_Model            , only : model_label_SetRunClock    => label_SetRunClock
  use NUOPC_Model            , only : model_label_Finalize       => label_Finalize
  use NUOPC_Model            , only : NUOPC_ModelGet
  use shr_kind_mod           , only : r8 => shr_kind_r8, cl=>shr_kind_cl
  use nuopc_shr_methods      , only : chkerr, state_setscalar, state_getscalar, state_diagnose, alarmInit
  use nuopc_shr_methods      , only : set_component_logging, get_component_instance, log_clock_advance
  use clm_varctl             , only : inst_index, inst_suffix, inst_name
  use lnd_import_export      , only : advertise_fields, realize_fields
  ! use lnd_import_export      , only : advertise_fields, realize_fields, import_fields, export_fields
  !use lnd_comp_shr           , only : mesh, model_meshfile !, model_clock
  use shr_sys_mod            , only : shr_sys_abort
  
  implicit none
  private ! except 

  ! Module public routines 
  public  :: SetServices
  public  :: SetVM

  ! Module private routines  
  private :: InitializeP0
  private :: InitializeAdvertise
  private :: InitializeRealize
  private :: ModelAdvance
  private :: ModelFinalize
  
  character(len=CL)      :: flds_scalar_name = ''
  integer                :: flds_scalar_num = 0

  character(*),parameter :: modName =  "(lnd_comp_nuopc)"
  character(len=*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================   
contains
!===============================================================================   

  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    character(len=*),parameter  :: subname=trim(modName)//':(SetServices) '

    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! the NUOPC gcomp component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! switching to IPD versions    
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         userRoutine=InitializeP0, phase=0, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! set entry point for methods that require specific implementation 
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

!!!!!!!!!!!! JP: just try to get this far for now
    ! attach specializing method(s)            
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
         specRoutine=ModelAdvance, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
         specRoutine=ModelFinalize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    
  end subroutine SetServices


  !=============================================================================== 
  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)

    ! input/output variables 
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    !-------------------------------------------------------------------------------           

    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries   
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, acceptStringList=(/"IPDv01p"/), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine InitializeP0

  !===============================================================================  

  !===============================================================================  
  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_VM)      :: vm
    integer            :: lmpicom
    integer            :: ierr
    integer            :: n
    integer            :: localpet
    integer            :: compid      ! component id    
    integer            :: shrlogunit  ! original log unit             
    character(len=CL)  :: cvalue
    character(len=CL)  :: logmsg
    logical            :: isPresent, isSet
    logical            :: cism_evolve
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------
    rc = ESMF_SUCCESS
    
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call get_component_instance(gcomp, inst_suffix, inst_index, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    inst_name = 'LND'

    !----------------------------------------------------------------------------      
    ! advertise fields 
    !---------------------------------------------------------------------------- 

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldName", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       flds_scalar_name = trim(cvalue)
       call ESMF_LogWrite(trim(subname)//' flds_scalar_name = '//trim(flds_scalar_name), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else
    !    call shr_sys_abort(subname//'Need to set attribute ScalarFieldName')
    endif

    call NUOPC_CompAttributeGet(gcomp, name="ScalarFieldCount", value=cvalue, isPresent=isPresent, isSet=isSet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (isPresent .and. isSet) then
       read(cvalue, *) flds_scalar_num
       write(logmsg,*) flds_scalar_num
       call ESMF_LogWrite(trim(subname)//' flds_scalar_num = '//trim(logmsg), ESMF_LOGMSG_INFO)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else
    !    call shr_sys_abort(subname//'Need to set attribute ScalarFieldCount')
    endif

    call advertise_fields(gcomp, flds_scalar_name, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    
  end subroutine InitializeAdvertise


  !===============================================================================  
  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    ! input/output variables      
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    character(len=CL) ::  meshfile_lnd = 'INPUT/C96_181018_ESMFmesh.nc'
    type(ESMF_Mesh)   ::  mesh
    
    !-------------------------------------------------------------------------------
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! Read in the land mesh from the file
    mesh = ESMF_MeshCreate(filename=trim(meshfile_lnd), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ! Read in mask meshfile if needed       
    ! if (trim(meshfile_mask) /= trim(meshfile_lnd)) then
    !    mesh_maskinput = ESMF_MeshCreate(filename=trim(meshfile_mask), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    !    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! end if

    ! ---------------------  
    ! Realize the actively coupled fields        
    ! ---------------------            
    call realize_fields(gcomp, mesh, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    
  end subroutine InitializeRealize

  !===============================================================================
  subroutine ModelAdvance(gcomp, rc)

    use lsm_noah, only: lsm_noah_run ! tmp test
    
    ! Arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! Local variables
    type(ESMF_Clock)           :: clock,dclock,mclock
    type(ESMF_Alarm)           :: alarm
    type(ESMF_Time)            :: startTime
    type(ESMF_Time)            :: currTime
    type(ESMF_Time)            :: nextTime
    type(ESMF_TimeInterval)    :: timeStep
    type(ESMF_State)           :: importState, exportState
    character(len=*),parameter :: subname=trim(modName)//':(ModelAdvance) '
    character(len=256)   :: msgString
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
         exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_ClockPrint(clock, options="currTime", &
         preString="------>Advancing LND from: ", unit=msgString, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(subname//trim(msgString), ESMF_LOGMSG_INFO)

    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, &
         timeStep=timeStep, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_TimePrint(currTime + timeStep, &
         preString="--------------------------------> to: ", unit=msgString, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO)

    ! call ESMF_ClockAdvance(clock,rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return


  end subroutine ModelAdvance

  !===============================================================================

  subroutine ModelFinalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname=trim(modName)//':(ModelFinalize) '

    ! begin
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)
  end subroutine ModelFinalize

    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!   use ESMF
!   use NUOPC
!   use NUOPC_Model, only: &
!        model_routine_SS      => SetServices, &
!        model_label_SetClock  => label_SetClock, &
!        model_label_Advance   => label_Advance
!   !use NUOPC_Model        , only : NUOPC_ModelGet
!   use NUOPC_Model 
!   implicit none

!   ! -- import fields
!   integer, parameter :: importFieldCount = 1
!   character(len=*), dimension(importFieldCount), parameter :: &
!        importFieldNames = (/  "mean_down_lw_flx"  /)
!   ! -- export fields
!   integer, parameter :: exportFieldCount = 1
!   character(len=*), dimension(exportFieldCount), parameter :: &
!        exportFieldNames = (/  "inst_tracer_mass_frac"  /)

!   ! -- verbosity
!   integer :: verbosity


!   private

!   public SetServices

!   !-----------------------------------------------------------------------------
! contains
!   !-----------------------------------------------------------------------------

!   subroutine SetServices(model, rc)
!     type(ESMF_GridComp)  :: model
!     integer, intent(out) :: rc

!     ! begin
!     rc = ESMF_SUCCESS

!     ! the NUOPC model component will register the generic methods
!     call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!     ! Provide InitializeP0 to switch from default IPDv00 to IPDv03
!     call ESMF_GridCompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
!          userRoutine=InitializeP0, phase=0, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!     ! set entry point for methods that require specific implementation
!     call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
!          phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeP1, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!     ! attach specializing method(s)
!     call NUOPC_CompSpecialize(model, specLabel=label_DataInitialize, &
!          specRoutine=DataInitialize, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
!     call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
!          specRoutine=ModelAdvance, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
!     call NUOPC_CompSpecialize(model, specLabel=label_Finalize, &
!          specRoutine=ModelFinalize, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!   end subroutine SetServices

!   !-----------------------------------------------------------------------------

!   subroutine InitializeP0(model, importState, exportState, clock, rc)
!     type(ESMF_GridComp)  :: model
!     type(ESMF_State)     :: importState, exportState
!     type(ESMF_Clock)     :: clock
!     integer, intent(out) :: rc

!     ! local variables
!     character(len=5)     :: value

!     ! begin
!     rc = ESMF_SUCCESS

!     ! get component verbosity
!     call ESMF_AttributeGet(model, name="Verbosity", value=value, &
!          defaultValue="min", convention="NUOPC", purpose="Instance", rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
!     ! convert value to verbosity
!     verbosity = ESMF_UtilString2Int(value, &
!          specialStringList=(/"min ","high", "max "/), &
!          specialValueList=(/0,255,255/), rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
!     ! write to log files
!     call ESMF_LogWrite("CHM: Verbosity = " // trim(value), &
!          ESMF_LOGMSG_INFO, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!     ! switch to IPDv03 by filtering all other phaseMap entries
!     call NUOPC_CompFilterPhaseMap(model, ESMF_METHOD_INITIALIZE, &
!          acceptStringList=(/"IPDv03p"/), rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!   end subroutine InitializeP0


!   !-----------------------------------------------------------------------------
!   subroutine InitializeP1(model, importState, exportState, clock, rc)
!     type(ESMF_GridComp)  :: model
!     type(ESMF_State)     :: importState, exportState
!     type(ESMF_Clock)     :: clock
!     integer, intent(out) :: rc

!     ! begin
!     rc = ESMF_SUCCESS

!     call ESMF_LogWrite("JP, begin InitializeP1",ESMF_LOGMSG_INFO, rc=rc)
    
!     ! -- advertise imported fields
!     if (importFieldCount > 0) then
!        call NUOPC_Advertise(importState, importFieldNames, &
!             TransferOfferGeomObject="cannot provide", &
!             SharePolicyField="share", &
!             rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out
!     end if

!     ! -- advertise exported fields
!     if (exportFieldCount > 0) then
!        call NUOPC_Advertise(exportState, exportFieldNames, &
!             TransferOfferGeomObject="cannot provide", &
!             SharePolicyField="share", &
!             rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out
!     end if

!     call ESMF_LogWrite("JP, end InitializeP1",ESMF_LOGMSG_INFO, rc=rc)
    
!   end subroutine InitializeP1

!   !-----------------------------------------------------------------------------


!   !-----------------------------------------------------------------------------

!   subroutine DataInitialize(model, rc)
!     type(ESMF_GridComp)  :: model
!     integer, intent(out) :: rc

!     ! local variables
!     type(ESMF_State)              :: importState, exportState
!     type(ESMF_Field)              :: field
!     type(ESMF_Clock)              :: clock
!     type(ESMF_Alarm)              :: alarm
!     type(ESMF_Grid)               :: grid
!     type(ESMF_VM)                 :: vm
!     type(ESMF_GeomType_flag)      :: geomtype
!     type(ESMF_DistGrid)           :: distgrid
!     type(ESMF_Array)              :: array
!     type(ESMF_Config)             :: cf
!     integer                       :: de, item, localrc, localDe, tile
!     integer                       :: comm
!     real(ESMF_KIND_R8), dimension(:,:), pointer :: coord

!     integer :: dimCount, tileCount, deCount, localDeCount
!     integer, dimension(:),   allocatable :: deToTileMap, localDeToDeMap
!     integer, dimension(:,:), allocatable :: minIndexPDe, maxIndexPDe, minIndexPTile, maxIndexPTile
!     integer, dimension(:,:), allocatable :: computationalLBound, computationalUBound

!     logical                    :: restart
!     integer                    :: restart_interval
!     integer                    :: yy, mm, dd, h, m
!     integer(ESMF_KIND_I8)      :: advanceCount
!     real(ESMF_KIND_R8)         :: dts
!     type(ESMF_Time)            :: startTime, currTime
!     type(ESMF_TimeInterval)    :: TimeStep
!     character(len=255)         :: msgString

!     character(len=*), parameter :: config_file = "model_configure"

!     ! begin
!     call ESMF_LogWrite("JP, begin DataInitialize",ESMF_LOGMSG_INFO, rc=rc)
!     rc = ESMF_SUCCESS

!     ! -- initialize chemistry model
!     call ESMF_GridCompGet(model, vm=vm, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!     ! -- check if import fields are defined
!     if (importFieldCount < 1) then 
!        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
!             msg="This component requires import fields to be defined.", &
!             line=__LINE__, file=__FILE__, &
!             rcToReturn=rc)
!        return  ! bail out
!     end if

!     ! -- check if export fields are defined
!     if (exportFieldCount < 1) then 
!        call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
!             msg="This component requires export fields to be defined.", &
!             line=__LINE__, file=__FILE__, &
!             rcToReturn=rc)
!        return  ! bail out
!     end if

!     ! -- query the Component for its clock, importState and exportState
!     call NUOPC_ModelGet(model, importState=importState, &
!          exportState=exportState, modelClock=clock, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!     ! get coordinates from Grid object
!     ! assume all fields on same grid
!     ! use first field 
!     call ESMF_StateGet(importState, field=field, &
!          itemName=trim(importFieldNames(1)), rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
!     call ESMF_FieldGet(field, geomtype=geomtype, localDeCount=localDeCount, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out

!     if (geomtype == ESMF_GEOMTYPE_GRID) then
!        call ESMF_FieldGet(field, grid=grid, array=array, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out
!        call ESMF_ArrayGet(array, deCount=deCount, dimCount=dimCount, &
!             tileCount=tileCount, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out
!        allocate(minIndexPDe(dimCount, deCount), maxIndexPDe(dimCount, deCount),  &
!             minIndexPTile(dimCount, tileCount), maxIndexPTile(dimCount, tileCount), &
!             computationalLBound(dimCount, localDeCount), computationalUBound(dimCount, localDeCount), &
!             deToTileMap(deCount), localDeToDeMap(localDeCount), stat=localrc)
!        if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__,  &
!             file=__FILE__,  &
!             rcToReturn=rc)) &
!             return  ! bail out
!        call ESMF_ArrayGet(array, distgrid=distgrid, &
!             deToTileMap=deToTileMap, localDeToDeMap=localDeToDeMap, &
!             computationalLBound=computationalLBound, &
!             computationalUBound=computationalUBound, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out
!        call ESMF_DistGridGet(distgrid, &
!             minIndexPDe=minIndexPDe, maxIndexPDe=maxIndexPDe, &
!             minIndexPTile=minIndexPTile, maxIndexPTile=maxIndexPTile, &
!             rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out
!     end if

!     call ESMF_LogWrite("JP, end DataInit",ESMF_LOGMSG_INFO, rc=rc)

!   end subroutine DataInitialize

!   !-----------------------------------------------------------------------------

!      subroutine ModelAdvance(model, rc)
!        type(ESMF_GridComp)  :: model
!        integer, intent(out) :: rc

!        ! local variables
!        type(ESMF_Clock)              :: clock
!        type(ESMF_State)              :: importState, exportState
!        type(ESMF_Time)               :: currTime
!        type(ESMF_TimeInterval)       :: timeStep
!        type(ESMF_Field)              :: field
!        type(ESMF_VM)                 :: vm
!        type(ESMF_Alarm), pointer     :: alarmList(:)
!        integer                       :: alarmCount, item, localrc
!        integer                       :: yy, mm, dd, h, m, s
!        character(len=ESMF_MAXSTR)    :: alarmName
!        character(len=15)             :: timeStamp

!        ! begin
!        rc = ESMF_SUCCESS

!        ! query the Component for its clock, importState and exportState
!        call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
!             exportState=exportState, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out

!        ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

!        ! Because of the way that the internal Clock was set in SetClock(),
!        ! its timeStep is likely smaller than the parent timeStep. As a consequence
!        ! the time interval covered by a single parent timeStep will result in 
!        ! multiple calls to the ModelAdvance() routine. Every time the currTime
!        ! will come in by one internal timeStep advanced. This goes until the
!        ! stopTime of the internal Clock has been reached.

!        call ESMF_ClockPrint(clock, options="currTime", &
!             preString="------>Advancing CHM from: ", rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out

!        call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out

!        call ESMF_TimePrint(currTime + timeStep, &
!             preString="---------------------> to: ", rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!             line=__LINE__, &
!             file=__FILE__)) &
!             return  ! bail out

!        ! print field diagnostics
!        if (btest(verbosity,0)) then
!           call ESMF_GridCompGet(model, vm=vm, rc=rc)
!           if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!                line=__LINE__, &
!                file=__FILE__)) &
!                return  ! bail out
!           do item = 1, importFieldCount
!              call ESMF_StateGet(importState, field=field, &
!                   itemName=trim(importFieldNames(item)), rc=rc)
!              if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!                   line=__LINE__, &
!                   file=__FILE__)) &
!                   return  ! bail out
!              ! call fieldPrintMinMax(field, vm=vm, rc=rc)
!              ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!              !      line=__LINE__, &
!              !      file=__FILE__)) &
!              !      return  ! bail out
!           end do
!        end if

!      end subroutine ModelAdvance

!      !-----------------------------------------------------------------------------

!      subroutine ModelFinalize(model, rc)
!        type(ESMF_GridComp)   :: model
!        integer, intent(out)  :: rc

!        ! begin
!        rc = ESMF_SUCCESS

!        ! finalize model

!      end subroutine ModelFinalize

!   !-----------------------------------------------------------------------------


  
    
! !   subroutine SetServices(gcomp, rc)
! !     type(ESMF_GridComp)  :: gcomp
! !     integer, intent(out) :: rc
    
! !     rc = ESMF_SUCCESS
    
! !     ! the NUOPC model component will register the generic methods
! !     write(*,*) "JP debug 1"
! !     call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
    
! !     ! set entry point for methods that require specific implementation
! !     write(*,*) "JP debug 2"
! !     call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
! !       phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     write(*,*) "JP debug 3"
! !     call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
! !       phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
    
! !     ! attach specializing method(s)
! !     write(*,*) "JP debug 4"
! !     call NUOPC_CompSpecialize(gcomp, specLabel=model_label_SetClock, &
! !       specRoutine=SetClock, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     write(*,*) "JP debug 5"    
! !     call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Advance, &
! !       specRoutine=ModelAdvance, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     write(*,*) "JP debug 5a"    
! !   end subroutine
  
! !   !-----------------------------------------------------------------------------

! !   subroutine InitializeP1(gcomp, importState, exportState, clock, rc)
! !     type(ESMF_GridComp)  :: gcomp
! !     type(ESMF_State)     :: importState, exportState
! !     type(ESMF_Clock)     :: clock
! !     integer, intent(out) :: rc

! !     write(*,*) "JP debug 6"
! !     rc = ESMF_SUCCESS

! ! !!!! add some fields
! !     ! query for importState and exportState
! !     ! call NUOPC_ModelGet(gcomp, importState=importState, &
! !     !   exportState=exportState, rc=rc)
! !     ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !     !   line=__LINE__, &
! !     !   file=__FILE__)) &
! !     !   return  ! bail out

! !     ! importable field: air_pressure_at_sea_level
! !     ! call NUOPC_Advertise(importState, &
! !     !   StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
! !     ! call NUOPC_Advertise(importState, &
! !     !   StandardName="soil_type", name="soil_type", rc=rc)

! !     ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !     !   line=__LINE__, &
! !     !   file=__FILE__)) &
! !     !   return  ! bail out

! !     ! ! importable field: surface_net_downward_shortwave_flux
! !     ! call NUOPC_Advertise(importState, &
! !     !   StandardName="mean_down_lw_flx", name="mean_down_lw_flx", rc=rc)
! !     call NUOPC_Advertise(importState, &
! !          StandardName="mean_down_lw_flx", name="mean_down_lw_flx", &
! !          TransferOfferGeomObject="cannot provide", rc=rc)

! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out

! !     ! exportable field: sea_surface_temperature
! !     ! call NUOPC_Advertise(exportState, &
! !     !   StandardName="inst_down_lw_flx", name="inst_down_lw_flx", rc=rc)
! !     ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !     !   line=__LINE__, &
! !     !   file=__FILE__)) &
! !     !   return  ! bail out
    
! !   end subroutine
  
! !   !-----------------------------------------------------------------------------

! !   subroutine InitializeP2(gcomp, importState, exportState, clock, rc)
! !     type(ESMF_GridComp)  :: gcomp
! !     type(ESMF_State)     :: importState, exportState
! !     type(ESMF_Clock)     :: clock
! !     integer, intent(out) :: rc

! !     type(ESMF_Field)        :: field
! !     type(ESMF_Grid)         :: gridIn
! !     type(ESMF_Grid)         :: gridOut

! !     integer :: dimCount ! tmp, debug
! !     integer :: coordDimCount(ESMF_MAXDIM) ! tmp, debug
! !     write(*,*) "JP debug 7"
    
! !     rc = ESMF_SUCCESS
! ! !!! add some fields and grid
    
! !     ! ! query for importState and exportState
! !     ! call NUOPC_ModelGet(gcomp, importState=importState, &
! !     !   exportState=exportState, rc=rc)
! !     ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !     !   line=__LINE__, &
! !     !   file=__FILE__)) &
! !     !   return  ! bail out

! !     ! create a Grid object for Fields
! !     ! gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/360, 180/), &
! !     !      minCornerCoord=(/0.0_ESMF_KIND_R8, -90.0_ESMF_KIND_R8/), &
! !     !      maxCornerCoord=(/360.0_ESMF_KIND_R8, 90.0_ESMF_KIND_R8/), &
! !     !      coordSys=ESMF_COORDSYS_CART,  &
! !     !      staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
! !     !      rc=rc)
! !     ! gridIn = NUOPC_GridCreateSimpleSph(0._ESMF_KIND_R8, -85._ESMF_KIND_R8, &             
! !     !      360._ESMF_KIND_R8, 85._ESMF_KIND_R8, 400, 200, &     
! !     !      scheme=ESMF_REGRID_SCHEME_FULL3D, rc=rc)                   

! !     gridIn = ESMF_GridCreate1PeriDimUfrm(maxIndex=(/360,180/), &
! !          minCornerCoord=(/0.0_ESMF_KIND_R8,-90.0_ESMF_KIND_R8/), &
! !          maxCornerCoord=(/360.0_ESMF_KIND_R8,90.0_ESMF_KIND_R8/), &
! !          staggerLocList=(/ESMF_STAGGERLOC_CORNER, ESMF_STAGGERLOC_CENTER/), rc=rc)

! !    call ESMF_GridGet(gridIn, dimCount=dimCount, coordDimCount=coordDimCount, &
! !         rc=rc)
! !    !! write out dimcount, coorddimcount
! !    write(*,*) "JP grid dimCount", dimCount
! !    write(*,*) "JP grid coordDimCount", coordDimCount

   
! !    !  !! from tech note
! !    !  gridIn = ESMF_GridCreate1PeriDim(          &
! !    !       ! Define a regular distribution
! !    !       maxIndex=(/360,180/), & ! define index space
! !    !       !regDecomp=(/2,3/),  & ! define how to divide among DEs
! !    !       ! Specify mapping of coords dim to Grid dim
! !    !       coordDep1=(/1/), & ! 1st coord is 1D and depends on 1st Grid dim
! !    !       coordDep2=(/2/), & ! 2nd coord is 1D and depends on 2nd Grid dim
! !    !       indexflag=ESMF_INDEX_GLOBAL, &rc=rc)
! !    ! !-------------------------------------------------------------------
! !    ! ! Allocate coordinate storage and associate it with the center
! !    ! ! stagger location.  Since no coordinate values are specified in
! !    ! ! this call no coordinate values are set yet.
! !    ! !-------------------------------------------------------------------
! !    ! call ESMF_GridAddCoord(gridIn,  &
! !    !        staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
! !    ! !-------------------------------------------------------------------
! !    ! ! Get the pointer to the first coordinate array and the bounds
! !    ! ! of its global indices on the local DE.
! !    ! !-------------------------------------------------------------------
! !    ! call ESMF_GridGetCoord(gridIn, coordDim=1, localDE=0, &
! !    !        staggerloc=ESMF_STAGGERLOC_CENTER, &
! !    !        computationalLBound=lbnd, computationalUBound=ubnd, &
! !    !        farrayPtr=coordX, rc=rc)
! !    ! !-------------------------------------------------------------------
! !    ! ! Calculate and set coordinates in the first dimension [10-100].
! !    ! !-------------------------------------------------------------------
! !    ! do i=lbnd(1),ubnd(1)
! !    !      coordX(i) = i*10.0
! !    ! enddo
! !    ! !-------------------------------------------------------------------
! !    ! ! Get the pointer to the second coordinate array and the bounds of
! !    ! ! its global indices on the local DE.
! !    ! !-------------------------------------------------------------------
! !    ! call ESMF_GridGetCoord(gridIn, coordDim=2, localDE=0, &
! !    !        staggerloc=ESMF_STAGGERLOC_CENTER, &
! !    !        computationalLBound=lbnd, computationalUBound=ubnd, &
! !    !        farrayPtr=coordY, rc=rc)
! !    ! !-------------------------------------------------------------------
! !    ! ! Calculate and set coordinates in the second dimension [10-200]
! !    ! !-------------------------------------------------------------------
! !    ! do j=lbnd(1),ubnd(1)
! !    !      coordY(j) = j*10.0
! !    ! enddo

    
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     gridOut = gridIn ! for now out same as in

! !     ! ! importable field: air_pressure_at_sea_level
! !     ! field = ESMF_FieldCreate(name="soil_type", grid=gridIn, &
! !     !   typekind=ESMF_TYPEKIND_R8, rc=rc)
! !     ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !     !   line=__LINE__, &
! !     !   file=__FILE__)) &
! !     !   return  ! bail out
! !     ! call NUOPC_Realize(importState, field=field, rc=rc)
! !     ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !     !   line=__LINE__, &
! !     !   file=__FILE__)) &
! !     !   return  ! bail out

! !     ! importable field: surface_down_downward_shortwave_flux
! !     field = ESMF_FieldCreate(name="mean_down_lw_flx", grid=gridIn, &
! !       typekind=ESMF_TYPEKIND_R8, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     call NUOPC_Realize(importState, field=field, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out

! !     ! ! exportable field: sea_surface_temperature
! !     ! field = ESMF_FieldCreate(name="land_mask", grid=gridOut, &
! !     !   typekind=ESMF_TYPEKIND_R8, rc=rc)
! !     ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !     !   line=__LINE__, &
! !     !   file=__FILE__)) &
! !     !   return  ! bail out
! !     ! call NUOPC_Realize(exportState, field=field, rc=rc)
! !     ! if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !     !   line=__LINE__, &
! !     !   file=__FILE__)) &
! !     !   return  ! bail out
    
! !   end subroutine
  
! !   !-----------------------------------------------------------------------------

! !   subroutine SetClock(gcomp, rc)
! !     type(ESMF_GridComp)  :: gcomp
! !     integer, intent(out) :: rc
    
! !     ! local variables
! !     type(ESMF_Clock)              :: clock
! !     type(ESMF_TimeInterval)       :: stabilityTimeStep

! !     rc = ESMF_SUCCESS

! !     write(*,*) "JP debug 8"    
! !     ! query the Component for its clock, importState and exportState
! !     call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
      
! !     ! initialize internal clock
! !     ! here: parent Clock and stability timeStep determine actual model timeStep
! !     !TODO: stabilityTimeStep should be read in from configuation
! !     !TODO: or computed from internal Grid information
! !     write(*,*) "JP debug 9"    
! !     call ESMF_TimeIntervalSet(stabilityTimeStep, m=5, rc=rc) ! 5 minute steps
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     write(*,*) "JP debug 10"        
! !     call NUOPC_CompSetClock(gcomp, clock, stabilityTimeStep, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     write(*,*) "JP debug 10a"            
! !   end subroutine

! !   !-----------------------------------------------------------------------------

! !   subroutine ModelAdvance(gcomp, rc)
! !     type(ESMF_GridComp)  :: gcomp
! !     integer, intent(out) :: rc
    
! !     ! local variables
! !     type(ESMF_Clock)            :: clock
! !     type(ESMF_State)            :: importState, exportState
! !     type(ESMF_Time)             :: currTime
! !     type(ESMF_TimeInterval)     :: timeStep
! !     character(len=160)          :: msgString

! !     rc = ESMF_SUCCESS
    
! !     ! query the Component for its clock, importState and exportState
! !     write(*,*) "JP debug 11"        
! !     call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
! !       exportState=exportState, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out

! !     ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
! !     ! Because of the way that the internal Clock was set in SetClock(),
! !     ! its timeStep is likely smaller than the parent timeStep. As a consequence
! !     ! the time interval covered by a single parent timeStep will result in 
! !     ! multiple calls to the ModelAdvance() routine. Every time the currTime
! !     ! will come in by one internal timeStep advanced. This goes until the
! !     ! stopTime of the internal Clock has been reached.
! !     write(*,*) "JP debug 12"    
! !     call ESMF_ClockPrint(clock, options="currTime", &
! !       preString="------>Advancing LND from: ", unit=msgString, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
    
! !     call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
    
! !     call ESMF_TimePrint(currTime + timeStep, &
! !       preString="---------------------> to: ", unit=msgString, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out
! !     call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
! !     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
! !       line=__LINE__, &
! !       file=__FILE__)) &
! !       return  ! bail out

! !   end subroutine

end module
