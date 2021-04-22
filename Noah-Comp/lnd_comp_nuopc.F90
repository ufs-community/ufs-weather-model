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

  use noah_driver            , only : init_driver

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

    use lnd_set_decomp_and_domain , only : lnd_set_decomp_and_domain_from_readmesh ! try out read mask
    
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    character(len=*),parameter :: subname=trim(modName)//':(InitializeRealize) '
    character(len=CL) ::  meshfile_lnd = 'INPUT/C96_181018_ESMFmesh.nc'
    character(len=CL) ::  meshfile_mask
    type(ESMF_Mesh)   ::  mesh_input

    real(r8), pointer   :: lndfrac_glob(:)
    integer,  pointer   :: lndmask_glob(:)
    type(ESMF_DistGrid) :: distgrid
    integer             :: lsize
    integer , pointer   :: lndmask_loc(:)
    integer , pointer   :: itemp_glob(:)
    type(ESMF_Array)    :: elemMaskArray
    integer , pointer   :: gindex(:)
    integer             :: gsize
    integer             :: n
    type(ESMF_VM)       :: vm
    !-------------------------------------------------------------------------------
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)


    ! assume mesh file has mask
    meshfile_mask = meshfile_lnd
    
    ! Read in the land mesh from the file
    
    mesh_input = ESMF_MeshCreate(filename=trim(meshfile_lnd), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine lsize and distgrid_lnd
    call ESMF_MeshGet(mesh_input, elementdistGrid=distgrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distgrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine lndmask_loc
    ! The call to ESMF_MeshGet fills in the values of lndmask_loc
    allocate(lndmask_loc(lsize))
    elemMaskArray = ESMF_ArrayCreate(distgrid, lndmask_loc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_MeshGet(mesh_input, elemMaskArray=elemMaskArray, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine global landmask_glob 
    ! how to get gsize? same as elementCount? No, that's lsize....
    gsize = lsize ! this is a (wrong) palceholder to test

    allocate(gindex(lsize))
    allocate(itemp_glob(gsize))
    ! List of sequence indices for the elements on localDe
    call ESMF_DistGridGet(distgrid, 0, seqIndexList=gindex, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    write(*,*) "JP A2.4. lsize = ", lsize
    
    ! do n = 1,lsize
    !    write(*,*) "JP n, gindex(n), lndmask_loc(n) =", n, gindex(n), lndmask_loc(n)
    !    lndmask_glob(gindex(n)) = lndmask_loc(n)
    ! end do

    ! call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! call ESMF_VMAllReduce(vm, sendData=lndmask_glob, recvData=itemp_glob, count=gsize, &
    !      reduceflag=ESMF_REDUCE_SUM, rc=rc)
    ! write(*,*) "JP A5"
    ! lndmask_glob(:) = int(itemp_glob(:))
    ! deallocate(itemp_glob)
    ! deallocate(gindex)
    ! deallocate(lndmask_loc)
    ! write(*,*) "JP A6"

    ! ! ASSUME that land fraction is identical to land mask in this case                                                                                       
    ! lndfrac_glob(:) = lndmask_glob(:)
    ! write(*,*) "JP A7"

    

    ! ---------------------
    ! Realize the actively coupled fields
    ! ---------------------
    call realize_fields(gcomp, mesh_input, flds_scalar_name, flds_scalar_num, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! --------------------- 
    ! initialize noah:
    ! ---------------------   
    call init_driver


  end subroutine InitializeRealize

  !===============================================================================
  subroutine ModelAdvance(gcomp, rc)

    !use lsm_noah, only: lsm_noah_run
    !use noah_loop, only: noah_loop_run
    use noah_driver, only: noah_loop_drv

    !use field_print_debug, only: field_foo ! test tmp
    !use import_2_export_test, only: import_2_export ! test tmp
    use import_fields, only: import_field
    
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


    !tmp debug
    real(r8), pointer          :: dataPtr_i(:)
    !character(len=*),parameter :: fldname_i="foo_atm2lndfield"
    character(len=*),parameter :: fldname_i="inst_land_sea_mask"
    real(r8), pointer          :: dataPtr_e(:)
    character(len=*),parameter :: fldname_e="foo_lnd2atmfield"
    !real(r8), allocatable      :: lnd_data(:)  ! change this
!    character(len=*),parameter :: fldname_i="Faxa_lwdn"
    character(len=80), allocatable :: flds(:)
    character(len=80)   :: fldname
    integer                   :: n
    ! query the Component for its clock, importState and exportState
    ! call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
    !      exportState=exportState, rc=rc)
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !--------------------------------
    ! Unpack import state
    !--------------------------------

    ! call import_fields(importState,rc)....

    ! test tmp import
    !call field_foo(importState, trim(fldname_i), dataPtr_i, rc=rc)
    !call import_2_export(importState, exportState, fldname_i, fldname_e, rc)

    allocate(flds(7))
    flds = (/'Faxa_lwdn  '    , 'Faxa_swndr '   , 'Faxa_swvdr '   , 'Faxa_swndf ' , 'Faxa_swvdf ', &
         'Faxa_rain  '    , 'Faxa_snow  ' /)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call import_field(importState, fldname, rc)
    end do
    deallocate(flds)

    allocate(flds(22))
    flds=(/ &
         'Faxa_soiltyp     ', &
         'Faxa_vegtype     ', &
         'Faxa_sigmaf      ', &
         'Faxa_sfcemis     ', &
         'Faxa_dlwflx      ', &
         'Faxa_dswsfc      ', &
         'Faxa_snet        ', &
         'Faxa_tg3         ', &
         'Faxa_cm          ', &
         'Faxa_ch          ', &
         'Faxa_prsl1       ', &
         'Faxa_prslki      ', &
         'Faxa_zf          ', &
         'Faxa_land        ', &
         'Faxa_slopetyp    ', &
         'Faxa_shdmin      ', &
         'Faxa_shdmax      ', &
         'Faxa_snoalb      ', &
         'Faxa_sfalb       ', &
         'Faxa_bexppert    ', &
         'Faxa_xlaipert    ', &
         'Faxa_vegfpert    ' &
         /)

    do n = 1,size(flds)
       fldname = trim(flds(n))
       call import_field(importState, fldname, rc)
    end do
    deallocate(flds)

    ! end test tmp

    ! ...end call import_fields
    

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


end module
