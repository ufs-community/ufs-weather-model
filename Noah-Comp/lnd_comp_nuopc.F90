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
  use import_fields_mod      , only : ie_set_geomtype
  ! use lnd_import_export      , only : advertise_fields, realize_fields, import_fields, export_fields
  !use lnd_comp_shr           , only : mesh, model_meshfile !, model_clock
  use shr_sys_mod            , only : shr_sys_abort

  use mpp_mod,            only: mpp_init, mpp_pe, mpp_root_pe
  use fms_mod,            only: fms_init
  use fms_io_mod,         only: read_data
  
  use noah_driver,        only: init_driver
  !use land_domain_mod,    only: domain_create 

  
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

  type(ESMF_GeomType_Flag) :: geomtype
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
    integer :: mype, ntasks, mpi_comm_land
    character(len=*), parameter :: subname=trim(modName)//':(InitializeAdvertise) '
    character(len=*), parameter :: format = "('("//trim(subname)//") :',A)"
    !-------------------------------------------------------------------------------
    rc = ESMF_SUCCESS

    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! start putting in fms/mpp stuff here
    call ESMF_VMGetCurrent(vm=VM,rc=RC)
    call ESMF_VMGet(vm=VM, localPet=mype, mpiCommunicator=mpi_comm_land, &
         petCount=ntasks, rc=rc)
    if (mype == 0) write(0,*) 'in lnd comp initadvert, ntasks=',ntasks
    
    ! fms_init will return without doing anything if already called on same proc layout
    call fms_init(mpi_comm_land)


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

    use proc_bounds, only : procbounds, control_init_type
    use fms_io_mod,         only: field_exist, read_data
    
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

    ! tmp, debugging
    integer :: dimCount     ! Number of dimensions (rank) ofdistgrid
    integer :: deCount      ! Number of DEs in the DELayout indistgrid
    integer :: localDeCount ! Number of local DEs in the DELayout in distgrid on this PET
    
    !integer :: localDe      ! Local DE for which information is requested.[0,..,localDeCount-1]
    integer :: de            ! The global DE associated with thelocalDe. DE indexing starts at 0.   
    integer :: elementCount  ! Number of elements in the localDe, i.e. identical to elementCountPDe(localDeToDeMap(localDe))
    integer , pointer   :: seqIndexList(:)  ! List of sequence indices for the elements onlocalDe

    integer :: i

    ! cube sphere mosaic
    type(ESMF_Decomp_Flag)  :: decompflagPTile(2,6)
    character(256)          :: gridfile
    integer                 :: tl
    integer,dimension(2,6)  :: decomptile                  !define delayout for the 6 cubed-sphere tiles
    type(ESMF_Grid)         :: lndGrid
    character(50)           :: gridchoice
    type (control_init_type)::   ctrl_init
    !-------------------------------------------------------------------------------
    rc = ESMF_SUCCESS
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO)

    ! moved to init_driver
    ! ! domain create with FMS:
    ! call domain_create()

    ! --------------------- 
    ! initialize noah:
    ! ---------------------   
    call init_driver(ctrl_init)

    if ( trim(ctrl_init%grid) == '1Dmesh' )  then ! mesh creation
       geomtype = ESMF_GEOMTYPE_MESH
       
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
       gsize = lsize ! this is a (wrong) placeholder to test

       allocate(gindex(lsize))
       allocate(itemp_glob(gsize))
       ! List of sequence indices for the elements on localDe
       !call ESMF_DistGridGet(distgrid, 0, seqIndexList=gindex, rc=rc)
       !if (chkerr(rc,__LINE__,u_FILE_u)) return
       !write(*,*) "JP A2.4. gindex = ", gindex

!!! debug prints
       call ESMF_DistGridGet(distgrid,dimCount=dimCount, deCount=deCount, localDeCount=localDeCount, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       !write(*,*) "DGG1: ", dimCount, deCount, localDeCount

       allocate(seqIndexList(lsize))
       call ESMF_DistGridGet(distgrid, localDe=0, de=de,seqIndexList=seqIndexList, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! do i = 1, lsize
       !    write(*,*) "DGG2: ", de, seqIndexList(i)
       ! end do

       procbounds%de = de
       ! think this will suffice for now for gc indexing:
       procbounds%gridbeg = seqIndexList(1)
       procbounds%gridend = seqIndexList(lsize)
       procbounds%im      = lsize

       write(*,*) "bounds: ", de, seqIndexList(1), seqIndexList(lsize), lsize
    
       ! ---------------------
       ! Realize the actively coupled fields
       ! ---------------------
       call realize_fields(gcomp, mesh=mesh_input, flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    elseif ( trim(ctrl_init%grid) == 'CS' ) then ! cube sphere mosaic

       geomtype = ESMF_GEOMTYPE_GRID
       
       gridfile = "grid_spec.nc" ! default                                                                                              
       if (field_exist("INPUT/grid_spec.nc", "atm_mosaic_file")) then
          call read_data("INPUT/grid_spec.nc", "atm_mosaic_file", gridfile)
       endif

       do tl=1,6
          decomptile(1,tl) = ctrl_init%layout(1)
          decomptile(2,tl) = ctrl_init%layout(2)
          decompflagPTile(:,tl) = (/ESMF_DECOMP_SYMMEDGEMAX,ESMF_DECOMP_SYMMEDGEMAX/)
       enddo
       lndGrid = ESMF_GridCreateMosaic(filename="INPUT/"//trim(gridfile),     &
            regDecompPTile=decomptile,tileFilePath="INPUT/",                   &
            decompflagPTile=decompflagPTile,                                   &
            staggerlocList=(/ESMF_STAGGERLOC_CENTER, ESMF_STAGGERLOC_CORNER/), &
            name='lnd_grid', rc=rc)

       ! ---------------------
       ! Realize the actively coupled fields
       ! ---------------------
       call realize_fields(gcomp, grid=lndGrid, flds_scalar_name=flds_scalar_name, &
            flds_scalar_num=flds_scalar_num, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return


       ! ! ---------------------
       ! ! Define Blocking
       ! ! ---------------------
       ! call define_blocks_packed ('land_model', Lnd_block, isc, iec, jsc, jec, nlev, &
       !      blocksize, block_message)

    else
       write(*,*) 'Do not know if imports/exports should be 1D or 2D from ctrl_init%grid'                                           
       stop

    end if


    ! set variable geomtype for import/export module
    call ie_set_geomtype(geomtype)
    
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


    !!! add for testing,
    ! call ESMF_MeshGet(mesh_input, elementCount=elementCount, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! write(*,*) "JP B. element count = ", elementCount


  end subroutine InitializeRealize

  !===============================================================================
  subroutine ModelAdvance(gcomp, rc)

    !use lsm_noah, only: lsm_noah_run
    !use noah_loop, only: noah_loop_run
    use noah_driver, only: noah_block_run, noah_model, ctrl_init, noah_loop_drv
    use import_fields_mod, only: write_import_field, import_allfields_am, export_allfields
    use proc_bounds, only : procbounds
    
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
    real(r8), pointer          :: dataPtr_e(:)
    character(len=*),parameter :: fldname_e="foo_lnd2atmfield"
    character(len=80), allocatable :: flds(:)
    character(len=80)   :: fldname
    integer                   :: n

    ! tmp debug
    integer                 :: i, de, gridbeg, gridend, im
    real(r8)                :: foodata(procbounds%im)

    ! query the Component for its clock, importState and exportState
    ! call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
    !      exportState=exportState, rc=rc)
    call NUOPC_ModelGet(gcomp, modelClock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    ! Commenting out, because won't currently work with tiles
    ! ! write out imports. Should make this optional, and put into a routine
    
    ! allocate(flds(7))
    ! flds = (/'Faxa_lwdn  '    , 'Faxa_swndr '   , 'Faxa_swvdr '   , 'Faxa_swndf ' , 'Faxa_swvdf ', &
    !      'Faxa_rain  '    , 'Faxa_snow  ' /)

    ! do n = 1,size(flds)
    !    fldname = trim(flds(n))
    !    call write_import_field(importState, fldname, rc)
    ! end do
    ! deallocate(flds)

    ! allocate(flds(27))
    ! flds=(/ &
    !      'Faxa_soiltyp     ', &
    !      'Faxa_vegtype     ', &
    !      'Faxa_sigmaf      ', &
    !      'Faxa_sfcemis     ', &
    !      'Faxa_dlwflx      ', &
    !      'Faxa_dswsfc      ', &
    !      'inst_down_sw_flx ', &
    !      'Faxa_snet        ', &
    !      'Faxa_tg3         ', &
    !      'Faxa_cm          ', &
    !      'Faxa_ch          ', &
    !      'Faxa_prsl1       ', &
    !      'Faxa_prslki      ', &
    !      'Faxa_zf          ', &
    !      'Faxa_land        ', &
    !      'Faxa_slopetyp    ', &
    !      'Faxa_shdmin      ', &
    !      'Faxa_shdmax      ', &
    !      'Faxa_snoalb      ', &
    !      'Faxa_sfalb       ', &
    !      'Faxa_bexppert    ', &
    !      'Faxa_xlaipert    ', &
    !      'Faxa_vegfpert    ', &
    !      'Faxa_tsurf       ', &
    !      'Faxa_wind        ', &
    !      'Faxa_ps          ', &
    !      'Faxa_t1          ', &
    !      'Faxa_q1          ', &
    !      'Faxa_z0rl        ', &
    !      'Faxa_canopy      ', &
    !      'Faxa_tprcp       ', &
    !      'Faxa_weasd       ', &
    !      'Faxa_ustar       ' &
    !      /)


    ! ! tmp
    ! !write(*,*) 'procbound test:', procbounds%de, procbounds%gridbeg, procbounds%gridend

    ! do n = 1,size(flds)
    !    fldname = trim(flds(n))
    !    call write_import_field(importState, fldname, rc)
    ! end do
    ! deallocate(flds)


    ! end test tmp

    call import_allfields_am(importState, procbounds, noah_model, ctrl_init, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! ! tmp workaround to run first time step with simple nems.configure sequence, without having land restart read
    if (noah_model%control%first_time) then
       write(*,*) 'lnd_comp: skipping first time step'
    else
    
       ! run model
       !call noah_block_run(procbounds, noah_model) !! IF using blocking
       call noah_loop_drv(procbounds, noah_model)
       
       ! tmp test
       ! im = procbounds%im
       ! gridbeg = procbounds%gridbeg
       ! gridend = procbounds%gridend
       ! foodata(1:im) = noah_model%model%foo_atm2lndfield(gridbeg:gridend)
       ! do i = 1,im
       !    write(*,*) 'MA1: ', de, gridbeg,gridend, size(foodata), foodata(i)
       ! end do

    end if ! first_time

    ! TODO: move within first_time step if else statement
    call export_allfields(exportState, procbounds, noah_model, ctrl_init, rc)

   
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

    ! tmp workaround to run first time step with simple nems.configure sequence, without having land restart read 
    if (noah_model%control%first_time) then
       noah_model%control%first_time = .false.
    endif
    
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
