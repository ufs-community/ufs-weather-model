module import_fields_mod

  use ESMF,                   only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF,                      only: ESMF_Field, ESMF_FieldGet, ESMF_FieldCreate
  use ESMF,                   only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF,                   only : ESMF_GEOMTYPE_FLAG, ESMF_GEOMTYPE_GRID, ESMF_GEOMTYPE_MESH
  use ESMF,                   only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF,                   only : operator(/=), operator(==)
  use ESMF,                   only : ESMF_FieldWrite
  use shr_kind_mod,           only : r8 => shr_kind_r8
  use nuopc_shr_methods,      only : chkerr

  use proc_bounds,            only: procbounds_type

  implicit none

  type(ESMF_GeomType_Flag) :: geomtype
  
  private

  ! interface import_field
  !    module procedure import_field_1Dreal
  !    module procedure import_field_1Dint
  !    module procedure import_field_1Dlog
  !    module procedure import_field_2Dreal
  !    module procedure import_field_2Dint
  !    module procedure import_field_2Dlog
  ! end interface import_field

  interface export_field
     module procedure export_field_1D
     module procedure export_field_2D     
  end interface export_field

  ! from MOM cap
  interface State_GetFldPtr
     module procedure State_GetFldPtr_1d
     module procedure State_GetFldPtr_2d
  end interface State_GetFldPtr

  ! This is needed until restarts are read in
  interface State_GetImport
     module procedure State_GetImport_Real
     module procedure State_GetImport_Int
     module procedure State_GetImport_Log
  end interface State_GetImport

  public write_import_field, import_allfields_am, export_allfields, ie_set_geomtype
  
    character(*),parameter :: u_FILE_u = &
       __FILE__
  
contains

  ! Sets module variable geometry type
  subroutine ie_set_geomtype(geomtype_in)
    type(ESMF_GeomType_Flag), intent(in)    :: geomtype_in !< ESMF type describing type of
    !! geometry (mesh or grid)

    geomtype = geomtype_in

  end subroutine ie_set_geomtype

  
  ! subroutine import_field_1Dreal(State_i, fldname_i, out_data, rc)

  !   use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
  !   use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
  !   use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

  !   type(ESMF_State),             intent(in)    :: State_i
  !   character(len=*),             intent(in)    :: fldname_i
  !   real(r8),                     intent(inout) :: out_data(:)
  !   integer,                      intent(out)   :: rc

  !   ! local
  !   real(r8), pointer           :: fldPtr1d_i(:)
  !   integer                     :: i

    
  !   if (fldchk(State_i,trim(fldname_i))) then
  !      call state_getfldptr(State_i,trim(fldname_i), fldptr1d=fldptr1d_i, rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
  !      do i=1,size(fldptr1d_i)
  !         out_data(i) = fldptr1d_i(i)
  !      enddo

  !   else
  !      write(*,*) 'Field ', trim(fldname_i), ' NOT found'

  !   endif

    
  ! end subroutine import_field_1Dreal

  ! subroutine import_field_2Dreal(State_i, fldname_i, out_data, rc)

  !   use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
  !   use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
  !   use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

  !   type(ESMF_State),             intent(in)    :: State_i
  !   character(len=*),             intent(in)    :: fldname_i
  !   real(r8),                     intent(inout) :: out_data(:,:)
  !   integer,                      intent(out)   :: rc

  !   ! local
  !   real(r8), pointer           :: fldPtr2d_i(:,:)
  !   integer                     :: i, j

    
  !   if (fldchk(State_i,trim(fldname_i))) then
  !      call state_getfldptr(State_i,trim(fldname_i), fldptr2d=fldptr2d_i, rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
  !      do i=1,size(fldptr2d_i, dim=1)
  !         do j=1,size(fldptr2d_i, dim=2)
  !            out_data(i,j) = fldptr2d_i(i,j)
  !         enddo
  !      enddo

  !   else
  !      write(*,*) 'Field ', trim(fldname_i), ' NOT found'

  !   endif

    
  ! end subroutine import_field_2Dreal

  ! subroutine import_field_1Dint(State_i, fldname_i, out_data, rc)

  !   use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
  !   use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
  !   use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

  !   type(ESMF_State),            intent(in)    :: State_i
  !   character(len=*),            intent(in)    :: fldname_i
  !   integer,                     intent(inout) :: out_data(:)
  !   integer,                     intent(out)   :: rc

  !   ! local
  !   real(r8), pointer           :: fldPtr1d_i(:)
  !   integer                     :: i

    
  !   if (fldchk(State_i,trim(fldname_i))) then
  !      call state_getfldptr(State_i,trim(fldname_i), fldptr1d_i, rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
  !      do i=1,size(fldptr1d_i)
  !         out_data(i) = fldptr1d_i(i)
  !      enddo

  !   else
  !      write(*,*) 'Field ', trim(fldname_i), ' NOT found'

  !   endif

    
  ! end subroutine import_field_1Dint

  ! subroutine import_field_2Dint(State_i, fldname_i, out_data, rc)

  !   use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
  !   use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
  !   use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

  !   type(ESMF_State),            intent(in)    :: State_i
  !   character(len=*),            intent(in)    :: fldname_i
  !   integer,                     intent(inout) :: out_data(:,:)
  !   integer,                     intent(out)   :: rc

  !   ! local
  !   real(r8), pointer           :: fldPtr2d_i(:,:)
  !   integer                     :: i,j

    
  !   if (fldchk(State_i,trim(fldname_i))) then
  !      call state_getfldptr(State_i,trim(fldname_i), fldptr2d=fldptr2d_i, rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
  !      do i=1,size(fldptr2d_i, dim=1)
  !         do j=1,size(fldptr2d_i, dim=2)
  !            out_data(i,j) = fldptr2d_i(i,j)
  !         enddo
  !      enddo

  !   else
  !      write(*,*) 'Field ', trim(fldname_i), ' NOT found'

  !   endif

    
  ! end subroutine import_field_2Dint
  
  ! subroutine import_field_1Dlog(State_i, fldname_i, out_data, rc)

  !   use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
  !   use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
  !   use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

  !   type(ESMF_State),            intent(in)    :: State_i
  !   character(len=*),            intent(in)    :: fldname_i
  !   logical,                     intent(inout) :: out_data(:)
  !   integer,                     intent(out)   :: rc

  !   ! local
  !   real(r8), pointer           :: fldPtr1d_i(:)
  !   integer                     :: i

    
  !   if (fldchk(State_i,trim(fldname_i))) then
  !      call state_getfldptr(State_i,trim(fldname_i), fldptr1d_i, rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
  !      do i=1,size(fldptr1d_i)
  !         out_data(i) = fldptr1d_i(i)
  !      enddo

  !   else
  !      write(*,*) 'Field ', trim(fldname_i), ' NOT found'

  !   endif

    
  ! end subroutine import_field_1Dlog

  ! subroutine import_field_2Dlog(State_i, fldname_i, out_data, rc)

  !   use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
  !   use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
  !   use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

  !   type(ESMF_State),            intent(in)    :: State_i
  !   character(len=*),            intent(in)    :: fldname_i
  !   logical,                     intent(inout) :: out_data(:,:)
  !   integer,                     intent(out)   :: rc

  !   ! local
  !   real(r8), pointer           :: fldPtr2d_i(:,:)
  !   integer                     :: i,j

    
  !   if (fldchk(State_i,trim(fldname_i))) then
  !      call state_getfldptr(State_i,trim(fldname_i), fldptr2d=fldptr2d_i, rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
  !      do i=1,size(fldptr2d_i, dim=1)
  !         do j=1,size(fldptr2d_i, dim=2)
  !            out_data(i,j) = fldptr2d_i(i,j)
  !         enddo
  !      enddo

  !   else
  !      write(*,*) 'Field ', trim(fldname_i), ' NOT found'

  !   endif

    
  ! end subroutine import_field_2Dlog
  

  subroutine write_import_field(State_i, fldname_i, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State_i
    !type(ESMF_State),             intent(in)    :: State_e
    character(len=*),             intent(in)    :: fldname_i
    !character(len=*),             intent(in)    :: fldname_e
    !real(r8), pointer, optional , intent(out)   :: fldptr1d_i(:)
    !real(r8), pointer, optional , intent(out)   :: fldptr1d_e(:)
    !real(r8),                     intent(inout) :: out_data(:)
    integer,                      intent(out)   :: rc

    ! local variables
    real(r8), pointer           :: fldPtr1d_i(:)
    !real(r8), pointer           :: fldPtr1d_e(:)
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    integer                     :: i,j


    if (fldchk(State_i,trim(fldname_i))) then
       !write(*,*) 'Field ', trim(fldname_i), ' found'
       call state_getfldptr(State_i,trim(fldname_i), fldptr1d_i, rc=rc)

       ! do j=1,size(fldptr1d_i)
       !    out_data(j) = fldptr1d_i(j)
       ! enddo

    else
       write(*,*) 'Field ', trim(fldname_i), ' NOT found'

    endif

    call ESMF_StateGet(State_i, itemName=trim(fldname_i), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Write out to netcdf
    call ESMF_FieldWrite(lfield, fileName='lnd_import.'//trim(fldname_i)//'.nc', &
         timeslice=1, overwrite=.true., rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    
  end subroutine write_import_field


  subroutine import_allfields_am(State_i, procbounds, noah_model, ctrl_init, rc)
    use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    use noah_type_mod, only: noah_type
    use proc_bounds,   only: control_init_type

    ! ----------------------------------------------
    ! Get pointers to internal land model variables, for multiple fields
    ! ----------------------------------------------
    type(ESMF_State),        intent(in)    :: State_i
    type(procbounds_type),   intent(in)    :: procbounds
    type(noah_type),         intent(inout) :: noah_model ! land model's variable type
    type(control_init_type), intent(in)    :: ctrl_init
    integer,                 intent(out)   :: rc

    ! local
    !type(ESMF_State)        :: importState
    character(len=25)       :: fldname_i
    real(r8), pointer       :: fldPtr1d_i(:)
    integer                 :: g, im, nb
    integer                 :: isc,iec, jsc, jec

    
    isc = noah_model%control%isc 
    iec = noah_model%control%iec 
    jsc = noah_model%control%jsc 
    jec = noah_model%control%jec
    
    ! if ( trim(ctrl_init%grid) == '1Dmesh' ) then ! 
    !    geomtype = ESMF_GEOMTYPE_MESH
    !    write (*,*) 'land import setting geomtype is MESH'
    ! elseif ( trim(ctrl_init%grid) == 'CS' ) then
    !    geomtype = ESMF_GEOMTYPE_GRID
    !    write (*,*) 'land import setting geomtype is GRID'
    ! end if

    call state_getimport(State_i, 'foo_atm2lndfield', isc, iec, jsc, jec, noah_model%model%foo_atm2lndfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! write(*,*) 'i isc, iec, jsc, jec: ', isc, iec, jsc, jec
    ! write(*,*) 'foo_atm2lndfield: ', noah_model%model%foo_atm2lndfield

    call state_getimport(State_i, 'Faxa_soiltyp', isc, iec, jsc, jec, noah_model%model%soiltyp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_vegtype', isc, iec, jsc, jec, noah_model%model%vegtype, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_sigmaf', isc, iec, jsc, jec, noah_model%model%sigmaf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_sfcemis', isc, iec, jsc, jec, noah_model%model%sfcemis, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_dlwflx', isc, iec, jsc, jec, noah_model%model%dlwflx, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_dswsfc', isc, iec, jsc, jec, noah_model%model%dswsfc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'inst_down_sw_flx', isc, iec, jsc, jec, noah_model%model%dswsfci, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_snet', isc, iec, jsc, jec, noah_model%model%snet, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_tg3', isc, iec, jsc, jec, noah_model%model%tg3, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_cm', isc, iec, jsc, jec, noah_model%model%cm, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_ch', isc, iec, jsc, jec, noah_model%model%ch, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_prsl1', isc, iec, jsc, jec, noah_model%model%prsl1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_prslki', isc, iec, jsc, jec, noah_model%model%prslki, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_zf', isc, iec, jsc, jec, noah_model%model%zf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_land', isc, iec, jsc, jec, noah_model%model%land, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_slopetyp', isc, iec, jsc, jec, noah_model%model%slopetyp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_shdmin', isc, iec, jsc, jec, noah_model%model%shdmin, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_shdmax', isc, iec, jsc, jec, noah_model%model%shdmax, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_snoalb', isc, iec, jsc, jec, noah_model%model%snoalb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_sfalb', isc, iec, jsc, jec, noah_model%model%sfalb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_bexppert', isc, iec, jsc, jec, noah_model%model%bexppert, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_xlaipert', isc, iec, jsc, jec, noah_model%model%xlaipert, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_vegfpert', isc, iec, jsc, jec, noah_model%model%vegfpert, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call state_getimport(State_i, 'Faxa_ps', isc, iec, jsc, jec, noah_model%model%ps, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_t1', isc, iec, jsc, jec, noah_model%model%t1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_q1', isc, iec, jsc, jec, noah_model%model%q1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    
    ! call state_getimport(State_i, 'Faxa_prsik1', isc, iec, jsc, jec, noah_model%model%prsik1, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_sfalb', isc, iec, jsc, jec, noah_model%model%sfalb, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_weasd', isc, iec, jsc, jec, noah_model%model%weasd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_snwdph', isc, iec, jsc, jec, noah_model%model%snwdph, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_tskin', isc, iec, jsc, jec, noah_model%model%tskin, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_tprcp', isc, iec, jsc, jec, noah_model%model%tprcp, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_srflag', isc, iec, jsc, jec, noah_model%model%srflag, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(State_i, 'Faxa_smc', isc, iec, jsc, jec, noah_model%model%smc, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(State_i, 'Faxa_stc', isc, iec, jsc, jec, noah_model%model%stc, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(State_i, 'Faxa_slc', isc, iec, jsc, jec, noah_model%model%slc, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_canopy', isc, iec, jsc, jec, noah_model%model%canopy, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_trans', isc, iec, jsc, jec, noah_model%model%trans, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_tsurf', isc, iec, jsc, jec, noah_model%model%tsurf, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_z0rl', isc, iec, jsc, jec, noah_model%model%z0rl, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(State_i, 'Faxa_z0pert', isc, iec, jsc, jec, noah_model%model%z0pert, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_getimport(State_i, 'Faxa_ztpert', isc, iec, jsc, jec, noah_model%model%ztpert, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_ustar', isc, iec, jsc, jec, noah_model%model%ustar, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_wind', isc, iec, jsc, jec, noah_model%model%wind, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    call state_getimport(State_i, 'Faxa_albdvis_lnd', isc, iec, jsc, jec, noah_model%model%albdvis_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_albdnir_lnd', isc, iec, jsc, jec, noah_model%model%albdnir_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_albivis_lnd', isc, iec, jsc, jec, noah_model%model%albivis_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_albinir_lnd', isc, iec, jsc, jec, noah_model%model%albinir_lnd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_adjvisbmd', isc, iec, jsc, jec, noah_model%model%adjvisbmd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_adjnirbmd', isc, iec, jsc, jec, noah_model%model%adjnirbmd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_adjvisdfd', isc, iec, jsc, jec, noah_model%model%adjvisdfd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_adjnirdfd', isc, iec, jsc, jec, noah_model%model%adjnirdfd, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_prslk1', isc, iec, jsc, jec, noah_model%model%prslk1, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_getimport(State_i, 'Faxa_garea', isc, iec, jsc, jec, noah_model%model%garea, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    

    
  end subroutine import_allfields_am


  ! subroutine import_allfields(State_i, procbounds, noah_model, ctrl_init, rc)

  !   use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
  !   use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
  !   use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

  !   use noah_type_mod, only: noah_type
  !   use proc_bounds,   only: control_init_type
    
  !   ! ----------------------------------------------
  !   ! Get pointers to internal land model variables, for multiple fields
  !   ! ----------------------------------------------
  !   type(ESMF_State),        intent(in)    :: State_i
  !   type(procbounds_type),   intent(in)    :: procbounds 
  !   type(noah_type),         intent(inout) :: noah_model(:) ! land model's variable type
  !   type(control_init_type), intent(in)    :: ctrl_init
  !   integer,                 intent(out)   :: rc

  !   ! local
  !   !type(ESMF_State)        :: importState
  !   character(len=25)       :: fldname_i
  !   real(r8), pointer       :: fldPtr1d_i(:)
  !   integer                 :: g, im, nb
    

    
  !   do nb=1, size(noah_model)

  !      ! TODO: move out of loop
  !      !if (trim(ctrl_init%grid) == '1Dmesh') then
  !      !   im  = procbounds%im
  !      !else if(trim(ctrl_init%grid) == 'CS') then ! blocking
  !         im = noah_model(nb)%control%blksz
  !      !else
  !      !   write(*,*) 'Do not know if imports/exports should be 1D or 2D from ctrl_init%grid'
  !      !   stop
  !      !endif

  !      ! tmp for testing
  !      call import_field(State_i, 'foo_atm2lndfield', noah_model(nb)%model%foo_atm2lndfield(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !      call import_field(State_i, 'Faxa_soiltyp', noah_model(nb)%model%soiltyp(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_vegtype', noah_model(nb)%model%vegtype(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_sigmaf', noah_model(nb)%model%sigmaf(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_sfcemis', noah_model(nb)%model%sfcemis(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_dlwflx', noah_model(nb)%model%dlwflx(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_dswsfc', noah_model(nb)%model%dswsfc(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return    
  !      call import_field(State_i, 'inst_down_sw_flx', noah_model(nb)%model%dswsfci(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return    
  !      call import_field(State_i, 'Faxa_snet', noah_model(nb)%model%snet(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_tg3', noah_model(nb)%model%tg3(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_cm', noah_model(nb)%model%cm(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_ch', noah_model(nb)%model%ch(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_prsl1', noah_model(nb)%model%prsl1(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_prslki', noah_model(nb)%model%prslki(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_zf', noah_model(nb)%model%zf(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_land', noah_model(nb)%model%land(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_slopetyp', noah_model(nb)%model%slopetyp(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_shdmin', noah_model(nb)%model%shdmin(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_shdmax', noah_model(nb)%model%shdmax(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_snoalb', noah_model(nb)%model%snoalb(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_sfalb', noah_model(nb)%model%sfalb(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_bexppert', noah_model(nb)%model%bexppert(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_xlaipert', noah_model(nb)%model%xlaipert(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_vegfpert', noah_model(nb)%model%vegfpert(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !      call import_field(State_i, 'Faxa_ps', noah_model(nb)%model%ps(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_t1', noah_model(nb)%model%t1(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_q1', noah_model(nb)%model%q1(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return

  !      ! call import_field(State_i, 'Faxa_prsik1', noah_model(nb)%model%prsik1(1:im), rc=rc)
  !      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_sfalb', noah_model(nb)%model%sfalb(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_weasd', noah_model(nb)%model%weasd(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_snwdph', noah_model(nb)%model%snwdph(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_tskin', noah_model(nb)%model%tskin(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_tprcp', noah_model(nb)%model%tprcp(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_srflag', noah_model(nb)%model%srflag(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      ! call import_field(State_i, 'Faxa_smc', noah_model(nb)%model%smc(1:im), rc=rc)
  !      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      ! call import_field(State_i, 'Faxa_stc', noah_model(nb)%model%stc(1:im), rc=rc)
  !      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      ! call import_field(State_i, 'Faxa_slc', noah_model(nb)%model%slc(1:im), rc=rc)
  !      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_canopy', noah_model(nb)%model%canopy(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_trans', noah_model(nb)%model%trans(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_tsurf', noah_model(nb)%model%tsurf(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_z0rl', noah_model(nb)%model%z0rl(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      ! call import_field(State_i, 'Faxa_z0pert', noah_model(nb)%model%z0pert(1:im), rc=rc)
  !      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      ! call import_field(State_i, 'Faxa_ztpert', noah_model(nb)%model%ztpert(1:im), rc=rc)
  !      ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_ustar', noah_model(nb)%model%ustar(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !      call import_field(State_i, 'Faxa_wind', noah_model(nb)%model%wind(1:im), rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return    

  !   end do
    
  ! end subroutine import_allfields
  

  subroutine export_allfields(State_o, procbounds, noah_model, ctrl_init, rc)

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    use noah_type_mod, only: noah_type
    use proc_bounds,   only: control_init_type
    
    ! ----------------------------------------------
    ! Get pointers to internal land model variables, for multiple fields
    ! ----------------------------------------------
    type(ESMF_State),        intent(in)    :: State_o
    type(procbounds_type),   intent(in)    :: procbounds
    type(noah_type),         intent(in)    :: noah_model ! land model's variable type
    type(control_init_type), intent(in)    :: ctrl_init
    integer,                 intent(out)   :: rc

    ! local
    !type(ESMF_State)        :: importState
    character(len=25)       :: fldname_i
    real(r8), pointer       :: fldPtr1d_i(:)
    integer                 :: g, gridbeg, gridend,im, nb


    ! do nb=1, size(noah_model)
    
    !    ! TODO: move out of loop
    !    if (trim(ctrl_init%grid) == '1Dmesh') then
          im  = procbounds%im
       ! else if(trim(ctrl_init%grid) == 'CS') then ! blocking
       !    im = noah_model(nb)%control%blksz
       ! else
       !    write(*,*) 'Do not know if imports/exports should be 1D or 2D from ctrl_init%grid'
       !    stop
       ! endif

       ! inouts   
       call export_field(State_o, 'Fall_weasd', noah_model%model%weasd(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_snwdph', noah_model%model%snwdph(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_tskin', noah_model%model%tskin(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_tprcp', noah_model%model%tprcp(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_srflag', noah_model%model%srflag(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call export_field(State_o, 'Fall_smc', noah_model%model%smc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call export_field(State_o, 'Fall_stc', noah_model%model%stc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call export_field(State_o, 'Fall_slc', noah_model%model%slc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_canopy', noah_model%model%canopy(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_trans', noah_model%model%trans(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_tsurf', noah_model%model%tsurf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_z0rl', noah_model%model%z0rl(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! noahouts
       call export_field(State_o, 'Fall_sncovr1', noah_model%model%sncovr1(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_qsurf', noah_model%model%qsurf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_gflux', noah_model%model%gflux(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_drain', noah_model%model%drain(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_evap', noah_model%model%evap(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_hflx', noah_model%model%hflx(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_ep', noah_model%model%ep(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_runoff', noah_model%model%runoff(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_cmm', noah_model%model%cmm(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_chh', noah_model%model%chh(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_evbs', noah_model%model%evbs(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_evcw', noah_model%model%evcw(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_sbsno', noah_model%model%sbsno(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_snowc', noah_model%model%snowc(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_stm', noah_model%model%stm(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_snohf', noah_model%model%snohf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_smcwlt2', noah_model%model%smcwlt2(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_smcref2', noah_model%model%smcref2(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_wet1', noah_model%model%wet1(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! diffouts
       call export_field(State_o, 'Fall_rb_lnd', noah_model%model%rb_lnd(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_fm_lnd', noah_model%model%fm_lnd(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_fh_lnd', noah_model%model%fh_lnd(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_fm10_lnd', noah_model%model%fm10_lnd(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_fh2_lnd', noah_model%model%fh2_lnd(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call export_field(State_o, 'Fall_stress', noah_model%model%stress(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !end do
  end subroutine export_allfields

  
!===============================================================================

  subroutine export_field_1D(state, fldname, outdata, rc)

    ! fill in export data for 1d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in) :: state
    character(len=*) , intent(in) :: fldname
    real(r8)         , intent(in) :: outdata(:)
    integer          , intent(out):: rc

    ! local variables
    real(r8), pointer :: fldPtr1d(:)
    integer           :: g
    ! ----------------------------------------------

    !nullify(fldPtr1d) ! do I need to do this?
    !allocate( fldPtr1d(size(outdata)) ) ! do I need to do this?
    call state_getfldptr(state, trim(fldname), fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    fldptr1d(:) = 0._r8
    do g = 1,size(outdata)
       fldptr1d(g) = outdata(g)
    end do

  end subroutine export_field_1D
  
  !===============================================================================

  subroutine export_field_2D(state, fldname, outdata, rc)

    ! fill in export data for 2d field

    use ESMF, only : ESMF_LOGERR_PASSTHRU, ESMF_END_ABORT, ESMF_LogFoundError
    use ESMF, only : ESMF_Finalize

    ! input/output variabes
    type(ESMF_State) , intent(in) :: state
    character(len=*) , intent(in) :: fldname
    real(r8)         , intent(in) :: outdata(:,:)
    integer          , intent(out):: rc

    ! local variables
    real(r8), pointer :: fldPtr2d(:,:)
    integer           :: i,j
    ! ----------------------------------------------

   
    call state_getfldptr(state, trim(fldname), fldptr2d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    fldptr2d(:,:) = 0._r8
    do i=1,size(fldptr2d, dim=1)
       do j=1,size(fldptr2d, dim=2)
          fldptr2d(i,j) = outdata(i,j)
       enddo
    enddo

  end subroutine export_field_2D

  ! Get field pointer 1D
  subroutine State_GetFldPtr_1d(State, fldname, fldptr, rc)
    type(ESMF_State)            , intent(in)  :: State    !< ESMF state
    character(len=*)            , intent(in)  :: fldname  !< Field name
    real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:)!< Pointer to the 1D field
    integer, optional           , intent(out) :: rc       !< Return code

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(Noah_cap:State_GetFldPtr)'

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr_1d

  ! Get field pointer 2D
  subroutine State_GetFldPtr_2d(State, fldname, fldptr, rc)
    type(ESMF_State)            , intent(in)  :: State      !< ESMF state
    character(len=*)            , intent(in)  :: fldname    !< Field name
    real(ESMF_KIND_R8), pointer , intent(in)  :: fldptr(:,:)!< Pointer to the 2D field
    integer, optional           , intent(out) :: rc         !< Return code

    ! local variables
    type(ESMF_Field) :: lfield
    integer :: lrc
    character(len=*),parameter :: subname='(Noah_cap:State_GetFldPtr)'

    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=lrc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldGet(lfield, farrayPtr=fldptr, rc=lrc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (present(rc)) rc = lrc

  end subroutine State_GetFldPtr_2d

  
  ! !===============================================================================
  ! subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

  !   ! ----------------------------------------------
  !   ! Get pointer to a state field
  !   ! ----------------------------------------------

  !   use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
  !   use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
  !   use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

  !   use shr_sys_mod, only : shr_sys_abort
    
  !   ! input/output variables
  !   type(ESMF_State),             intent(in)    :: State
  !   character(len=*),             intent(in)    :: fldname
  !   real(r8), pointer, optional , intent(out)   :: fldptr1d(:)
  !   real(r8), pointer, optional , intent(out)   :: fldptr2d(:,:)
  !   integer,                      intent(out)   :: rc

  !   ! local variables
  !   type(ESMF_FieldStatus_Flag) :: status
  !   type(ESMF_Field)            :: lfield
  !   character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
  !   type(ESMF_StateItem_Flag)   :: itemFlag
  !   ! ----------------------------------------------

  !   rc = ESMF_SUCCESS
  !   call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
  !   !if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !   if (present(fldptr1d)) then
  !      call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !   else if (present(fldptr2d)) then
  !      call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
  !      if (ChkErr(rc,__LINE__,u_FILE_u)) return
  !   else
  !      call shr_sys_abort("either fldptr1d or fldptr2d must be an input argument")
  !   end if
  ! end subroutine state_getfldptr

  !===============================================================================
  logical function fldchk(state, fldname)
    ! ----------------------------------------------
    ! Determine if field with fldname is in the input state
    ! ----------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(in)  :: state
    character(len=*), intent(in)  :: fldname

    ! local variables
    type(ESMF_StateItem_Flag)   :: itemFlag
    ! ----------------------------------------------
    call ESMF_StateGet(state, trim(fldname), itemFlag)
    if (itemflag /= ESMF_STATEITEM_NOTFOUND) then
       fldchk = .true.
    else
       fldchk = .false.
    endif
  end function fldchk


  ! From MOM: Map import state field to output array
  subroutine State_GetImport_Real(state, fldname, isc, iec, jsc, jec, output, do_sum, areacor, rc)
    type(ESMF_State)    , intent(in)    :: state   !< ESMF state
    character(len=*)    , intent(in)    :: fldname !< Field name
    integer             , intent(in)    :: isc     !< The start i-index of cell centers within
    !! the computational domain
    integer             , intent(in)    :: iec     !< The end i-index of cell centers within the
    !! computational domain
    integer             , intent(in)    :: jsc     !< The start j-index of cell centers within
    !! the computational domain
    integer             , intent(in)    :: jec     !< The end j-index of cell centers within
    !! the computational domain
    !real (ESMF_KIND_R8) , intent(inout) :: output(isc:iec,jsc:jec)!< Output 2D array
    real (ESMF_KIND_R8) , intent(inout) :: output((iec-isc+1)*(jec-isc+1))!< Output 1D array, size im
    logical, optional   , intent(in)    :: do_sum  !< If true, sums the data
    real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
    !! applicable to meshes
    integer             , intent(out)   :: rc      !< Return code

    ! local variables
    type(ESMF_StateItem_Flag)     :: itemFlag
    integer                       :: n, i, j, i1, j1
    integer                       :: lbnd1,lbnd2
    real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
    character(len=*)  , parameter :: subname='(Noah import_fields:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, trim(fldname), itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

       if (geomtype == ESMF_GEOMTYPE_MESH) then

          ! get field pointer
          call state_getfldptr(state, trim(fldname), dataptr1d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! determine output array and apply area correction if present
          n = 0
          do j = jsc,jec
             do i = isc,iec
                n = n + 1
                if (present(do_sum)) then
                   if (present(areacor)) then
                      !output(i,j)  = output(i,j) + dataPtr1d(n) * areacor(n)
                      output(n)  = output(n) + dataPtr1d(n) * areacor(n)
                   else
                      !output(i,j)  = output(i,j) + dataPtr1d(n)
                      output(n)  = output(n) + dataPtr1d(n)
                   end if
                else
                   if (present(areacor)) then
                      !output(i,j)  = dataPtr1d(n) * areacor(n)
                      output(n)  = dataPtr1d(n) * areacor(n)
                   else
                      !output(i,j)  = dataPtr1d(n)
                      output(n)  = dataPtr1d(n)
                   end if
                endif
             enddo
          enddo

       else if (geomtype == ESMF_GEOMTYPE_GRID) then

          call state_getfldptr(state, trim(fldname), dataptr2d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          lbnd1 = lbound(dataPtr2d,1)
          lbnd2 = lbound(dataPtr2d,2)

          n = 0 
          do j = jsc, jec
             j1 = j + lbnd2 - jsc
             do i = isc, iec
                i1 = i + lbnd1 - isc
                n = n +1
                if (present(do_sum)) then
                   !output(i,j) = output(i,j) + dataPtr2d(i1,j1)
                   output(n) = output(n) + dataPtr2d(i1,j1)
                else
                   output(n) = dataPtr2d(i1,j1)
                   !output(i,j) = dataPtr2d(i1,j1)
                endif
             enddo
          enddo
       endif
    endif

  end subroutine State_GetImport_Real

!!!!!!!!!!!!!!! Following aren't needed after restarts are read in:
  ! From MOM: Map import state field to output array
  subroutine State_GetImport_Int(state, fldname, isc, iec, jsc, jec, output, do_sum, areacor, rc)
    type(ESMF_State)    , intent(in)    :: state   !< ESMF state
    character(len=*)    , intent(in)    :: fldname !< Field name
    integer             , intent(in)    :: isc     !< The start i-index of cell centers within
    !! the computational domain
    integer             , intent(in)    :: iec     !< The end i-index of cell centers within the
    !! computational domain
    integer             , intent(in)    :: jsc     !< The start j-index of cell centers within
    !! the computational domain
    integer             , intent(in)    :: jec     !< The end j-index of cell centers within
    !! the computational domain
    !real (ESMF_KIND_R8) , intent(inout) :: output(isc:iec,jsc:jec)!< Output 2D array
    integer             , intent(inout) :: output((iec-isc+1)*(jec-isc+1))!< Output 1D array
    logical, optional   , intent(in)    :: do_sum  !< If true, sums the data
    real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
    !! applicable to meshes
    integer             , intent(out)   :: rc      !< Return code

    ! local variables
    type(ESMF_StateItem_Flag)     :: itemFlag
    integer                       :: n, i, j, i1, j1
    integer                       :: lbnd1,lbnd2
    real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
    character(len=*)  , parameter :: subname='(Noah import_fields:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, trim(fldname), itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

       if (geomtype == ESMF_GEOMTYPE_MESH) then

          ! get field pointer
          call state_getfldptr(state, trim(fldname), dataptr1d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! determine output array and apply area correction if present
          n = 0
          do j = jsc,jec
             do i = isc,iec
                n = n + 1
                if (present(do_sum)) then
                   if (present(areacor)) then
                      !output(i,j)  = output(i,j) + dataPtr1d(n) * areacor(n)
                      output(n)  = output(n) + dataPtr1d(n) * areacor(n)
                   else
                      !output(i,j)  = output(i,j) + dataPtr1d(n)
                      output(n)  = output(n) + dataPtr1d(n)
                   end if
                else
                   if (present(areacor)) then
                      !output(i,j)  = dataPtr1d(n) * areacor(n)
                      output(n)  = dataPtr1d(n) * areacor(n)
                   else
                      !output(i,j)  = dataPtr1d(n)
                      output(n)  = dataPtr1d(n)
                   end if
                endif
             enddo
          enddo

       else if (geomtype == ESMF_GEOMTYPE_GRID) then

          call state_getfldptr(state, trim(fldname), dataptr2d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          lbnd1 = lbound(dataPtr2d,1)
          lbnd2 = lbound(dataPtr2d,2)

          n = 0 
          do j = jsc, jec
             j1 = j + lbnd2 - jsc
             do i = isc, iec
                i1 = i + lbnd1 - isc
                n = n +1
                if (present(do_sum)) then
                   !output(i,j) = output(i,j) + dataPtr2d(i1,j1)
                   output(n) = output(n) + dataPtr2d(i1,j1)
                else
                   output(n) = dataPtr2d(i1,j1)
                   !output(i,j) = dataPtr2d(i1,j1)
                endif
             enddo
          enddo
       endif
    endif

  end subroutine State_GetImport_Int

    ! From MOM: Map import state field to output array
  subroutine State_GetImport_Log(state, fldname, isc, iec, jsc, jec, output, do_sum, areacor, rc)
    type(ESMF_State)    , intent(in)    :: state   !< ESMF state
    character(len=*)    , intent(in)    :: fldname !< Field name
    integer             , intent(in)    :: isc     !< The start i-index of cell centers within
    !! the computational domain
    integer             , intent(in)    :: iec     !< The end i-index of cell centers within the
    !! computational domain
    integer             , intent(in)    :: jsc     !< The start j-index of cell centers within
    !! the computational domain
    integer             , intent(in)    :: jec     !< The end j-index of cell centers within
    !! the computational domain
    !real (ESMF_KIND_R8) , intent(inout) :: output(isc:iec,jsc:jec)!< Output 2D array
    logical , intent(inout) :: output((iec-isc+1)*(jec-isc+1))!< Output 1D array
    logical, optional   , intent(in)    :: do_sum  !< If true, sums the data
    real (ESMF_KIND_R8), optional,  intent(in) :: areacor(:) !< flux area correction factors
    !! applicable to meshes
    integer             , intent(out)   :: rc      !< Return code

    ! local variables
    type(ESMF_StateItem_Flag)     :: itemFlag
    integer                       :: n, i, j, i1, j1
    integer                       :: lbnd1,lbnd2
    real(ESMF_KIND_R8), pointer   :: dataPtr1d(:)
    real(ESMF_KIND_R8), pointer   :: dataPtr2d(:,:)
    character(len=*)  , parameter :: subname='(Noah import_fields:state_getimport)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_StateGet(State, trim(fldname), itemFlag, rc=rc)
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then

       if (geomtype == ESMF_GEOMTYPE_MESH) then

          ! get field pointer
          call state_getfldptr(state, trim(fldname), dataptr1d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! determine output array and apply area correction if present
          n = 0
          do j = jsc,jec
             do i = isc,iec
                n = n + 1
                if (present(do_sum)) then
                   if (present(areacor)) then
                      !output(i,j)  = output(i,j) + dataPtr1d(n) * areacor(n)
                      output(n)  = output(n) + dataPtr1d(n) * areacor(n)
                   else
                      !output(i,j)  = output(i,j) + dataPtr1d(n)
                      output(n)  = output(n) + dataPtr1d(n)
                   end if
                else
                   if (present(areacor)) then
                      !output(i,j)  = dataPtr1d(n) * areacor(n)
                      output(n)  = dataPtr1d(n) * areacor(n)
                   else
                      !output(i,j)  = dataPtr1d(n)
                      output(n)  = dataPtr1d(n)
                   end if
                endif
             enddo
          enddo

       else if (geomtype == ESMF_GEOMTYPE_GRID) then

          call state_getfldptr(state, trim(fldname), dataptr2d, rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          lbnd1 = lbound(dataPtr2d,1)
          lbnd2 = lbound(dataPtr2d,2)

          n = 0 
          do j = jsc, jec
             j1 = j + lbnd2 - jsc
             do i = isc, iec
                i1 = i + lbnd1 - isc
                n = n +1
                if (present(do_sum)) then
                   !output(i,j) = output(i,j) + dataPtr2d(i1,j1)
                   output(n) = output(n) + dataPtr2d(i1,j1)
                else
                   output(n) = dataPtr2d(i1,j1)
                   !output(i,j) = dataPtr2d(i1,j1)
                endif
             enddo
          enddo
       endif
    endif

  end subroutine State_GetImport_Log


  
end module import_fields_mod
