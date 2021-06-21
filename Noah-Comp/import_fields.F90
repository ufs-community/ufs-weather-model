module import_fields

  use ESMF,                   only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF,                   only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use ESMF,                   only : ESMF_STATEITEM_NOTFOUND, ESMF_StateItem_Flag
  use ESMF,                   only : operator(/=), operator(==)
  use ESMF,                   only : ESMF_FieldWrite
  use shr_kind_mod,           only : r8 => shr_kind_r8
  use nuopc_shr_methods,      only : chkerr

  use proc_bounds,            only: procbounds_type
  
  implicit none

  private

  interface import_fieldv2
     module procedure import_field_1Dreal
     module procedure import_field_1Dint
     module procedure import_field_1Dlog
     module procedure import_field_2Dreal
     module procedure import_field_2Dint
     module procedure import_field_2Dlog
  end interface import_fieldv2
  
  public write_import_field, import_allfields, export_allfields

  character(*),parameter :: u_FILE_u = &
       __FILE__

contains

  subroutine import_field_1Dreal(State_i, fldname_i, out_data, rc)

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),             intent(in)    :: State_i
    character(len=*),             intent(in)    :: fldname_i
    real(r8),                     intent(inout) :: out_data(:)
    integer,                      intent(out)   :: rc

    ! local
    real(r8), pointer           :: fldPtr1d_i(:)
    integer                     :: i

    
    if (fldchk(State_i,trim(fldname_i))) then
       call state_getfldptr(State_i,trim(fldname_i), fldptr1d_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
       do i=1,size(fldptr1d_i)
          out_data(i) = fldptr1d_i(i)
       enddo

    else
       write(*,*) 'Field ', trim(fldname_i), ' NOT found'

    endif

    
  end subroutine import_field_1Dreal

  subroutine import_field_2Dreal(State_i, fldname_i, out_data, rc)

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),             intent(in)    :: State_i
    character(len=*),             intent(in)    :: fldname_i
    real(r8),                     intent(inout) :: out_data(:,:)
    integer,                      intent(out)   :: rc

    ! local
    real(r8), pointer           :: fldPtr2d_i(:,:)
    integer                     :: i, j

    
    if (fldchk(State_i,trim(fldname_i))) then
       call state_getfldptr(State_i,trim(fldname_i), fldptr2d=fldptr2d_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
       do i=1,size(fldptr2d_i, dim=1)
          do j=1,size(fldptr2d_i, dim=2)
             out_data(i,j) = fldptr2d_i(i,j)
          enddo
       enddo

    else
       write(*,*) 'Field ', trim(fldname_i), ' NOT found'

    endif

    
  end subroutine import_field_2Dreal

  subroutine import_field_1Dint(State_i, fldname_i, out_data, rc)

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),            intent(in)    :: State_i
    character(len=*),            intent(in)    :: fldname_i
    integer,                     intent(inout) :: out_data(:)
    integer,                     intent(out)   :: rc

    ! local
    real(r8), pointer           :: fldPtr1d_i(:)
    integer                     :: i

    
    if (fldchk(State_i,trim(fldname_i))) then
       call state_getfldptr(State_i,trim(fldname_i), fldptr1d_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
       do i=1,size(fldptr1d_i)
          out_data(i) = fldptr1d_i(i)
       enddo

    else
       write(*,*) 'Field ', trim(fldname_i), ' NOT found'

    endif

    
  end subroutine import_field_1Dint

  subroutine import_field_2Dint(State_i, fldname_i, out_data, rc)

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),            intent(in)    :: State_i
    character(len=*),            intent(in)    :: fldname_i
    integer,                     intent(inout) :: out_data(:,:)
    integer,                     intent(out)   :: rc

    ! local
    real(r8), pointer           :: fldPtr2d_i(:,:)
    integer                     :: i,j

    
    if (fldchk(State_i,trim(fldname_i))) then
       call state_getfldptr(State_i,trim(fldname_i), fldptr2d=fldptr2d_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
       do i=1,size(fldptr2d_i, dim=1)
          do j=1,size(fldptr2d_i, dim=2)
             out_data(i,j) = fldptr2d_i(i,j)
          enddo
       enddo

    else
       write(*,*) 'Field ', trim(fldname_i), ' NOT found'

    endif

    
  end subroutine import_field_2Dint
  
  subroutine import_field_1Dlog(State_i, fldname_i, out_data, rc)

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),            intent(in)    :: State_i
    character(len=*),            intent(in)    :: fldname_i
    logical,                     intent(inout) :: out_data(:)
    integer,                     intent(out)   :: rc

    ! local
    real(r8), pointer           :: fldPtr1d_i(:)
    integer                     :: i

    
    if (fldchk(State_i,trim(fldname_i))) then
       call state_getfldptr(State_i,trim(fldname_i), fldptr1d_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
       do i=1,size(fldptr1d_i)
          out_data(i) = fldptr1d_i(i)
       enddo

    else
       write(*,*) 'Field ', trim(fldname_i), ' NOT found'

    endif

    
  end subroutine import_field_1Dlog

  subroutine import_field_2Dlog(State_i, fldname_i, out_data, rc)

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    type(ESMF_State),            intent(in)    :: State_i
    character(len=*),            intent(in)    :: fldname_i
    logical,                     intent(inout) :: out_data(:,:)
    integer,                     intent(out)   :: rc

    ! local
    real(r8), pointer           :: fldPtr2d_i(:,:)
    integer                     :: i,j

    
    if (fldchk(State_i,trim(fldname_i))) then
       call state_getfldptr(State_i,trim(fldname_i), fldptr2d=fldptr2d_i, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       
       do i=1,size(fldptr2d_i, dim=1)
          do j=1,size(fldptr2d_i, dim=2)
             out_data(i,j) = fldptr2d_i(i,j)
          enddo
       enddo

    else
       write(*,*) 'Field ', trim(fldname_i), ' NOT found'

    endif

    
  end subroutine import_field_2Dlog
  

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


  subroutine import_allfields(State_i, procbounds, noah_model, ctrl_init, rc)

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
    integer                 :: g, im
    

    if (trim(ctrl_init%grid) == '1Dmesh') then
       
       im      = procbounds%im

       ! tmp for testing
       call import_fieldv2(State_i, 'foo_atm2lndfield', noah_model%model%foo_atm2lndfield(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call import_fieldv2(State_i, 'Faxa_soiltyp', noah_model%model%soiltyp(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_vegtype', noah_model%model%vegtype(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_sigmaf', noah_model%model%sigmaf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_sfcemis', noah_model%model%sfcemis(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_dlwflx', noah_model%model%dlwflx(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_dswsfc', noah_model%model%dswsfc(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return    
       call import_fieldv2(State_i, 'inst_down_sw_flx', noah_model%model%dswsfci(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return    
       call import_fieldv2(State_i, 'Faxa_snet', noah_model%model%snet(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_tg3', noah_model%model%tg3(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_cm', noah_model%model%cm(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_ch', noah_model%model%ch(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_prsl1', noah_model%model%prsl1(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_prslki', noah_model%model%prslki(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_zf', noah_model%model%zf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_land', noah_model%model%land(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_slopetyp', noah_model%model%slopetyp(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_shdmin', noah_model%model%shdmin(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_shdmax', noah_model%model%shdmax(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_snoalb', noah_model%model%snoalb(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_sfalb', noah_model%model%sfalb(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_bexppert', noah_model%model%bexppert(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_xlaipert', noah_model%model%xlaipert(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_vegfpert', noah_model%model%vegfpert(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call import_fieldv2(State_i, 'Faxa_ps', noah_model%model%ps(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_t1', noah_model%model%t1(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_q1', noah_model%model%q1(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! call import_fieldv2(State_i, 'Faxa_prsik1', noah_model%model%prsik1(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_sfalb', noah_model%model%sfalb(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_weasd', noah_model%model%weasd(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_snwdph', noah_model%model%snwdph(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_tskin', noah_model%model%tskin(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_tprcp', noah_model%model%tprcp(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_srflag', noah_model%model%srflag(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_smc', noah_model%model%smc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_stc', noah_model%model%stc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_slc', noah_model%model%slc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_canopy', noah_model%model%canopy(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_trans', noah_model%model%trans(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_tsurf', noah_model%model%tsurf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_z0rl', noah_model%model%z0rl(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_z0pert', noah_model%model%z0pert(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_ztpert', noah_model%model%ztpert(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_ustar', noah_model%model%ustar(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_wind', noah_model%model%wind(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return    

    else if(trim(ctrl_init%grid) == 'CS') then ! blocking


       ! tmp for testing
       call import_fieldv2(State_i, 'foo_atm2lndfield', noah_model%model%foo_atm2lndfield(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call import_fieldv2(State_i, 'Faxa_soiltyp', noah_model%model%soiltyp(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_vegtype', noah_model%model%vegtype(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_sigmaf', noah_model%model%sigmaf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_sfcemis', noah_model%model%sfcemis(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_dlwflx', noah_model%model%dlwflx(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_dswsfc', noah_model%model%dswsfc(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return    
       call import_fieldv2(State_i, 'inst_down_sw_flx', noah_model%model%dswsfci(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return    
       call import_fieldv2(State_i, 'Faxa_snet', noah_model%model%snet(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_tg3', noah_model%model%tg3(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_cm', noah_model%model%cm(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_ch', noah_model%model%ch(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_prsl1', noah_model%model%prsl1(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_prslki', noah_model%model%prslki(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_zf', noah_model%model%zf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_land', noah_model%model%land(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_slopetyp', noah_model%model%slopetyp(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_shdmin', noah_model%model%shdmin(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_shdmax', noah_model%model%shdmax(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_snoalb', noah_model%model%snoalb(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_sfalb', noah_model%model%sfalb(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_bexppert', noah_model%model%bexppert(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_xlaipert', noah_model%model%xlaipert(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_vegfpert', noah_model%model%vegfpert(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call import_fieldv2(State_i, 'Faxa_ps', noah_model%model%ps(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_t1', noah_model%model%t1(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_q1', noah_model%model%q1(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! call import_fieldv2(State_i, 'Faxa_prsik1', noah_model%model%prsik1(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_sfalb', noah_model%model%sfalb(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_weasd', noah_model%model%weasd(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_snwdph', noah_model%model%snwdph(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_tskin', noah_model%model%tskin(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_tprcp', noah_model%model%tprcp(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_srflag', noah_model%model%srflag(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_smc', noah_model%model%smc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_stc', noah_model%model%stc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_slc', noah_model%model%slc(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_canopy', noah_model%model%canopy(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_trans', noah_model%model%trans(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_tsurf', noah_model%model%tsurf(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_z0rl', noah_model%model%z0rl(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_z0pert', noah_model%model%z0pert(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       ! call import_fieldv2(State_i, 'Faxa_ztpert', noah_model%model%ztpert(1:im), rc=rc)
       ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_ustar', noah_model%model%ustar(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call import_fieldv2(State_i, 'Faxa_wind', noah_model%model%wind(1:im), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return    

    else
       write(*,*) 'Do not know if imports/exports should be 1D or 2D from ctrl_init%grid'
       stop
       
    endif
    
  end subroutine import_allfields
  

  subroutine export_allfields(State_o, procbounds, noah_model, rc)

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    use noah_type_mod, only: noah_type

    ! ----------------------------------------------
    ! Get pointers to internal land model variables, for multiple fields
    ! ----------------------------------------------
    type(ESMF_State),       intent(in)    :: State_o
    type(procbounds_type),  intent(in)    :: procbounds
    type(noah_type),        intent(in)    :: noah_model ! land model's variable type
    integer,                intent(out)   :: rc

    ! local
    !type(ESMF_State)        :: importState
    character(len=25)       :: fldname_i
    real(r8), pointer       :: fldPtr1d_i(:)
    integer                 :: g, gridbeg, gridend,im

    gridbeg = procbounds%gridbeg
    gridend = procbounds%gridend
    im      = procbounds%im

    ! inouts   
    call state_setexport_1d(State_o, 'Fall_weasd', noah_model%model%weasd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_snwdph', noah_model%model%snwdph(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_tskin', noah_model%model%tskin(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_tprcp', noah_model%model%tprcp(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_srflag', noah_model%model%srflag(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_setexport_1d(State_o, 'Fall_smc', noah_model%model%smc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_setexport_1d(State_o, 'Fall_stc', noah_model%model%stc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_setexport_1d(State_o, 'Fall_slc', noah_model%model%slc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_canopy', noah_model%model%canopy(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_trans', noah_model%model%trans(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_tsurf', noah_model%model%tsurf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_z0rl', noah_model%model%z0rl(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! noahouts
    call state_setexport_1d(State_o, 'Fall_sncovr1', noah_model%model%sncovr1(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_qsurf', noah_model%model%qsurf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_gflux', noah_model%model%gflux(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_drain', noah_model%model%drain(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_evap', noah_model%model%evap(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_hflx', noah_model%model%hflx(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_ep', noah_model%model%ep(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_runoff', noah_model%model%runoff(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_cmm', noah_model%model%cmm(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_chh', noah_model%model%chh(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_evbs', noah_model%model%evbs(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_evcw', noah_model%model%evcw(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_sbsno', noah_model%model%sbsno(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_snowc', noah_model%model%snowc(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_stm', noah_model%model%stm(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_snohf', noah_model%model%snohf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_smcwlt2', noah_model%model%smcwlt2(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_smcref2', noah_model%model%smcref2(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_wet1', noah_model%model%wet1(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! diffouts
    call state_setexport_1d(State_o, 'Fall_rb_lnd', noah_model%model%rb_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_fm_lnd', noah_model%model%fm_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_fh_lnd', noah_model%model%fh_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_fm10_lnd', noah_model%model%fm10_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_fh2_lnd', noah_model%model%fh2_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_stress', noah_model%model%stress(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
  end subroutine export_allfields

  
!===============================================================================

  subroutine state_setexport_1d(state, fldname, outdata, rc)

    ! fill in ctsm export data for 1d field

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

    call state_getfldptr(state, trim(fldname), fldptr1d=fldptr1d, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    fldptr1d(:) = 0._r8
    do g = 1,size(outdata)
       fldptr1d(g) = outdata(g)
    end do

  end subroutine state_setexport_1d
  
!===============================================================================
  subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    use shr_sys_mod, only : shr_sys_abort
    
    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fldname
    real(r8), pointer, optional , intent(out)   :: fldptr1d(:)
    real(r8), pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    character(len=*), parameter :: subname='(lnd_import_export:state_getfldptr)'
    type(ESMF_StateItem_Flag)   :: itemFlag
    ! ----------------------------------------------

    rc = ESMF_SUCCESS
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (present(fldptr1d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else if (present(fldptr2d)) then
       call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call shr_sys_abort("either fldptr1d or fldptr2d must be an input argument")
    end if
  end subroutine state_getfldptr

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


end module import_fields
