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
     module procedure import_field_real
     module procedure import_field_int
     module procedure import_field_log
  end interface import_fieldv2
  
  public write_import_field, import_allfields, export_allfields

  character(*),parameter :: u_FILE_u = &
       __FILE__

contains

  subroutine import_field_real(State_i, fldname_i, out_data, rc)

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

    
  end subroutine import_field_real

  subroutine import_field_int(State_i, fldname_i, out_data, rc)

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

    
  end subroutine import_field_int

  subroutine import_field_log(State_i, fldname_i, out_data, rc)

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

    
  end subroutine import_field_log

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


  subroutine import_allfields(State_i, procbounds, noah_pubinst, rc)

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    use noah_type_mod, only: noah_type
    
    ! ----------------------------------------------
    ! Get pointers to internal land model variables, for multiple fields
    ! ----------------------------------------------
    type(ESMF_State),       intent(in)    :: State_i
    type(procbounds_type),  intent(in)    :: procbounds 
    type(noah_type),        intent(inout) :: noah_pubinst ! land model's variable type
    integer,                intent(out)   :: rc

    ! local
    !type(ESMF_State)        :: importState
    character(len=25)       :: fldname_i
    real(r8), pointer       :: fldPtr1d_i(:)
    integer                 :: g, gridbeg, gridend,im
    
    gridbeg = procbounds%gridbeg
    gridend = procbounds%gridend
    im      = procbounds%im

    ! tmp debug
    !write(*,*) 'IF1 debug: ', im, noah_pubinst%model%foo_atm2lndfield(1:im)
    
    call import_field_int(State_i, 'soiltyp', noah_pubinst%model%soiltyp(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    call import_fieldv2(State_i, 'foo_atm2lndfield', noah_pubinst%model%foo_atm2lndfield(1:im), rc=rc)
    !call import_fieldv2(State_i, 'foo_atm2lndfield', fldPtr1d_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


    call import_fieldv2(State_i, 'Faxa_soiltyp', noah_pubinst%model%soiltyp(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_vegtype', noah_pubinst%model%vegtype(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_sigmaf', noah_pubinst%model%sigmaf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_sfcemis', noah_pubinst%model%sfcemis(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_dlwflx', noah_pubinst%model%dlwflx(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_snet', noah_pubinst%model%snet(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_tg3', noah_pubinst%model%tg3(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_cm', noah_pubinst%model%cm(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_ch', noah_pubinst%model%ch(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_prsl1', noah_pubinst%model%prsl1(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_prslki', noah_pubinst%model%prslki(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_zf', noah_pubinst%model%zf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_land', noah_pubinst%model%land(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_slopetyp', noah_pubinst%model%slopetyp(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_shdmin', noah_pubinst%model%shdmin(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_shdmax', noah_pubinst%model%shdmax(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_snoalb', noah_pubinst%model%snoalb(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_sfalb', noah_pubinst%model%sfalb(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_bexppert', noah_pubinst%model%bexppert(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_xlaipert', noah_pubinst%model%xlaipert(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_vegfpert', noah_pubinst%model%vegfpert(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    call import_fieldv2(State_i, 'Faxa_ps', noah_pubinst%model%ps(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_t1', noah_pubinst%model%t1(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_q1', noah_pubinst%model%q1(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    
    ! call import_fieldv2(State_i, 'Faxa_prsik1', noah_pubinst%model%prsik1(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_sfalb', noah_pubinst%model%sfalb(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_weasd', noah_pubinst%model%weasd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_snwdph', noah_pubinst%model%snwdph(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_tskin', noah_pubinst%model%tskin(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_tprcp', noah_pubinst%model%tprcp(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_srflag', noah_pubinst%model%srflag(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call import_fieldv2(State_i, 'Faxa_smc', noah_pubinst%model%smc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call import_fieldv2(State_i, 'Faxa_stc', noah_pubinst%model%stc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call import_fieldv2(State_i, 'Faxa_slc', noah_pubinst%model%slc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_canopy', noah_pubinst%model%canopy(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_trans', noah_pubinst%model%trans(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call import_fieldv2(State_i, 'Faxa_tsurf', noah_pubinst%model%tsurf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call import_fieldv2(State_i, 'Faxa_z0rl', noah_pubinst%model%z0rl(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call import_fieldv2(State_i, 'Faxa_z0pert', noah_pubinst%model%z0pert(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call import_fieldv2(State_i, 'Faxa_ztpert', noah_pubinst%model%ztpert(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call import_fieldv2(State_i, 'Faxa_ustar', noah_pubinst%model%ustar(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine import_allfields
  

  subroutine export_allfields(State_o, procbounds, noah_pubinst, rc)

    use ESMF, only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF, only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF, only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    use noah_type_mod, only: noah_type

    ! ----------------------------------------------
    ! Get pointers to internal land model variables, for multiple fields
    ! ----------------------------------------------
    type(ESMF_State),       intent(in)    :: State_o
    type(procbounds_type),  intent(in)    :: procbounds
    type(noah_type),        intent(in)    :: noah_pubinst ! land model's variable type
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
    call state_setexport_1d(State_o, 'Fall_weasd', noah_pubinst%model%weasd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_snwdph', noah_pubinst%model%snwdph(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_tskin', noah_pubinst%model%tskin(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_tprcp', noah_pubinst%model%tprcp(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_srflag', noah_pubinst%model%srflag(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_setexport_1d(State_o, 'Fall_smc', noah_pubinst%model%smc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_setexport_1d(State_o, 'Fall_stc', noah_pubinst%model%stc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call state_setexport_1d(State_o, 'Fall_slc', noah_pubinst%model%slc(1:im), rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_canopy', noah_pubinst%model%canopy(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_trans', noah_pubinst%model%trans(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_tsurf', noah_pubinst%model%tsurf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_z0rl', noah_pubinst%model%zorl(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! noahouts
    call state_setexport_1d(State_o, 'Fall_sncovr1', noah_pubinst%model%sncovr1(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_qsurf', noah_pubinst%model%qsurf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_gflux', noah_pubinst%model%gflux(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_drain', noah_pubinst%model%drain(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_evap', noah_pubinst%model%evap(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_hflx', noah_pubinst%model%hflx(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_ep', noah_pubinst%model%ep(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_runoff', noah_pubinst%model%runoff(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_cmm', noah_pubinst%model%cmm(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_chh', noah_pubinst%model%chh(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_evbs', noah_pubinst%model%evbs(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_evcw', noah_pubinst%model%evcw(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_sbsno', noah_pubinst%model%sbsno(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_snowc', noah_pubinst%model%snowc(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_stm', noah_pubinst%model%stm(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_snohf', noah_pubinst%model%snohf(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_smcwlt2', noah_pubinst%model%smcwlt2(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_smcref2', noah_pubinst%model%smcref2(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_wet1', noah_pubinst%model%wet1(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! diffouts
    call state_setexport_1d(State_o, 'Fall_rb_lnd', noah_pubinst%model%rb_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_fm_lnd', noah_pubinst%model%fm_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_fh_lnd', noah_pubinst%model%fh_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_fm10_lnd', noah_pubinst%model%fm10_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_fh2_lnd', noah_pubinst%model%fh2_lnd(1:im), rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call state_setexport_1d(State_o, 'Fall_stress', noah_pubinst%model%stress(1:im), rc=rc)
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
!  subroutine state_getfldptr(State, fldname, fldptr1d, fldptr2d, rc)
  subroutine state_getfldptr(State, fldname, fldptr1d, rc)

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
    !real(r8), pointer, optional , intent(out)   :: fldptr2d(:,:)
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
       !if (ChkErr(rc,__LINE__,u_FILE_u)) return
    ! else if (present(fldptr2d)) then
    !    call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
    !    !if (ChkErr(rc,__LINE__,u_FILE_u)) return
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
