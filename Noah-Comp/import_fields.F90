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
  end interface import_fieldv2
  
  public write_import_field, import_allfields

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
    
    !call import_field_int(State_i, fldname_i, noah_pubinst%model%soiltyp(1:im), rc=rc)
    ! call import_field_int(State_i, fldname_i, fldPtr1d_i, rc=rc)
    ! if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call import_fieldv2(State_i, 'foo_atm2lndfield', noah_pubinst%model%foo_atm2lndfield(1:im), rc=rc)
    !call import_fieldv2(State_i, 'foo_atm2lndfield', fldPtr1d_i, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return


  end subroutine import_allfields
  
 
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
