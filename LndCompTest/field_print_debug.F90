module field_print_debug

  use ESMF                    , only : ESMF_GridComp, ESMF_State, ESMF_Mesh, ESMF_StateGet
  use ESMF                    , only : ESMF_KIND_R8, ESMF_SUCCESS, ESMF_MAXSTR, ESMF_LOGMSG_INFO
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use nuopc_shr_methods       , only : chkerr
  implicit none

  private

  public field_foo

  character(*),parameter :: u_FILE_u = &
       __FILE__

contains

  subroutine field_foo(State, fldname, fldptr2d, fldptr1d, rc)

    ! ----------------------------------------------
    ! Get pointer to a state field
    ! ----------------------------------------------

    use ESMF , only : ESMF_State, ESMF_Field, ESMF_Mesh, ESMF_FieldStatus_Flag
    use ESMF , only : ESMF_StateGet, ESMF_FieldGet, ESMF_MeshGet
    use ESMF , only : ESMF_FIELDSTATUS_COMPLETE, ESMF_FAILURE

    ! input/output variables
    type(ESMF_State),             intent(in)    :: State
    character(len=*),             intent(in)    :: fldname
    real(r8), pointer, optional , intent(out)   :: fldptr1d(:)
    real(r8), pointer, optional , intent(out)   :: fldptr2d(:,:)
    integer,                      intent(out)   :: rc

    ! local variables
    type(ESMF_FieldStatus_Flag) :: status
    type(ESMF_Field)            :: lfield
    integer                     :: i,j

    write(*,*) 'JP foo1'
    call ESMF_StateGet(State, itemName=trim(fldname), field=lfield, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    write(*,*) 'JP foo2'

    if (present(fldptr1d)) then
       write(*,*) 'JP foo3'
       call ESMF_FieldGet(lfield, farrayPtr=fldptr1d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(*,*) 'JP foo3b'
       do j=1,size(fldptr1d)
          write(*,*) 'JP j ', j, fldptr1d(j)
       enddo
       
    else if (present(fldptr2d)) then
       write(*,*) 'JP foo4'
       call ESMF_FieldGet(lfield, farrayPtr=fldptr2d, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       write(*,*) 'JP foo4b'

       do j=lbound(fldptr2d,2),ubound(fldptr2d,2)
          do i=lbound(fldptr2d,1),ubound(fldptr2d,1)
             write(*,*) 'JP i,j ', i,j, fldptr2d(i,j)
          enddo
       enddo

    else
       call shr_sys_abort("JP: either fldptr1d or fldptr2d must be an input argument")
    end if



  end subroutine field_foo

end module field_print_debug
