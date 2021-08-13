module land_restart_mod

  !-----------------------------------------------------------------------
  ! Generally templating from FV3GFS_io_mod
  ! For now, without blocking
  !-----------------------------------------------------------------------

  use machine,            only: kind_phys
  use mpp_mod,            only: mpp_error,  mpp_pe, mpp_root_pe, &
                                mpp_chksum, mpp_sync, NOTE,   FATAL

  use fms_mod,            only: file_exist, stdout
  use fms_io_mod,         only: restart_file_type, free_restart_type, &
                                register_restart_field, restore_state, save_restart
  use mpp_domains_mod,    only: domain1d, domain2d, domainUG

  use noah_type_mod,      only: noah_type

  !-----------------------------------------------------------------------
  implicit none
  private


  public sfc_prop_restart_read, sfc_prop_restart_write


  !--- GFDL filenames
  character(len=32)  :: fn_oro = 'oro_data.nc'

  character(len=32)  :: fn_srf = 'sfc_data.nc'
  character(len=32)  :: fn_srf_TEST = 'sfc_data_TEST.nc'

  !--- GFDL FMS netcdf restart data types
  type(restart_file_type) :: Oro_restart, Sfc_restart

  !--- GFDL FMS restart containers
  character(len=32),    allocatable,         dimension(:)       :: oro_name2, sfc_name2
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: oro_var2, sfc_var2

contains
  !-----------------------------------------------------------------------

  subroutine sfc_prop_restart_read(noah_model, land_domain, warm_start)

    type (noah_type),          intent(inout) :: noah_model
    type (domain2d),           intent(in)    :: land_domain
    logical,                   intent(in)    :: warm_start

    !--- local variables
    integer :: i, j, k, ix, lsoil, num, j1, i1
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar_o2
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()

    nvar_o2  = 1

    isc = noah_model%control%isc
    iec = noah_model%control%iec
    jsc = noah_model%control%jsc
    jec = noah_model%control%jec
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)

    !--- OROGRAPHY FILE
    if (.not. allocated(oro_name2)) then
       !--- allocate the various containers needed for orography data
       allocate(oro_name2(nvar_o2))
       allocate(oro_var2(nx,ny,nvar_o2))
       oro_var2 = -9999._kind_phys

       oro_name2(1) = 'land_frac'  ! land fraction [0:1]
       do num = 1,nvar_o2
          var2_p => oro_var2(:,:,num)
          if (trim(oro_name2(num)) == 'lake_frac' .or. trim(oro_name2(num)) == 'lake_depth') then
             id_restart = register_restart_field(Oro_restart, fn_oro, oro_name2(num), var2_p, domain=land_domain, mandatory=.false.)
          else
             id_restart = register_restart_field(Oro_restart, fn_oro, oro_name2(num), var2_p, domain=land_domain)
          endif
       enddo
       nullify(var2_p)
    endif

    !--- read the orography restart/data
    call mpp_error(NOTE,'reading topographic/orographic information from INPUT/oro_data.tile*.nc')
    call restore_state(Oro_restart)

    call mpp_sync() !test sync to debug
    
    !--- copy data into GFS containers

    ix = 0
    write(*,*) 'In rIO: ', jsc,jec,isc,iec
    do j = jsc, jec
       j1 = j - jsc + 1
       do i = isc, iec
          i1 = i - isc +1
          ix = ix + 1

          noah_model%sfcprop%landfrac(ix)  = oro_var2(i1,j1,1) !land frac [0:1]

       enddo
    enddo

    write(*,*) 'Restart read test: ', oro_var2(i1,j1,1)
    
    !--- deallocate containers and free restart container
    deallocate(oro_name2, oro_var2)
    call free_restart_type(Oro_restart)


    
  end subroutine sfc_prop_restart_read

  !----------------------------------------------------------------------

  subroutine sfc_prop_restart_write (noah_model, land_domain, timestamp)
    !--- interface variable definitions
    type (noah_type),            intent(inout) :: noah_model
    type (domain2d),             intent(in)    :: land_domain
    character(len=32), optional, intent(in)    :: timestamp


    !--- local variables
    integer :: i, j, k, ix, lsoil, num, j1, i1
    integer :: isc, iec, jsc, jec, npz, nx, ny
    integer :: id_restart
    integer :: nvar_o2

    integer :: nvar2m, nvar2o, nvar3
    integer :: nvar2r, nvar2mp, nvar3mp

    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()

    nvar2m  = 1
    ! copied over from FV3GFS_io, can clean up unused
    nvar2o  = 0
    nvar3   = 0
    nvar2r  = 0
    nvar2mp = 0
    nvar3mp = 0

    isc = noah_model%control%isc
    iec = noah_model%control%iec
    jsc = noah_model%control%jsc
    jec = noah_model%control%jec
    nx  = (iec - isc + 1)
    ny  = (jec - jsc + 1)


    if (.not. allocated(sfc_name2)) then
       !--- allocate the various containers needed for restarts
       allocate(sfc_name2(nvar2m+nvar2o+nvar2mp+nvar2r))
       !allocate(sfc_name3(0:nvar3+nvar3mp))
       allocate(sfc_var2(nx,ny,nvar2m+nvar2o+nvar2mp+nvar2r))
       ! if (Model%lsm == Model%lsm_noah .or. Model%lsm == Model%lsm_noahmp .or. Model%lsm == Model%lsm_noah_wrfv4) then
       !    allocate(sfc_var3(nx,ny,Model%lsoil,nvar3))
       ! elseif (Model%lsm == Model%lsm_ruc) then
       !    allocate(sfc_var3(nx,ny,Model%lsoil_lsm,nvar3))
       ! endif
       sfc_var2   = -9999.0_kind_phys
       !sfc_var3   = -9999.0_kind_phys
       ! if (Model%lsm == Model%lsm_noahmp) then
       !    allocate(sfc_var3sn(nx,ny,-2:0,4:6))
       !    allocate(sfc_var3eq(nx,ny,1:4,7:7))
       !    allocate(sfc_var3zn(nx,ny,-2:4,8:8))

       !    sfc_var3sn = -9999.0_kind_phys
       !    sfc_var3eq = -9999.0_kind_phys
       !    sfc_var3zn = -9999.0_kind_phys
       ! endif

       !--- names of the 2D variables to save
       sfc_name2(1)  = 'landfrac'
       !sfc_name2(2)  = 'tsea'    !tsfc

       !--- register the 2D fields
       do num = 1,nvar2m
          var2_p => sfc_var2(:,:,num)
             id_restart = register_restart_field(Sfc_restart, fn_srf_TEST, sfc_name2(num), var2_p, domain=land_domain)
       enddo

       nullify(var2_p)
    end if

       ix = 0
       do j = jsc, jec
          j1 = j - jsc + 1
          do i = isc, iec
             i1 = i - isc +1
             ix = ix + 1

             !noah_model%sfcprop%landfrac(ix)  = oro_var2(i1,j1,1) !land frac [0:1]

             sfc_var2(i1,j1,1) = noah_model%sfcprop%landfrac(ix) !--- slmsk
             !sfc_var2(i1,j1,2) = noah_model%sfcprop%ltsfco(ix) !--- tsfc (tsea in sfc file)
          enddo
       enddo


       call mpp_sync() !test sync to debug 
       call save_restart(Sfc_restart, timestamp)
       call mpp_sync() !test sync to debug 
       
     end subroutine sfc_prop_restart_write


   end module land_restart_mod
