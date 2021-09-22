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


  public sfc_prop_restart_read, sfc_prop_restart_write, sfc_prop_transfer


  !--- GFDL filenames
  character(len=32)  :: fn_oro = 'oro_data.nc'

  character(len=32)  :: fn_srf = 'sfc_data.nc'
  character(len=32)  :: fn_srf_TEST = 'sfc_data_TEST.nc'

  !--- GFDL FMS netcdf restart data types
  type(restart_file_type) :: Oro_restart, Sfc_restart

  !--- GFDL FMS restart containers
  character(len=32),    allocatable,         dimension(:)       :: oro_name2, sfc_name2, sfc_name3
  real(kind=kind_phys), allocatable, target, dimension(:,:,:)   :: oro_var2, sfc_var2
  real(kind=kind_phys), allocatable, target, dimension(:,:,:,:) :: sfc_var3

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
    integer :: nvar_o2, nvar_s2, nvar_s3
    real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
    real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p  => NULL()
    
    nvar_o2  = 1   ! 2d oro fields
    nvar_s2  = 32  ! 2d surface data fields
    nvar_s3  = 3   ! 3d surface data fields

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
    do j = jsc, jec
       j1 = j - jsc + 1
       do i = isc, iec
          i1 = i - isc +1
          ix = ix + 1

          noah_model%sfcprop%landfrac(ix)  = oro_var2(i1,j1,1) !land frac [0:1]
       enddo
    enddo

    !--- deallocate containers and free restart container
    deallocate(oro_name2, oro_var2)
    call free_restart_type(Oro_restart)


    !--- SURFACE DATA FILE
    if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for sfc data
      allocate(sfc_name2(nvar_s2))
      allocate(sfc_name3(nvar_s3))
      allocate(sfc_var2(nx,ny,nvar_s2))
      allocate(sfc_var3(nx,ny,noah_model%static%km,nvar_s3))
      sfc_var2 = -9999._kind_phys
      sfc_var3 = -9999._kind_phys
      
      sfc_name2(1)  = 'slmsk'
      sfc_name2(2)  = 'tsfcl'   
      sfc_name2(3)  = 'sheleg'  !weasd
      sfc_name2(4)  = 'tg3'
      sfc_name2(5)  = 'zorll'
      sfc_name2(6)  = 'alvsf'
      sfc_name2(7)  = 'alvwf'
      sfc_name2(8)  = 'alnsf'
      sfc_name2(9)  = 'alnwf'
      sfc_name2(10) = 'facsf'
      sfc_name2(11) = 'facwf'
      sfc_name2(12) = 'vfrac'
      sfc_name2(13) = 'canopy'
      sfc_name2(14) = 'f10m'
      sfc_name2(15) = 't2m'
      sfc_name2(16) = 'q2m'
      sfc_name2(17) = 'vtype'
      sfc_name2(18) = 'stype'
      sfc_name2(19) = 'uustar'
      sfc_name2(20) = 'ffmm'
      sfc_name2(21) = 'ffhh'
      sfc_name2(22) = 'hice'
      sfc_name2(23) = 'fice'
      sfc_name2(24) = 'tisfc'
      sfc_name2(25) = 'tprcp'
      sfc_name2(26) = 'srflag'
      sfc_name2(27) = 'snwdph'  !snowd
      sfc_name2(28) = 'shdmin'
      sfc_name2(29) = 'shdmax'
      sfc_name2(30) = 'slope'
      sfc_name2(31) = 'snoalb'
      !--- variables below here are optional
      sfc_name2(32) = 'sncovr'


      !--- register the 2D fields
      do num = 1,nvar_s2
        var2_p => sfc_var2(:,:,num)
        if (trim(sfc_name2(num)) == 'sncovr'.or. trim(sfc_name2(num)) == 'tsfcl' .or. trim(sfc_name2(num)) == 'zorll' &
                                            .or. trim(sfc_name2(num)) == 'zorli' .or. trim(sfc_name2(num)) == 'zorlwav') then
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=land_domain, mandatory=.false.)
        else
          id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name2(num), var2_p, domain=land_domain)
        endif
      enddo

      nullify(var2_p)

      !--- names of the 3D variables to save
      sfc_name3(1) = 'stc'
      sfc_name3(2) = 'smc'
      sfc_name3(3) = 'slc'
      
      do num = 1,nvar_s3
         var3_p => sfc_var3(:,:,:,num)
         id_restart = register_restart_field(Sfc_restart, fn_srf, sfc_name3(num), var3_p, domain=land_domain)
      enddo

      nullify(var3_p)
      
   endif

   !--- read the surface restart/data
   call mpp_error(NOTE,'Land reading surface properties data from INPUT/sfc_data.tile*.nc')
   call restore_state(Sfc_restart)
  
   !--- copy data into GFS containers
   ix = 0
   do j = jsc, jec
      j1 = j - jsc + 1
      do i = isc, iec
         i1 = i - isc +1
         ix = ix + 1

         noah_model%sfcprop%slmsk(ix)  = sfc_var2(i1,j1,1)    !--- slmsk
         noah_model%sfcprop%tsfcl(ix)  = sfc_var2(i1,j1,2)    !--- tsfcl
         noah_model%sfcprop%weasd(ix)  = sfc_var2(i1,j1,3)    !--- weasd (sheleg in sfc file)
         noah_model%sfcprop%tg3(ix)    = sfc_var2(i1,j1,4)    !--- tg3
         noah_model%sfcprop%zorll(ix)  = sfc_var2(i1,j1,5)    !--- zorl on land
         noah_model%sfcprop%alvsf(ix)  = sfc_var2(i1,j1,6)    !--- alvsf
         noah_model%sfcprop%alvwf(ix)  = sfc_var2(i1,j1,7)    !--- alvwf
         noah_model%sfcprop%alnsf(ix)  = sfc_var2(i1,j1,8)    !--- alnsf
         noah_model%sfcprop%alnwf(ix)  = sfc_var2(i1,j1,9)    !--- alnwf
         noah_model%sfcprop%facsf(ix)  = sfc_var2(i1,j1,10)   !--- facsf
         noah_model%sfcprop%facwf(ix)  = sfc_var2(i1,j1,11)   !--- facwf
         noah_model%sfcprop%vfrac(ix)  = sfc_var2(i1,j1,12)   !--- vfrac
         noah_model%sfcprop%canopy(ix) = sfc_var2(i1,j1,13)   !--- canopy
         noah_model%sfcprop%f10m(ix)   = sfc_var2(i1,j1,14)   !--- f10m
         noah_model%sfcprop%t2m(ix)    = sfc_var2(i1,j1,15)   !--- t2m
         noah_model%sfcprop%q2m(ix)    = sfc_var2(i1,j1,16)   !--- q2m
         noah_model%sfcprop%vtype(ix)  = sfc_var2(i1,j1,17)   !--- vtype
         noah_model%sfcprop%stype(ix)  = sfc_var2(i1,j1,18)   !--- stype
         noah_model%sfcprop%uustar(ix) = sfc_var2(i1,j1,19)   !--- uustar
         noah_model%sfcprop%ffmm(ix)   = sfc_var2(i1,j1,20)   !--- ffmm
         noah_model%sfcprop%ffhh(ix)   = sfc_var2(i1,j1,21)   !--- ffhh
         noah_model%sfcprop%hice(ix)   = sfc_var2(i1,j1,22)   !--- hice
         noah_model%sfcprop%fice(ix)   = sfc_var2(i1,j1,23)   !--- fice
         noah_model%sfcprop%tisfc(ix)  = sfc_var2(i1,j1,24)   !--- tisfc
         noah_model%sfcprop%tprcp(ix)  = sfc_var2(i1,j1,25)   !--- tprcp
         noah_model%sfcprop%srflag(ix) = sfc_var2(i1,j1,26)   !--- srflag
         noah_model%sfcprop%snowd(ix)  = sfc_var2(i1,j1,27)   !--- snowd (snwdph in the file)
         noah_model%sfcprop%shdmin(ix) = sfc_var2(i1,j1,28)   !--- shdmin
         noah_model%sfcprop%shdmax(ix) = sfc_var2(i1,j1,29)   !--- shdmax
         noah_model%sfcprop%slope(ix)  = sfc_var2(i1,j1,30)   !--- slope
         noah_model%sfcprop%snoalb(ix) = sfc_var2(i1,j1,31)   !--- snoalb
         noah_model%sfcprop%sncovr(ix) = sfc_var2(i1,j1,32)   !--- sncovr

         do lsoil = 1,noah_model%static%km
            noah_model%sfcprop%stc(ix,lsoil) = sfc_var3(i1,j1,lsoil,1)   !--- stc
            noah_model%sfcprop%smc(ix,lsoil) = sfc_var3(i1,j1,lsoil,2)   !--- smc
            noah_model%sfcprop%slc(ix,lsoil) = sfc_var3(i1,j1,lsoil,3)   !--- slc
         enddo

      enddo
   enddo   

   !--- deallocate containers and free restart container
   deallocate(sfc_name2, sfc_name3, sfc_var2)
   call free_restart_type(Sfc_restart)

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

   integer :: nvar2o, nvar2s, nvar3

   real(kind=kind_phys), pointer, dimension(:,:)   :: var2_p  => NULL()
   real(kind=kind_phys), pointer, dimension(:,:,:) :: var3_p => NULL()
   
   nvar2o  = 1  ! 2D from oro
   nvar2s  = 32 ! 2D from sfc data
   nvar3   = 3  ! 3D
   
   ! copied over from FV3GFS_io, can clean up unused
   ! nvar2r  = 0
   ! nvar2mp = 0
   ! nvar3mp = 0

   isc = noah_model%control%isc
   iec = noah_model%control%iec
   jsc = noah_model%control%jsc
   jec = noah_model%control%jec
   nx  = (iec - isc + 1)
   ny  = (jec - jsc + 1)


   if (.not. allocated(sfc_name2)) then
      !--- allocate the various containers needed for restarts
      allocate(sfc_name2(nvar2o+nvar2s))
      allocate(sfc_name3(nvar3))
      allocate(sfc_var2(nx,ny,nvar2o+nvar2s))
      sfc_var2   = -9999.0_kind_phys
      sfc_var3   = -9999.0_kind_phys

      !--- names of the 2D variables to save
      sfc_name2(1)         = 'landfrac'
      sfc_name2(nvar2o+1)  = 'slmsk'
      sfc_name2(nvar2o+2)  = 'tsfcl'
      sfc_name2(nvar2o+3)  = 'sheleg'  !weasd
      sfc_name2(nvar2o+4)  = 'tg3'
      sfc_name2(nvar2o+5)  = 'zorll'
      sfc_name2(nvar2o+6)  = 'alvsf'
      sfc_name2(nvar2o+7)  = 'alvwf'
      sfc_name2(nvar2o+8)  = 'alnsf'
      sfc_name2(nvar2o+9)  = 'alnwf'
      sfc_name2(nvar2o+10) = 'facsf'
      sfc_name2(nvar2o+11) = 'facwf'
      sfc_name2(nvar2o+12) = 'vfrac'
      sfc_name2(nvar2o+13) = 'canopy'
      sfc_name2(nvar2o+14) = 'f10m'
      sfc_name2(nvar2o+15) = 't2m'
      sfc_name2(nvar2o+16) = 'q2m'
      sfc_name2(nvar2o+17) = 'vtype'
      sfc_name2(nvar2o+18) = 'stype'
      sfc_name2(nvar2o+19) = 'uustar'
      sfc_name2(nvar2o+20) = 'ffmm'
      sfc_name2(nvar2o+21) = 'ffhh'
      sfc_name2(nvar2o+22) = 'hice'
      sfc_name2(nvar2o+23) = 'fice'
      sfc_name2(nvar2o+24) = 'tisfc'
      sfc_name2(nvar2o+25) = 'tprcp'
      sfc_name2(nvar2o+26) = 'srflag'
      sfc_name2(nvar2o+27) = 'snwdph'  !snowd
      sfc_name2(nvar2o+28) = 'shdmin'
      sfc_name2(nvar2o+29) = 'shdmax'
      sfc_name2(nvar2o+30) = 'slope'
      sfc_name2(nvar2o+31) = 'snoalb'
      !--- variables below here are optional
      sfc_name2(nvar2o+32) = 'sncovr'


      !--- register the 2D fields
      do num = 1,(nvar2o+nvar2s)
         var2_p => sfc_var2(:,:,num)
         id_restart = register_restart_field(Sfc_restart, fn_srf_TEST, sfc_name2(num), var2_p, domain=land_domain)
      enddo
      nullify(var2_p)
      
      !--- names of the 3D variables to save
      sfc_name3(1) = 'stc'
      sfc_name3(2) = 'smc'
      sfc_name3(3) = 'slc'

      !--- register the 3D fields
      do num = 1,nvar3
         var3_p => sfc_var3(:,:,:,num)
         id_restart = register_restart_field(Sfc_restart, fn_srf_TEST, sfc_name3(num), var3_p, domain=land_domain)
      enddo
      nullify(var3_p)
        
        
   end if

   ix = 0
   do j = jsc, jec
      j1 = j - jsc + 1
      do i = isc, iec
         i1 = i - isc +1
         ix = ix + 1
         
         sfc_var2(i1,j1,1) = noah_model%sfcprop%landfrac(ix) !--- slmsk
         sfc_var2(i1,j1,1)         = noah_model%sfcprop%landfrac(ix)
         sfc_var2(i1,j1,nvar2o+1)  = noah_model%sfcprop%slmsk(ix)
         sfc_var2(i1,j1,nvar2o+2)  = noah_model%sfcprop%tsfcl(ix)   
         sfc_var2(i1,j1,nvar2o+3)  = noah_model%sfcprop%weasd(ix)  !weasd, sheleg
         sfc_var2(i1,j1,nvar2o+4)  = noah_model%sfcprop%tg3(ix)
         sfc_var2(i1,j1,nvar2o+5)  = noah_model%sfcprop%zorll(ix)
         sfc_var2(i1,j1,nvar2o+6)  = noah_model%sfcprop%alvsf(ix)
         sfc_var2(i1,j1,nvar2o+7)  = noah_model%sfcprop%alvwf(ix)
         sfc_var2(i1,j1,nvar2o+8)  = noah_model%sfcprop%alnsf(ix)
         sfc_var2(i1,j1,nvar2o+9)  = noah_model%sfcprop%alnwf(ix)
         sfc_var2(i1,j1,nvar2o+10) = noah_model%sfcprop%facsf(ix)
         sfc_var2(i1,j1,nvar2o+11) = noah_model%sfcprop%facwf(ix)
         sfc_var2(i1,j1,nvar2o+12) = noah_model%sfcprop%vfrac(ix)
         sfc_var2(i1,j1,nvar2o+13) = noah_model%sfcprop%canopy(ix)
         sfc_var2(i1,j1,nvar2o+14) = noah_model%sfcprop%f10m(ix)
         sfc_var2(i1,j1,nvar2o+15) = noah_model%sfcprop%t2m(ix)
         sfc_var2(i1,j1,nvar2o+16) = noah_model%sfcprop%q2m(ix)
         sfc_var2(i1,j1,nvar2o+17) = noah_model%sfcprop%vtype(ix)
         sfc_var2(i1,j1,nvar2o+18) = noah_model%sfcprop%stype(ix)
         sfc_var2(i1,j1,nvar2o+19) = noah_model%sfcprop%uustar(ix)
         sfc_var2(i1,j1,nvar2o+20) = noah_model%sfcprop%ffmm(ix)
         sfc_var2(i1,j1,nvar2o+21) = noah_model%sfcprop%ffhh(ix)
         sfc_var2(i1,j1,nvar2o+22) = noah_model%sfcprop%hice(ix)
         sfc_var2(i1,j1,nvar2o+23) = noah_model%sfcprop%fice(ix)
         sfc_var2(i1,j1,nvar2o+24) = noah_model%sfcprop%tisfc(ix)
         sfc_var2(i1,j1,nvar2o+25) = noah_model%sfcprop%tprcp(ix)
         sfc_var2(i1,j1,nvar2o+26) = noah_model%sfcprop%srflag(ix)
         sfc_var2(i1,j1,nvar2o+27) = noah_model%sfcprop%snowd(ix)  !snowd,snwdph
         sfc_var2(i1,j1,nvar2o+28) = noah_model%sfcprop%shdmin(ix)
         sfc_var2(i1,j1,nvar2o+29) = noah_model%sfcprop%shdmax(ix)
         sfc_var2(i1,j1,nvar2o+30) = noah_model%sfcprop%slope(ix)
         sfc_var2(i1,j1,nvar2o+31) = noah_model%sfcprop%snoalb(ix)
         !--- variables below here are optional
         sfc_var2(i1,j1,nvar2o+32) = noah_model%sfcprop%sncovr(ix)
         
         !--- 3D variables
         do lsoil = 1,noah_model%static%km
            sfc_var3(i1,j1,lsoil,1) = noah_model%sfcprop%stc(ix,lsoil) !--- stc
            sfc_var3(i1,j1,lsoil,2) = noah_model%sfcprop%smc(ix,lsoil) !--- smc
            sfc_var3(i1,j1,lsoil,3) = noah_model%sfcprop%slc(ix,lsoil) !--- slc
         enddo

      enddo
   enddo

   call save_restart(Sfc_restart, timestamp)

 end subroutine sfc_prop_restart_write

 !-----------------------------------------------------------------------

 subroutine sfc_prop_transfer(noah_model)

   type (noah_type),          intent(inout) :: noah_model
   
   ! local
   ! ----------------------------------------
   integer :: i
   
   integer         :: isot
   integer         :: soiltyp   (noah_model%static%im) ! soil type (integer index)      
   integer         :: vegtype   (noah_model%static%im) ! vegetation type (integer index)    
   integer         :: slopetyp  (noah_model%static%im) ! class of sfc slope (integer index)

   real(kind_phys) :: vfrac     (noah_model%static%im)
   real(kind_phys) :: stype     (noah_model%static%im)
   real(kind_phys) :: vtype     (noah_model%static%im)
   real(kind_phys) :: slope     (noah_model%static%im)


   ! Some vars need more work than just copying. Associate these
   associate(                                      &
        ! static
        im         => noah_model%static%im        ,&
        isot       => noah_model%static%isot      ,&
        ivegsrc    => noah_model%static%ivegsrc   ,&
        ! sfcprop
        islmsk     => noah_model%sfcprop%slmsk    ,&
        stype      => noah_model%sfcprop%stype    ,&
        vtype      => noah_model%sfcprop%vtype    ,&
        slope      => noah_model%sfcprop%slope    ,&
        ! model
        soiltyp    => noah_model%model%soiltyp    ,&
        vegtype    => noah_model%model%vegtype    ,&
        slopetyp   => noah_model%model%slopetyp   ,&
        sigmaf     => noah_model%model%sigmaf     ,&        
        land       => noah_model%model%land        &        
        )
     
   
    !noah_model%model%slmsk  = noah_model%sfcprop%slmsk  
    noah_model%model%tskin   = noah_model%sfcprop%tsfcl  ! surface_skin_temperature_over_land
    noah_model%model%weasd  = noah_model%sfcprop%weasd   
    noah_model%model%tg3    = noah_model%sfcprop%tg3    
    noah_model%model%z0rl   = noah_model%sfcprop%zorll  ! surface_roughness_length_over_land
    
    ! noah_model%model%alvsf  = noah_model%sfcprop%alvsf  
    ! noah_model%model%alvwf  = noah_model%sfcprop%alvwf  
    ! noah_model%model%alnsf  = noah_model%sfcprop%alnsf  
    ! noah_model%model%alnwf  = noah_model%sfcprop%alnwf  
    ! noah_model%model%facsf  = noah_model%sfcprop%facsf  
    ! noah_model%model%facwf  = noah_model%sfcprop%facwf
    
    !noah_model%model%vfrac  = noah_model%sfcprop%vfrac

    !! This is copied from ccpp's GFS_surface_generic_pre. Be cafeful of code drift.
    !! land from islmsk is from GFS_surface_composites_pre
    !! TODO: Common code should be a shared module with CCPP's GFS_surface_generic
    
    do i=1,im
       sigmaf(i) = max(vfrac(i), 0.01_kind_phys)
       if (islmsk(i) == 2) then
          if (isot == 1) then
             soiltyp(i)  = 16
          else
             soiltyp(i)  = 9
          endif
          if (ivegsrc == 0 .or. ivegsrc == 4) then
             vegtype(i)  = 24
          elseif (ivegsrc == 1) then
             vegtype(i)  = 15
          elseif (ivegsrc == 2) then
             vegtype(i)  = 13
          elseif (ivegsrc == 3 .or. ivegsrc == 5) then
             vegtype(i)  = 15
          endif
          slopetyp(i) = 9
       else
          soiltyp(i)  = int( stype(i)+0.5_kind_phys )
          vegtype(i)  = int( vtype(i)+0.5_kind_phys )
          slopetyp(i) = int( slope(i)+0.5_kind_phys )    !! clu: slope -> slopetyp
          if (soiltyp(i)  < 1) soiltyp(i)  = 14
          if (vegtype(i)  < 1) vegtype(i)  = 17
          if (slopetyp(i) < 1) slopetyp(i) = 1
       endif

       ! set land mask
       if (islmsk(i) == 1) then
          land(i)    = .true.
       end if
    end do
    
    noah_model%model%canopy = noah_model%sfcprop%canopy 
    !noah_model%model%f10m   = noah_model%sfcprop%f10m   ! ratio_of_wind_at_lowest_model_layer_and_wind_at_10m
    !noah_model%model%t2m    = noah_model%sfcprop%t2m    
    !noah_model%model%q2m    = noah_model%sfcprop%q2m    

    !noah_model%model%vtype  = noah_model%sfcprop%vtype  
    !noah_model%model%stype  = noah_model%sfcprop%stype
    !noah_model%model%slope  = noah_model%sfcprop%slope  

    ! TODO: These vars from restarts are not pure land. Ex, in GFS_surface_composites_post_run:
    ! ffmm(i)   = txl*ffmm_lnd(i)   + txi*ffmm_ice(i)   + txo*ffmm_wat(i)
    ! ffhh(i)   = txl*ffhh_lnd(i)   + txi*ffhh_ice(i)   + txo*ffhh_wat(i)
    ! uustar(i) = txl*uustar_lnd(i) + txi*uustar_ice(i) + txo*uustar_wat(i)
    ! But are they needed anyways? Not if only ouptuts of stability
    
    noah_model%model%ustar = noah_model%sfcprop%uustar  ! note GLOBAL surface friction velocity
    !noah_model%model%ffmm   = noah_model%sfcprop%ffmm   ! note GLOBAL Monin-Obukhov similarity function for momentum
    !noah_model%model%ffhh   = noah_model%sfcprop%ffhh   ! note GLOBAL

    ! These are not needed
    ! TODO: delete from restart reading
    !noah_model%model%hice   = noah_model%sfcprop%hice   
    !noah_model%model%fice   = noah_model%sfcprop%fice   
    !noah_model%model%tisfc  = noah_model%sfcprop%tisfc
    
    noah_model%model%tprcp  = noah_model%sfcprop%tprcp  
    noah_model%model%srflag = noah_model%sfcprop%srflag 
    noah_model%model%snwdph = noah_model%sfcprop%snowd  ! note GLOBAL surface_snow_thickness_water_equivalent
    noah_model%model%shdmin = noah_model%sfcprop%shdmin 
    noah_model%model%shdmax = noah_model%sfcprop%shdmax 
    noah_model%model%snoalb = noah_model%sfcprop%snoalb 
    noah_model%model%sncovr1= noah_model%sfcprop%sncovr 

    noah_model%model%stc    = noah_model%sfcprop%stc
    noah_model%model%smc    = noah_model%sfcprop%smc
    noah_model%model%slc    = noah_model%sfcprop%slc
    
  end associate
  
 end subroutine sfc_prop_transfer

end module land_restart_mod
