!>@brief The module 'get_stochy_pattern_mod' contains the subroutines to retrieve the random pattern in the cubed-sphere grid
module get_stochy_pattern_mod
 use kinddef, only : kind_dbl_prec, kind_evod
 use spectral_transforms, only : len_trie_ls,                       &
                                 len_trio_ls, ls_dim, stochy_la2ga,          &
                                 coslat_a, latg, levs, lonf, skeblevs,&
                                 four_to_grid, spec_to_four, dezouv_stochy,dozeuv_stochy
 use stochy_namelist_def, only : n_var_lndp, ntrunc, stochini,n_var_spp
 use stochy_data_mod, only : gg_lats, gg_lons, inttyp, nskeb, nshum, nsppt, &
                             nocnsppt,nepbl,nlndp,                          &
                             rnlat, rpattern_sfc, rpattern_skeb,            &
                             rpattern_shum, rpattern_sppt, rpattern_ocnsppt,&
                             rpattern_epbl1, rpattern_epbl2, skebu_save,    &
                             nspp,rpattern_spp,                             &
                             skebv_save, skeb_vwts, skeb_vpts, wlon
 use stochy_patterngenerator_mod, only: random_pattern, ndimspec,           &
                                        patterngenerator_advance
 use stochy_internal_state_mod, only: stochy_internal_state
 use mpi_wrapper, only : mp_reduce_sum,is_rootpe
 use mersenne_twister, only: random_seed
 implicit none
 private

 public  get_random_pattern_vector,get_random_pattern_spp 
 public  get_random_pattern_sfc,get_random_pattern_scalar
 public  write_stoch_restart_atm,write_stoch_restart_ocn
 logical :: first_call=.true.
 contains

!>@brief The subroutine 'get_random_pattern_sfc' converts spherical harmonics to the gaussian grid then interpolates to the target grid
!>@details This subroutine is for a 2-D (lat-lon) scalar field
subroutine get_random_pattern_sfc(rpattern,npatterns,&
           gis_stochy,pattern_3d)
!\callgraph

! generate a random pattern for stochastic physics
 implicit none
 type(random_pattern), intent(inout) :: rpattern(npatterns)
 type(stochy_internal_state), intent(in) :: gis_stochy
 integer,intent(in)::   npatterns

 integer i,j,lat,n,k
 real(kind=kind_dbl_prec), dimension(lonf,gis_stochy%lats_node_a,1):: wrk2d

! logical lprint

 real(kind=kind_dbl_prec), allocatable, dimension(:,:) :: workg
 real (kind=kind_dbl_prec)   glolal(lonf,gis_stochy%lats_node_a)
 integer kmsk0(lonf,gis_stochy%lats_node_a)
 real(kind=kind_dbl_prec),intent(out) :: pattern_3d(gis_stochy%nx,gis_stochy%ny,n_var_lndp)
 real(kind=kind_dbl_prec) :: pattern_1d(gis_stochy%nx)

 do k=1,n_var_lndp
   kmsk0 = 0
   glolal = 0.
   do n=1,npatterns
     call patterngenerator_advance(rpattern(n),k,.false.)
!     if (is_rootpe()) print *, 'Random pattern for LNDP PERTS in get_random_pattern_fv3_sfc: k, min, max ',k,minval(rpattern_sfc(n)%spec_o(:,:,k)), maxval(rpattern_sfc(n)%spec_o(:,:,k))
     call scalarspect_to_gaugrid(rpattern(n),gis_stochy,wrk2d,k)
     glolal = glolal + wrk2d(:,:,1)
   enddo

   allocate(workg(lonf,latg))
   workg = 0.
   do j=1,gis_stochy%lats_node_a
     lat=gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+j)
     do i=1,lonf
        workg(i,lat) = glolal(i,j)
     enddo
   enddo

   call mp_reduce_sum(workg,lonf,latg)
!   if (is_rootpe()) print *, 'workg after mp_reduce_sum for LNDP PERTS in get_random_pattern_fv3_sfc: k, min, max ',k,minval(workg), maxval(workg)

! interpolate to cube grid

   do j=1,gis_stochy%ny
      pattern_1d = 0
      associate( tlats=>gis_stochy%parent_lats(1:gis_stochy%len(j),j),&
                 tlons=>gis_stochy%parent_lons(1:gis_stochy%len(j),j))
      call stochy_la2ga(workg,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                        pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
      pattern_3d(:,j,k)=pattern_1d(:)
      end associate
   enddo
!   if (is_rootpe()) print *, '3D pattern for LNDP PERTS in get_random_pattern_fv3_sfc: k, min, max ',k,minval(pattern_3d(:,:,k)), maxval(pattern_3d(:,:,k))
   deallocate(workg)

 enddo  ! loop over k, n_var_lndp

end subroutine get_random_pattern_sfc


!>@brief The subroutine 'get_random_pattern_fv3_vect' converts spherical harmonics to a vector on gaussian grid then interpolates to the target grid 
!>@details This subroutine is for a 2-D (lat-lon) vector field
subroutine get_random_pattern_vector(rpattern,npatterns,&
           gis_stochy,upattern_3d,vpattern_3d)
!\callgraph

! generate a random pattern for stochastic physics
 implicit none
 type(stochy_internal_state), intent(in) :: gis_stochy
 type(random_pattern), intent(inout) :: rpattern(npatterns)

 real(kind=kind_evod), dimension(len_trie_ls,2) ::  vrtspec_e,divspec_e
 real(kind=kind_evod), dimension(len_trio_ls,2) ::  vrtspec_o,divspec_o
 integer::   npatterns

 real(kind=kind_dbl_prec) :: upattern_3d(gis_stochy%nx,gis_stochy%ny,levs)
 real(kind=kind_dbl_prec) :: vpattern_3d(gis_stochy%nx,gis_stochy%ny,levs)
 real(kind=kind_dbl_prec) :: pattern_1d(gis_stochy%nx)
 integer i,j,lat,n,nn,k
 real(kind_dbl_prec), dimension(lonf,gis_stochy%lats_node_a,1):: wrk2du,wrk2dv

! logical lprint

 real, allocatable, dimension(:,:) :: workgu,workgv
 integer kmsk0(lonf,gis_stochy%lats_node_a)
 kmsk0 = 0
 allocate(workgu(lonf,latg))
 allocate(workgv(lonf,latg))
 divspec_e = 0; divspec_o = 0.
 if (first_call) then
    allocate(skebu_save(gis_stochy%nx,gis_stochy%ny,skeblevs))
    allocate(skebv_save(gis_stochy%nx,gis_stochy%ny,skeblevs))
    do k=2,skeblevs
       workgu = 0.
       workgv = 0.
       do n=1,npatterns
          if (.not. stochini) call patterngenerator_advance(rpattern(n),k,first_call)
      !   ke norm (convert streamfunction forcing to vorticity forcing)
          do nn=1,2
             vrtspec_e(:,nn) = gis_stochy%kenorm_e*rpattern(n)%spec_e(:,nn,k)
             vrtspec_o(:,nn) = gis_stochy%kenorm_o*rpattern(n)%spec_o(:,nn,k)
          enddo
        ! convert to winds
          call vrtdivspect_to_uvgrid( divspec_e,divspec_o,vrtspec_e,vrtspec_o,&
                 wrk2du,wrk2dv, gis_stochy)
          do i=1,lonf
             do j=1,gis_stochy%lats_node_a
                lat=gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+j)
                workgu(i,lat) = workgu(i,lat) + wrk2du(i,j,1)
                workgv(i,lat) = workgv(i,lat) + wrk2dv(i,j,1)
             enddo
          enddo
       enddo
       call mp_reduce_sum(workgu,lonf,latg)
       call mp_reduce_sum(workgv,lonf,latg)
   ! interpolate to cube grid
       do j=1,gis_stochy%ny
          pattern_1d = 0
          associate( tlats=>gis_stochy%parent_lats(1:gis_stochy%len(j),j),&
                     tlons=>gis_stochy%parent_lons(1:gis_stochy%len(j),j))
          call stochy_la2ga(workgu,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                           pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
          skebu_save(:,j,k)=pattern_1d(:)
          call stochy_la2ga(workgv,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                           pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
          skebv_save(:,j,k)=-1*pattern_1d(:)
          end associate
       enddo
    enddo
 endif
 do k=1,skeblevs-1
    skebu_save(:,:,k)=skebu_save(:,:,k+1)
    skebv_save(:,:,k)=skebv_save(:,:,k+1)
    do n=1,npatterns
       rpattern(n)%spec_e(:,:,k)=rpattern(n)%spec_e(:,:,k+1)
       rpattern(n)%spec_o(:,:,k)=rpattern(n)%spec_o(:,:,k+1)
    enddo
 enddo
 ! get pattern for last level
 workgu = 0.
 workgv = 0.
 do n=1,npatterns
       call patterngenerator_advance(rpattern(n),skeblevs,first_call)
  ! ke norm (convert streamfunction forcing to vorticity forcing)
    divspec_e = 0; divspec_o = 0.
    do nn=1,2
       vrtspec_e(:,nn) = gis_stochy%kenorm_e*rpattern(n)%spec_e(:,nn,skeblevs)
       vrtspec_o(:,nn) = gis_stochy%kenorm_o*rpattern(n)%spec_o(:,nn,skeblevs)
    enddo
  ! convert to winds
    call vrtdivspect_to_uvgrid(&
           divspec_e,divspec_o,vrtspec_e,vrtspec_o,&
           wrk2du,wrk2dv, gis_stochy)
    do i=1,lonf
       do j=1,gis_stochy%lats_node_a
          lat=gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+j)
          workgu(i,lat) = workgu(i,lat) + wrk2du(i,j,1)
          workgv(i,lat) = workgv(i,lat) + wrk2dv(i,j,1)
       enddo
    enddo
 enddo
 call mp_reduce_sum(workgu,lonf,latg)
 call mp_reduce_sum(workgv,lonf,latg)
 ! interpolate to cube grid
 do j=1,gis_stochy%ny
    pattern_1d = 0
    associate( tlats=>gis_stochy%parent_lats(:,j),&
               tlons=>gis_stochy%parent_lons(:,j))
    call stochy_la2ga(workgu,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                      pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
    skebu_save(:,j,skeblevs)=pattern_1d(:)
    call stochy_la2ga(workgv,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                      pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
    skebv_save(:,j,skeblevs)=-1*pattern_1d(:)
    end associate
  enddo
  deallocate(workgu)
  deallocate(workgv)
   ! interpolate in the vertical  ! consider moving to cubed sphere side,  more memory, but less interpolations
  do k=1,levs
     do j=1,gis_stochy%ny
        upattern_3d(:,j,k) = skeb_vwts(k,1)*skebu_save(:,j,skeb_vpts(k,1))+skeb_vwts(k,2)*skebu_save(:,j,skeb_vpts(k,2))
        vpattern_3d(:,j,k) = skeb_vwts(k,1)*skebv_save(:,j,skeb_vpts(k,1))+skeb_vwts(k,2)*skebv_save(:,j,skeb_vpts(k,2))
     enddo
  enddo
  first_call=.false.

end subroutine get_random_pattern_vector   

!>@brief The subroutine 'get_random_pattern_scalar' converts spherical harmonics to the gaussian grid then interpolates to the target grid
!>@details This subroutine is for a 2-D (lat-lon) scalar field
subroutine get_random_pattern_scalar(rpattern,npatterns,&
           gis_stochy,pattern_2d)

! generate a random pattern for stochastic physics
 implicit none
 type(random_pattern), intent(inout)  :: rpattern(npatterns)
 type(stochy_internal_state)          :: gis_stochy
 integer,intent(in)::   npatterns

 integer i,j,lat,n
 real(kind=kind_dbl_prec), dimension(lonf,gis_stochy%lats_node_a,1):: wrk2d

! logical lprint

 real(kind=kind_dbl_prec), allocatable, dimension(:,:) :: workg
 real (kind=kind_dbl_prec)   glolal(lonf,gis_stochy%lats_node_a)
 integer kmsk0(lonf,gis_stochy%lats_node_a)
 real(kind=kind_dbl_prec) :: pattern_2d(gis_stochy%nx,gis_stochy%ny)
 real(kind=kind_dbl_prec) :: pattern_1d(gis_stochy%nx)

 kmsk0 = 0
 glolal = 0.
 do n=1,npatterns
    call patterngenerator_advance(rpattern(n),1,.false.)
    call scalarspect_to_gaugrid(rpattern(n),gis_stochy,   &
         wrk2d,1)
    glolal = glolal + wrk2d(:,:,1)
 enddo

 allocate(workg(lonf,latg))
 workg = 0.
  do j=1,gis_stochy%lats_node_a
     lat=gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+j)
     do i=1,lonf
        workg(i,lat) = glolal(i,j)
     enddo
  enddo

   call mp_reduce_sum(workg,lonf,latg)

! interpolate to cube grid
   do j=1,gis_stochy%ny
      pattern_1d = 0
      associate( tlats=>gis_stochy%parent_lats(1:gis_stochy%len(j),j),&
                 tlons=>gis_stochy%parent_lons(1:gis_stochy%len(j),j))
      call stochy_la2ga(workg,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                        pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
      pattern_2d(:,j)=pattern_1d(:)
      end associate
   enddo
   deallocate(workg)

end subroutine get_random_pattern_scalar

!>@brief The subroutine 'get_random_pattern_spp' converts spherical harmonics
!to the gaussian grid then interpolates to the target grid
!>@details This subroutine is for a 2-D (lat-lon) scalar field
subroutine get_random_pattern_spp(rpattern,npatterns,&
           gis_stochy,pattern_3d)

! generate a random pattern for stochastic physics
 implicit none
 type(random_pattern), intent(inout)  :: rpattern(npatterns)
 type(stochy_internal_state)          :: gis_stochy
 integer,intent(in)::   npatterns

 integer i,j,lat,n

! logical lprint

 real(kind=kind_dbl_prec), allocatable, dimension(:,:) :: workg
 real (kind=kind_dbl_prec)   glolal(lonf,gis_stochy%lats_node_a)
 integer kmsk0(lonf,gis_stochy%lats_node_a)
 real(kind=kind_dbl_prec) :: pattern_3d(gis_stochy%nx,gis_stochy%ny,npatterns)
 real(kind=kind_dbl_prec) :: pattern_1d(gis_stochy%nx)

 allocate(workg(lonf,latg))
 do n=1,npatterns
    kmsk0 = 0
    glolal = 0.
    call patterngenerator_advance(rpattern(n),1,.false.)
    call scalarspect_to_gaugrid(rpattern(n),gis_stochy,   &
         glolal,1)

 workg = 0.
  do j=1,gis_stochy%lats_node_a
     lat=gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+j)
     do i=1,lonf
        workg(i,lat) = glolal(i,j)
     enddo
  enddo

   call mp_reduce_sum(workg,lonf,latg)

! interpolate to cube grid
   do j=1,gis_stochy%ny
      pattern_1d = 0
      associate( tlats=>gis_stochy%parent_lats(1:gis_stochy%len(j),j),&
                 tlons=>gis_stochy%parent_lons(1:gis_stochy%len(j),j))
      call stochy_la2ga(workg,lonf,latg,gg_lons,gg_lats,wlon,rnlat,&
                        pattern_1d(1:gis_stochy%len(j)),gis_stochy%len(j),tlats,tlons)
      pattern_3d(:,j,n)=pattern_1d(:)
      end associate
   enddo
 enddo
 deallocate(workg)

end subroutine get_random_pattern_spp

!>@brief The subroutine 'scalarspect_to_gaugrid' converts scalar spherical harmonics to a scalar on a gaussian grid
!>@details This subroutine is for a 2-D (lat-lon) scalar field
subroutine scalarspect_to_gaugrid(rpattern,gis_stochy,datag,n)
!\callgraph

      implicit none
      type(random_pattern),        intent(in)  :: rpattern
      type(stochy_internal_state), intent(in)  :: gis_stochy
      integer                 ,    intent(in)  :: n
      real(kind=kind_dbl_prec),    intent(out) :: datag(lonf,gis_stochy%lats_node_a)
! local vars
      real(kind=kind_dbl_prec) for_gr_a_1(gis_stochy%lon_dim_a,1,gis_stochy%lats_dim_a)
      real(kind=kind_dbl_prec) for_gr_a_2(lonf,1,gis_stochy%lats_dim_a)
      integer              i,k
      integer              lan,lat
      call spec_to_four(rpattern%spec_e(:,:,n), rpattern%spec_o(:,:,n), &
                  gis_stochy%plnev_a,gis_stochy%plnod_a,&
                  gis_stochy%ls_node, &
                  gis_stochy%lats_dim_a,for_gr_a_1,&
                  gis_stochy%ls_nodes,gis_stochy%max_ls_nodes,&
                  gis_stochy%lats_nodes_a,gis_stochy%global_lats_a,&
                  gis_stochy%lats_node_a,gis_stochy%ipt_lats_node_a,1)
      do lan=1,gis_stochy%lats_node_a
         lat = gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+lan)
         call four_to_grid(for_gr_a_1(:,:,lan),for_gr_a_2(:,:,lan),&
                           gis_stochy%lon_dim_a,1)
      enddo

      datag = 0.
      do lan=1,gis_stochy%lats_node_a
        lat      = gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+lan)
          do i=1,lonf
            datag(i,lan) = for_gr_a_2(i,1,lan)
        enddo
      enddo

      return
      end subroutine scalarspect_to_gaugrid


!>@brief The subroutine 'write_patterns' writes out a single pattern and the seed associated with the random number sequence in netcdf
!>@brief The subroutine 'write_stoch_restart_atm' writes out the speherical harmonics to a file, controlled by restart_interval
!>@details Only the active patterns are written out
subroutine write_stoch_restart_atm(sfile)
!\callgraph
    use netcdf
    use stochy_namelist_def, only : do_sppt,do_shum,do_skeb,lndp_type,do_spp
    implicit none
    character(len=*) :: sfile
    integer :: stochlun,k,n,isize,ierr
    integer :: ncid,varid1a,varid1b,varid2a,varid2b,varid3a,varid3b,varid4a,varid4b,varid5a,varid5b
    integer :: seed_dim_id,spec_dim_id,zt_dim_id,ztsfc_dim_id,np_dim_id,npsfc_dim_id
    integer :: ztspp_dim_id,npspp_dim_id

    include 'netcdf.inc'

    if ( ( .NOT. do_sppt) .AND. (.NOT. do_shum) .AND. (.NOT. do_skeb) .AND. (lndp_type==0 ) .AND. (.NOT. do_spp)) return
    stochlun=99
    if (is_rootpe()) then
       if (nsppt > 0 .OR. nshum > 0 .OR. nskeb > 0 .OR. nlndp>0 .OR. nspp>0 ) then
          ierr=nf90_create(trim(sfile),cmode=NF90_CLOBBER,ncid=ncid)
          ierr=NF90_PUT_ATT(ncid,NF_GLOBAL,"ntrunc",ntrunc)
          call random_seed(size=isize) ! get seed size
          ierr=NF90_DEF_DIM(ncid,"len_seed",isize,seed_dim_id)
          ierr=NF90_PUT_ATT(ncid,seed_dim_id,"long_name","length of random seed")
          ierr=NF90_DEF_DIM(ncid,"num_patterns",NF_UNLIMITED,np_dim_id) !  should be 5
          ierr=NF90_PUT_ATT(ncid,np_dim_id,"long_name","number of random patterns (max of 5)")
          if (lndp_type .NE. 0) then
             ierr=NF90_DEF_DIM(ncid,"num_patterns_sfc",nlndp,npsfc_dim_id) !  should be 5
             ierr=NF90_PUT_ATT(ncid,npsfc_dim_id,"long_name","number of random patterns for surface)")
             ierr=NF90_DEF_DIM(ncid,"n_var_lndp",n_var_lndp,ztsfc_dim_id)
             ierr=NF90_PUT_ATT(ncid,ztsfc_dim_id,"long_name","number of sfc perturbation types")
          endif
          if (nspp .GT. 0) then
             ierr=NF90_DEF_DIM(ncid,"num_patterns_spp",nspp,npspp_dim_id) !  should be 5
             ierr=NF90_PUT_ATT(ncid,npspp_dim_id,"long_name","number of random patterns for spp)")
             ierr=NF90_DEF_DIM(ncid,"n_var_spp",n_var_spp,ztspp_dim_id)
             ierr=NF90_PUT_ATT(ncid,ztspp_dim_id,"long_name","number of spp perturbation types")
          endif
          ierr=NF90_DEF_DIM(ncid,"ndimspecx2",2*ndimspec,spec_dim_id)
          ierr=NF90_PUT_ATT(ncid,spec_dim_id,"long_name","number of spectral cofficients")
          if (do_sppt) then
             ierr=NF90_DEF_VAR(ncid,"sppt_seed",NF90_DOUBLE,(/seed_dim_id, np_dim_id/), varid1a)
             ierr=NF90_PUT_ATT(ncid,varid1a,"long_name","random number seed for SPPT")
             ierr=NF90_DEF_VAR(ncid,"sppt_spec",NF90_DOUBLE,(/spec_dim_id, np_dim_id/), varid1b)
             ierr=NF90_PUT_ATT(ncid,varid1b,"long_name","spectral cofficients SPPT")
          endif
          if (do_shum) then
             ierr=NF90_DEF_VAR(ncid,"shum_seed",NF90_DOUBLE,(/seed_dim_id, np_dim_id/), varid2a)
             ierr=NF90_PUT_ATT(ncid,varid2a,"long_name","random number seed for SHUM")
             ierr=NF90_DEF_VAR(ncid,"shum_spec",NF90_DOUBLE,(/spec_dim_id, np_dim_id/), varid2b)
             ierr=NF90_PUT_ATT(ncid,varid2b,"long_name","spectral cofficients SHUM")
          endif
          if (do_skeb) then
             ierr=NF90_DEF_DIM(ncid,"skeblevs",skeblevs,zt_dim_id)
             ierr=NF90_PUT_ATT(ncid,zt_dim_id,"long_name","number of vertical levels for SKEB")
             ierr=NF90_DEF_VAR(ncid,"skeb_seed",NF90_DOUBLE,(/seed_dim_id, zt_dim_id,np_dim_id/), varid3a)
             ierr=NF90_PUT_ATT(ncid,varid3a,"long_name","random number seed for SKEB")
             ierr=NF90_DEF_VAR(ncid,"skeb_spec",NF90_DOUBLE,(/spec_dim_id, zt_dim_id,np_dim_id/), varid3b)
             ierr=NF90_PUT_ATT(ncid,varid3b,"long_name","spectral cofficients SKEB")
          endif
          if (nlndp>0) then
             ierr=NF90_DEF_VAR(ncid,"sfcpert_seed",NF90_DOUBLE,(/seed_dim_id, ztsfc_dim_id, npsfc_dim_id/), varid4a)
             ierr=NF90_PUT_ATT(ncid,varid4a,"long_name","random number seed for SHUM")
             ierr=NF90_DEF_VAR(ncid,"sfcpert_spec",NF90_DOUBLE,(/spec_dim_id, ztsfc_dim_id, npsfc_dim_id/), varid4b)
             ierr=NF90_PUT_ATT(ncid,varid4b,"long_name","spectral cofficients SHUM")
          endif
          if (nspp>0) then
             ierr=NF90_DEF_VAR(ncid,"spp_seed",NF90_DOUBLE,(/seed_dim_id, ztspp_dim_id, npspp_dim_id/), varid5a)
             ierr=NF90_PUT_ATT(ncid,varid5a,"long_name","random number seed for SPP")
             ierr=NF90_DEF_VAR(ncid,"spp_spec",NF90_DOUBLE,(/spec_dim_id, ztspp_dim_id, npspp_dim_id/), varid5b)
             ierr=NF90_PUT_ATT(ncid,varid5b,"long_name","spectral cofficients SPP")
          endif
          ierr=NF90_ENDDEF(ncid)
          if (ierr .NE. 0) then
             write(0,*) 'error creating stochastic restart file'
             return
         end if
       endif
    endif
    if (nsppt > 0) then
       do n=1,nsppt
          call write_pattern(rpattern_sppt(n),ncid,1,n,varid1a,varid1b,.false.,ierr)
       enddo
    endif
    if (nshum > 0) then
       do n=1,nshum
          call write_pattern(rpattern_shum(n),ncid,1,n,varid2a,varid2b,.false.,ierr)
       enddo
    endif
    if (nskeb > 0) then
       do n=1,nskeb
          do k=1,skeblevs
             call write_pattern(rpattern_skeb(n),ncid,k,n,varid3a,varid3b,.true.,ierr)
          enddo
       enddo
    endif
    if (lndp_type .NE. 0 .AND. nlndp>0) then
       do n=1,nlndp
          do k=1,n_var_lndp
             call write_pattern(rpattern_sfc(n),ncid,k,n,varid4a,varid4b,.true.,ierr)
          enddo
       enddo
    endif
    if (nspp > 0) then
       do n=1,nspp
          call write_pattern(rpattern_spp(n),ncid,1,n,varid5a,varid5b,.true.,ierr)
       enddo
    endif
    if (is_rootpe() ) then
       ierr=NF90_CLOSE(ncid)
       if (ierr .NE. 0) then
           write(0,*) 'error writing patterns and closing file'
           return
       endif
    endif
 end subroutine write_stoch_restart_atm


!>@brief The subroutine 'write_stoch_restart_ocn' writes out the speherical harmonics to a file, controlled by restart_interval
!>@details Only the active patterns are written out
subroutine write_stoch_restart_ocn(sfile)
!\callgraph
    use netcdf
    use stochy_namelist_def, only : do_ocnsppt,pert_epbl
    implicit none
    character(len=*) :: sfile
    integer :: stochlun,k,n,isize,ierr
    integer :: ncid,varid1a,varid1b,varid2a,varid2b,varid3a,varid3b
    integer :: seed_dim_id,spec_dim_id,np_dim_id
    include 'netcdf.inc'
    if ( ( .NOT. do_ocnsppt) .AND. (.NOT. pert_epbl) ) return
    stochlun=99
    if (is_rootpe()) then
       ierr=nf90_create(trim(sfile),cmode=NF90_CLOBBER,ncid=ncid)
       ierr=NF90_PUT_ATT(ncid,NF_GLOBAL,"ntrunc",ntrunc)
       call random_seed(size=isize) ! get seed size
       ierr=NF90_DEF_DIM(ncid,"len_seed",isize,seed_dim_id)
       ierr=NF90_PUT_ATT(ncid,seed_dim_id,"long_name","length of random seed")
       ierr=NF90_DEF_DIM(ncid,"num_patterns",NF_UNLIMITED,np_dim_id) !  should be 5
       ierr=NF90_PUT_ATT(ncid,np_dim_id,"long_name","number of random patterns (max of 5)")
       ierr=NF90_DEF_DIM(ncid,"ndimspecx2",2*ndimspec,spec_dim_id)
       ierr=NF90_PUT_ATT(ncid,spec_dim_id,"long_name","number of spectral cofficients")
       if (do_ocnsppt) then
          ierr=NF90_DEF_VAR(ncid,"ocnsppt_seed",NF90_DOUBLE,(/seed_dim_id, np_dim_id/), varid1a)
          ierr=NF90_PUT_ATT(ncid,varid1a,"long_name","random number seed for SPPT")
          ierr=NF90_DEF_VAR(ncid,"ocnsppt_spec",NF90_DOUBLE,(/spec_dim_id, np_dim_id/), varid1b)
          ierr=NF90_PUT_ATT(ncid,varid1b,"long_name","spectral cofficients SPPT")
       endif
       if (pert_epbl) then
          ierr=NF90_DEF_VAR(ncid,"epbl1_seed",NF90_DOUBLE,(/seed_dim_id, np_dim_id/), varid2a)
          ierr=NF90_PUT_ATT(ncid,varid2a,"long_name","random number seed for EPBL1")
          ierr=NF90_DEF_VAR(ncid,"epbl1_spec",NF90_DOUBLE,(/spec_dim_id, np_dim_id/), varid2b)
          ierr=NF90_PUT_ATT(ncid,varid2b,"long_name","spectral cofficients EPBL1")
          ierr=NF90_DEF_VAR(ncid,"epbl2_seed",NF90_DOUBLE,(/seed_dim_id, np_dim_id/), varid3a)
          ierr=NF90_PUT_ATT(ncid,varid3a,"long_name","random number seed for EPBL2")
          ierr=NF90_DEF_VAR(ncid,"epbl2_spec",NF90_DOUBLE,(/spec_dim_id, np_dim_id/), varid3b)
          ierr=NF90_PUT_ATT(ncid,varid3b,"long_name","spectral cofficients EPBL2")
       endif
       ierr=NF90_ENDDEF(ncid)
       if (ierr .NE. 0) then
          write(0,*) 'error creating stochastic restart file'
          return
       end if
    endif
    if (nocnsppt > 0) then
       do n=1,nocnsppt
          call write_pattern(rpattern_ocnsppt(n),ncid,1,n,varid1a,varid1b,.false.,ierr)
       enddo
    endif
    if (nepbl > 0) then
       do n=1,nepbl
          call write_pattern(rpattern_epbl1(n),ncid,1,n,varid2a,varid2b,.false.,ierr)
          call write_pattern(rpattern_epbl2(n),ncid,1,n,varid3a,varid3b,.false.,ierr)
       enddo
    endif
    if (is_rootpe() ) then
       ierr=NF90_CLOSE(ncid)
       if (ierr .NE. 0) then
           write(0,*) 'error writing patterns and closing file'
           return
       endif
    endif
 end subroutine write_stoch_restart_ocn
!>@brief The subroutine 'write_patterns' writes out a single pattern and the seed associated with the random number sequence
!>@details Spherical harminoncs are stored with trianglular truncation
 subroutine write_pattern(rpattern,outlun,lev,np,varid1,varid2,slice_of_3d,iret)
!\callgraph
   use netcdf
   implicit none
   type(random_pattern), intent(inout) :: rpattern
   integer, intent(in) :: outlun,lev
   integer, intent(in) :: np,varid1,varid2
   logical, intent(in) :: slice_of_3d
   integer, intent(out) :: iret
   real(kind_dbl_prec), allocatable  :: pattern2d(:)
   integer nm,nn,arrlen,isize,ierr
   integer,allocatable :: isave(:)
   include 'netcdf.inc'
   arrlen=2*ndimspec
   iret=0
   allocate(pattern2d(arrlen))
   pattern2d=0.0
   ! fill in apprpriate pieces of array
   do nn=1,len_trie_ls
      nm = rpattern%idx_e(nn)
      if (nm == 0) cycle
      pattern2d(nm)          = rpattern%spec_e(nn,1,lev)
      pattern2d(ndimspec+nm) = rpattern%spec_e(nn,2,lev)
   enddo
   do nn=1,len_trio_ls
      nm = rpattern%idx_o(nn)
      if (nm == 0) cycle
      pattern2d(nm)          = rpattern%spec_o(nn,1,lev)
      pattern2d(ndimspec+nm) = rpattern%spec_o(nn,2,lev)
   enddo
   call mp_reduce_sum(pattern2d,arrlen)
  !  write only on root process
   if (is_rootpe()) then
      print*,'writing out random pattern (min/max/size)',&
      minval(pattern2d),maxval(pattern2d),size(pattern2d)
      call random_seed(size=isize) ! get seed size
      allocate(isave(isize)) ! get seed
      call random_seed(get=isave,stat=rpattern%rstate) ! write seed
      ierr=NF90_PUT_VAR(outlun,varid1,isave,(/1,np/))
      if (slice_of_3d) then
         ierr=NF90_PUT_VAR(outlun,varid2,pattern2d,(/1,lev,np/))
      else
         ierr=NF90_PUT_VAR(outlun,varid2,pattern2d,(/1,np/))
      endif
      if (ierr .NE. 0) then
          write(0,*) 'error writing to stochastic restart file'
          iret = ierr
          return
      end if
   endif
   deallocate(pattern2d)
 end subroutine write_pattern
!>@brief The subroutine 'vrtdivspect_to_uvgrid' converts vorticty and divergence spherical harmonics to 
! zonal and meridional winds on the gaussian grid
!>@details This subroutine is for a 2-D (lat-lon) vector field
 subroutine vrtdivspect_to_uvgrid(&
           trie_di,trio_di,trie_ze,trio_ze,&
           uug,vvg, gis_stochy)
!\callgraph

      implicit none
      type(stochy_internal_state), intent(in) :: gis_stochy
      real(kind=kind_dbl_prec), intent(in) :: trie_di(len_trie_ls,2)
      real(kind=kind_dbl_prec), intent(in) :: trio_di(len_trio_ls,2)
      real(kind=kind_dbl_prec), intent(in) :: trie_ze(len_trie_ls,2)
      real(kind=kind_dbl_prec), intent(in) :: trio_ze(len_trio_ls,2)
      real(kind=kind_dbl_prec), intent(out) :: uug(lonf,gis_stochy%lats_node_a)
      real(kind=kind_dbl_prec), intent(out) :: vvg(lonf,gis_stochy%lats_node_a)
! local vars
      real(kind=kind_dbl_prec) trie_ls(len_trie_ls,2,2)
      real(kind=kind_dbl_prec) trio_ls(len_trio_ls,2,2)
      real(kind=kind_dbl_prec) for_gr_a_1(gis_stochy%lon_dim_a,2,gis_stochy%lats_dim_a)
      real(kind=kind_dbl_prec) for_gr_a_2(lonf,2,gis_stochy%lats_dim_a)
      integer              i,k
      integer              lan,lat
      real (kind=kind_dbl_prec) tx1

      call dezouv_stochy(trie_di(:,:),       trio_ze(:,:), &
                  trie_ls(:,:,1), trio_ls(:,:,2), gis_stochy%epsedn,gis_stochy%epsodn, &
                  gis_stochy%snnp1ev,gis_stochy%snnp1od,gis_stochy%ls_node)
      call dozeuv_stochy(trio_di(:,:),       trie_ze(:,:), &
                  trio_ls(:,:,1), trie_ls(:,:,2), gis_stochy%epsedn,gis_stochy%epsodn, &
                  gis_stochy%snnp1ev,gis_stochy%snnp1od,gis_stochy%ls_node)

      call spec_to_four(trie_ls, trio_ls, &
                  gis_stochy%plnev_a,gis_stochy%plnod_a,&
                  gis_stochy%ls_node,&
                  gis_stochy%lats_dim_a,for_gr_a_1,&
                  gis_stochy%ls_nodes,gis_stochy%max_ls_nodes,&
                  gis_stochy%lats_nodes_a,gis_stochy%global_lats_a,&
                  gis_stochy%lats_node_a,gis_stochy%ipt_lats_node_a,2)

      do lan=1,gis_stochy%lats_node_a
         lat = gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+lan)
         call four_to_grid(for_gr_a_1(:,:,lan),for_gr_a_2(:,:,lan),&
                           gis_stochy%lon_dim_a,2)
      enddo

      uug = 0.; vvg = 0.
      do lan=1,gis_stochy%lats_node_a
        lat      = gis_stochy%global_lats_a(gis_stochy%ipt_lats_node_a-1+lan)
        tx1      = 1. / coslat_a(lat)
        do i=1,lonf
          uug(i,lan) = for_gr_a_2(i,1,lan) * tx1
          vvg(i,lan) = for_gr_a_2(i,2,lan) * tx1
        enddo
      enddo

      return
 end subroutine vrtdivspect_to_uvgrid
end module get_stochy_pattern_mod
