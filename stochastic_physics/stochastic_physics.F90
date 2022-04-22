!>@brief The module 'stochastic_physics' is for initialization and running of
!! the stochastic physics random pattern generators
module stochastic_physics

use kinddef, only : kind_dbl_prec

implicit none

private

public :: init_stochastic_physics,init_stochastic_physics_ocn
public :: run_stochastic_physics,run_stochastic_physics_ocn
public :: finalize_stochastic_physics

contains

!>@brief The subroutine 'init_stochastic_physics' initializes the stochastic
!!pattern genertors
!>@details It reads the stochastic physics namelist (nam_stoch and nam_sfcperts)
!allocates and polulates the necessary arrays

subroutine init_stochastic_physics(levs, blksz, dtp, sppt_amp, input_nml_file_in, fn_nml, nlunit, &
    xlon,xlat, &
    do_sppt_in, do_shum_in, do_skeb_in, lndp_type_in, n_var_lndp_in, use_zmtnblck_out, skeb_npass_out,    &
    lndp_var_list_out, lndp_prt_list_out,                         &
    n_var_spp_in, spp_var_list_out, spp_prt_list_out, spp_stddev_cutoff_out, do_spp_in,     &
    ak, bk, nthreads, mpiroot, mpicomm, iret) 
!\callgraph
!use stochy_internal_state_moa
use stochy_data_mod, only : init_stochdata,gg_lats,gg_lons,nsppt, &
                            rad2deg,INTTYP,wlon,rnlat,gis_stochy, &
                            vfact_skeb,vfact_sppt,vfact_shum,skeb_vpts,skeb_vwts,sl, &
                            nspp, vfact_spp
use stochy_namelist_def
use spectral_transforms,only:colrad_a,latg,lonf,skeblevs
use mpi_wrapper, only : mpi_wrapper_initialize,mype,npes,is_rootpe

implicit none
integer, intent(out)                    :: iret

! Interface variables

integer,                  intent(in)    :: levs, nlunit, nthreads, mpiroot, mpicomm
integer,                  intent(in)    :: blksz(:)
real(kind=kind_dbl_prec), intent(in)    :: dtp
real(kind=kind_dbl_prec), intent(out)   :: sppt_amp
character(len=*),         intent(in)    :: input_nml_file_in(:)
character(len=*),         intent(in)    :: fn_nml
real(kind=kind_dbl_prec), intent(in)    :: xlon(:,:)
real(kind=kind_dbl_prec), intent(in)    :: xlat(:,:)
logical,                  intent(in)    :: do_sppt_in, do_shum_in, do_skeb_in ,do_spp_in
integer,                  intent(in)    :: lndp_type_in, n_var_lndp_in
integer,                  intent(in)    :: n_var_spp_in
real(kind=kind_dbl_prec), intent(in)    :: ak(:), bk(:) 
logical,                  intent(out)   :: use_zmtnblck_out
integer,                  intent(out)   :: skeb_npass_out
character(len=3),         dimension(:), intent(out) :: lndp_var_list_out
real(kind=kind_dbl_prec), dimension(:), intent(out) :: lndp_prt_list_out
character(len=3),         dimension(:), intent(out) :: spp_var_list_out
real(kind=kind_dbl_prec), dimension(:), intent(out) :: spp_prt_list_out
real(kind=kind_dbl_prec), dimension(:), intent(out) :: spp_stddev_cutoff_out


! Local variables
real(kind=kind_dbl_prec), parameter     :: con_pi =4.0d0*atan(1.0d0)
integer :: nblks,len
real*8 :: PRSI(levs),PRSL(levs),dx
real, allocatable :: skeb_vloc(:)
integer :: k,kflip,latghf,blk,k2,v,i
character*2::proc

! Initialize MPI and OpenMP
call mpi_wrapper_initialize(mpiroot,mpicomm)
gis_stochy%nodes = npes
gis_stochy%mype=mype
gis_stochy%nx=maxval(blksz)
nblks = size(blksz)
gis_stochy%ny=nblks
rad2deg=180.0/con_pi

! ------------------------------------------

nblks = size(blksz)
allocate(gis_stochy%len(nblks))
allocate(gis_stochy%parent_lons(gis_stochy%nx,gis_stochy%ny))
allocate(gis_stochy%parent_lats(gis_stochy%nx,gis_stochy%ny))
do blk=1,nblks
   len=blksz(blk)
   gis_stochy%parent_lons(1:len,blk)=xlon(blk,1:len)*rad2deg
   gis_stochy%parent_lats(1:len,blk)=xlat(blk,1:len)*rad2deg
   gis_stochy%len(blk)=len
enddo

! replace
INTTYP=0 ! bilinear interpolation
call init_stochdata(levs,dtp,input_nml_file_in,fn_nml,nlunit,iret)
if (iret .ne. 0) return
! check namelist entries for consistency
if (do_sppt_in.neqv.do_sppt) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_sppt and sppt'
   iret = 20 
   return
else if (do_shum_in.neqv.do_shum) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_shum and shum'
   iret = 20 
   return
else if (do_skeb_in.neqv.do_skeb) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_skeb and skeb'
   iret = 20 
   return
else if (lndp_type_in /= lndp_type) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings lndp_type in physics and nam_sfcperts'
   iret = 20 
   return
else if (n_var_lndp_in /=  n_var_lndp) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings n_var_lndp in physics nml, and lndp_* in nam_sfcperts'
   iret = 20 
   return
else if (n_var_spp_in .ne. n_var_spp) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings n_var_spp in physics nml, and spp_* in nam_sppperts'
   write(0,*) 'n_var_spp, n_var_spp_in', n_var_spp, n_var_spp_in
   iret = 20
   return
else if (do_spp_in.neqv.do_spp) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_init: incompatible', &
                   & ' namelist settings do_spp and spp'
   iret = 20
   return
end if
! update remaining model configuration parameters from namelist
use_zmtnblck_out=use_zmtnblck
skeb_npass_out=skeb_npass
if (n_var_lndp>0) then
   lndp_var_list_out=lndp_var_list(1:n_var_lndp)
   lndp_prt_list_out=lndp_prt_list(1:n_var_lndp)
endif
if (n_var_spp>0) then
   spp_var_list_out=spp_var_list(1:n_var_spp)
   spp_prt_list_out=spp_prt_list(1:n_var_spp)
   spp_stddev_cutoff_out=spp_stddev_cutoff(1:n_var_spp)
endif
if ( (.NOT. do_sppt) .AND. (.NOT. do_shum) .AND. (.NOT. do_skeb)  .AND. (lndp_type==0) .AND. (.NOT. do_spp)) return
allocate(sl(levs))
do k=1,levs
   sl(k)= 0.5*(ak(k)/101300.+bk(k)+ak(k+1)/101300.0+bk(k+1)) ! si are now sigmas
enddo
if (do_sppt) then
   allocate(vfact_sppt(levs))
   do k=1,levs
      if (sl(k) .lt. sppt_sigtop1 .and. sl(k) .gt. sppt_sigtop2) then
         vfact_sppt(k) = (sl(k)-sppt_sigtop2)/(sppt_sigtop1-sppt_sigtop2)
      else if (sl(k) .lt. sppt_sigtop2) then
          vfact_sppt(k) = 0.0
      else
          vfact_sppt(k) = 1.0
      endif
   enddo
   if (sppt_sfclimit) then
       vfact_sppt(2)=vfact_sppt(3)*0.5
       vfact_sppt(1)=0.0
   endif
   if (is_rootpe()) then
      do k=1,levs
         print *,'sppt vert profile',k,sl(k),vfact_sppt(k)
      enddo
   endif
   sppt_amp=sqrt(SUM(sppt(1:nsppt)**2))
endif
if (do_skeb) then
   allocate(vfact_skeb(levs))
   allocate(skeb_vloc(skeblevs)) ! local
   allocate(skeb_vwts(levs,2)) ! save for later
   allocate(skeb_vpts(levs,2)) ! save for later
   do k=1,levs
      if (sl(k) .lt. skeb_sigtop1 .and. sl(k) .gt. skeb_sigtop2) then
         vfact_skeb(k) = (sl(k)-skeb_sigtop2)/(skeb_sigtop1-skeb_sigtop2)
      else if (sl(k) .lt. skeb_sigtop2) then
          vfact_skeb(k) = 0.0
      else
          vfact_skeb(k) = 1.0
      endif
      if (is_rootpe())  print *,'skeb vert profile',k,sl(k),vfact_skeb(k)
   enddo
! calculate vertical interpolation weights
   do k=1,skeblevs
      skeb_vloc(k)=sl(1)-real(k-1)/real(skeblevs-1.0)*(sl(1)-sl(levs))
   enddo
! surface
skeb_vwts(1,2)=0
skeb_vpts(1,1)=1
! top
skeb_vwts(levs,2)=1
skeb_vpts(levs,1)=skeblevs-2
! internal
DO k=2,levs-1
   DO k2=1,skeblevs-1
      IF (sl(k) .LE. skeb_vloc(k2) .AND. sl(k) .GT. skeb_vloc(k2+1)) THEN
        skeb_vpts(k,1)=k2
        skeb_vwts(k,2)=(skeb_vloc(k2)-sl(k))/(skeb_vloc(k2)-skeb_vloc(k2+1))
      ENDIF
   ENDDO
ENDDO
deallocate(skeb_vloc)
if (is_rootpe()) then
DO k=1,levs
   print*,'skeb vpts ',skeb_vpts(k,1),skeb_vwts(k,2)
ENDDO
endif
skeb_vwts(:,1)=1.0-skeb_vwts(:,2)
skeb_vpts(:,2)=skeb_vpts(:,1)+1.0
endif

if (do_shum) then
   allocate(vfact_shum(levs))
   do k=1,levs
      vfact_shum(k) = exp((sl(k)-1.)/shum_sigefold)
      if (sl(k).LT. 2*shum_sigefold) then
         vfact_shum(k)=0.0
      endif
      if (is_rootpe())  print *,'shum vert profile',k,sl(k),vfact_shum(k)
   enddo
endif
if (do_spp) then
   allocate(vfact_spp(levs))
   do k=1,levs
      if (sl(k) .lt. spp_sigtop1(1) .and. sl(k) .gt. spp_sigtop2(1)) then
         vfact_spp(k) = (sl(k)-spp_sigtop2(1))/(spp_sigtop1(1)-spp_sigtop2(1))
      else if (sl(k) .lt. spp_sigtop2(1)) then
          vfact_spp(k) = 0.0
      else
          vfact_spp(k) = 1.0
      endif
      if (is_rootpe())  print *,'spp vert profile',k,sl(k),vfact_spp(k)
   enddo
endif
! get interpolation weights
! define gaussian grid lats and lons
latghf=latg/2
allocate(gg_lats(latg))
allocate(gg_lons(lonf))
do k=1,latghf
   gg_lats(k)=-1.0*colrad_a(latghf-k+1)*rad2deg
   gg_lats(latg-k+1)=-1*gg_lats(k)
enddo
dx=360.0/lonf
do k=1,lonf
  gg_lons(k)=dx*(k-1)
enddo
WLON=gg_lons(1)-(gg_lons(2)-gg_lons(1))
RNLAT=gg_lats(1)*2-gg_lats(2)

end subroutine init_stochastic_physics

!!!!!!!!!!!!!!!!!!!!
subroutine init_stochastic_physics_ocn(delt,geoLonT,geoLatT,nx,ny,nz,pert_epbl_in,do_sppt_in, &
                                       mpiroot, mpicomm, iret)
use stochy_data_mod, only : init_stochdata_ocn,gg_lats,gg_lons,&
                            rad2deg,INTTYP,wlon,rnlat,gis_stochy_ocn
use spectral_transforms , only : latg,lonf,colrad_a
!use MOM_grid, only : ocean_grid_type   
use stochy_namelist_def
use mersenne_twister, only: random_gauss
use mpi_wrapper, only : mpi_wrapper_initialize,mype,npes,is_rootpe

implicit none
real,intent(in)  :: delt
integer,intent(in) :: nx,ny,nz
real,intent(in) :: geoLonT(nx,ny),geoLatT(nx,ny)
logical,intent(in) :: pert_epbl_in,do_sppt_in
integer,intent(in)    :: mpiroot, mpicomm
integer, intent(out) :: iret
real(kind=kind_dbl_prec), parameter     :: con_pi =4.0d0*atan(1.0d0)

real :: dx
integer :: k,latghf,km
rad2deg=180.0/con_pi
call mpi_wrapper_initialize(mpiroot,mpicomm)
gis_stochy_ocn%nodes = npes
gis_stochy_ocn%mype = mype
gis_stochy_ocn%nx=nx  
gis_stochy_ocn%ny=ny
allocate(gis_stochy_ocn%len(ny))
allocate(gis_stochy_ocn%parent_lons(nx,ny))
allocate(gis_stochy_ocn%parent_lats(nx,ny))
gis_stochy_ocn%len(:)=nx
gis_stochy_ocn%parent_lons=geoLonT
gis_stochy_ocn%parent_lats=geoLatT

INTTYP=0 ! bilinear interpolation
km=nz
call init_stochdata_ocn(km,delt,iret)
if (do_sppt_in.neqv.do_ocnsppt) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_ocn_init: incompatible', &
                   & ' namelist settings do_sppt and sppt'
   iret = 20 
   return
else if (pert_epbl_in.neqv.pert_epbl) then
   write(0,'(*(a))') 'Logic error in stochastic_physics_ocn_init: incompatible', &
                   & ' namelist settings pert_epbl and epbl'
   iret = 20 
   return
end if

! get interpolation weights
! define gaussian grid lats and lons
latghf=latg/2
allocate(gg_lats(latg))
allocate(gg_lons(lonf))
do k=1,latghf
   gg_lats(k)=-1.0*colrad_a(latghf-k+1)*rad2deg
   gg_lats(latg-k+1)=-1*gg_lats(k)
enddo
dx=360.0/lonf
do k=1,lonf
  gg_lons(k)=dx*(k-1)
enddo
WLON=gg_lons(1)-(gg_lons(2)-gg_lons(1))
RNLAT=gg_lats(1)*2-gg_lats(2)
end subroutine init_stochastic_physics_ocn

!!!!!!!!!!!!!!!!!!!!


!>@brief The subroutine 'run_stochastic_physics' updates the random patterns if
!!necessary
!>@details It updates the AR(1) in spectral space
!allocates and polulates the necessary arrays

subroutine run_stochastic_physics(levs, kdt, fhour, blksz, sppt_wts, shum_wts, skebu_wts,  & 
                                  skebv_wts, sfc_wts, spp_wts, nthreads)

!\callgraph
!use stochy_internal_state_mod
use stochy_data_mod, only : nshum,rpattern_shum,rpattern_sppt,nsppt,rpattern_skeb,nskeb,&
                            gis_stochy,vfact_sppt,vfact_shum,vfact_skeb, rpattern_sfc, nlndp, &
                            rpattern_spp, nspp, vfact_spp
use get_stochy_pattern_mod,only : get_random_pattern_scalar,get_random_pattern_vector, & 
                                  get_random_pattern_sfc,get_random_pattern_spp
use stochy_namelist_def, only : do_shum,do_sppt,do_skeb,nssppt,nsshum,nsskeb,sppt_logit,    & 
                                lndp_type, n_var_lndp, n_var_spp, do_spp, spp_stddev_cutoff, spp_prt_list
use mpi_wrapper, only: is_rootpe
implicit none

! Interface variables
integer,                  intent(in) :: levs, kdt
real(kind=kind_dbl_prec), intent(in) :: fhour
integer,                  intent(in) :: blksz(:)
real(kind=kind_dbl_prec), intent(inout) :: sppt_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: shum_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: skebu_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: skebv_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: sfc_wts(:,:,:)
real(kind=kind_dbl_prec), intent(inout) :: spp_wts(:,:,:,:)
integer,                  intent(in)    :: nthreads

real,allocatable :: tmp_wts(:,:),tmpu_wts(:,:,:),tmpv_wts(:,:,:),tmpl_wts(:,:,:),tmp_spp_wts(:,:,:)
!D-grid
integer :: k,v
integer j,ierr,i
integer :: nblks, blk, len, maxlen
character*120 :: sfile
character*6   :: STRFH
logical :: do_advance_pattern

if ( (.NOT. do_sppt) .AND. (.NOT. do_shum) .AND. (.NOT. do_skeb) .AND. (lndp_type==0 ) .AND. (n_var_spp .le. 0)) return

! Update number of threads in shared variables in spectral_layout_mod and set block-related variables
nblks = size(blksz)
maxlen = maxval(blksz(:))


if ( (lndp_type==1) .and. (kdt==0) ) then ! old land pert scheme called once at start
        write(0,*) 'calling get_random_pattern_sfc'  
        allocate(tmpl_wts(nblks,maxlen,n_var_lndp))
        call get_random_pattern_sfc(rpattern_sfc,nlndp,gis_stochy,tmpl_wts)
        DO blk=1,nblks
           len=blksz(blk)
           ! for perturbing vars or states, saved value is N(0,1)  and apply scaling later.
           DO k=1,n_var_lndp
               !sfc_wts(blk,1:len,k) = tmpl_wts(blk,1:len,k)
               sfc_wts(blk,1:len,k) = tmpl_wts(1:len,blk,k)
           ENDDO
        ENDDO
        deallocate(tmpl_wts)
endif
allocate(tmp_wts(gis_stochy%nx,gis_stochy%ny))
allocate(tmpu_wts(gis_stochy%nx,gis_stochy%ny,levs))
allocate(tmpv_wts(gis_stochy%nx,gis_stochy%ny,levs))
if (do_sppt) then
   if (mod(kdt,nssppt) == 1 .or. nssppt == 1) then
      call get_random_pattern_scalar(rpattern_sppt,nsppt,gis_stochy,tmp_wts)
      DO blk=1,nblks
         len=blksz(blk)
         DO k=1,levs
            !sppt_wts(blk,1:len,k)=tmp_wts(blk,1:len)*vfact_sppt(k)
            sppt_wts(blk,1:len,k)=tmp_wts(1:len,blk)*vfact_sppt(k)
         ENDDO
         if (sppt_logit) sppt_wts(blk,:,:) = (2./(1.+exp(sppt_wts(blk,:,:))))-1.
         sppt_wts(blk,:,:) = sppt_wts(blk,:,:)+1.0
      ENDDO
   endif
endif
if (do_shum) then
   if (mod(kdt,nsshum) == 1 .or. nsshum == 1) then
      call get_random_pattern_scalar(rpattern_shum,nshum,gis_stochy,tmp_wts)
      DO blk=1,nblks
         len=blksz(blk)
         DO k=1,levs
            shum_wts(blk,1:len,k)=tmp_wts(1:len,blk)*vfact_shum(k)
         ENDDO
      ENDDO
   endif
endif
if (do_skeb) then
   if (mod(kdt,nsskeb) == 1 .or. nsskeb == 1) then
      call get_random_pattern_vector(rpattern_skeb,nskeb,gis_stochy,tmpu_wts,tmpv_wts)
      DO blk=1,nblks
         len=blksz(blk)
         DO k=1,levs
            skebu_wts(blk,1:len,k)=tmpu_wts(1:len,blk,k)*vfact_skeb(k)
            skebv_wts(blk,1:len,k)=tmpv_wts(1:len,blk,k)*vfact_skeb(k)
         ENDDO
      ENDDO
   endif
endif
if ( lndp_type .EQ. 2  ) then 
    ! add time check?
    allocate(tmpl_wts(gis_stochy%nx,gis_stochy%ny,n_var_lndp))
    call get_random_pattern_sfc(rpattern_sfc,nlndp,gis_stochy,tmpl_wts)
    DO blk=1,nblks
       len=blksz(blk)
       ! for perturbing vars or states, saved value is N(0,1)  and apply scaling later.
       DO k=1,n_var_lndp
           sfc_wts(blk,1:len,k) = tmpl_wts(1:len,blk,k)
       ENDDO
    ENDDO
    deallocate(tmpl_wts)
endif
if (n_var_spp .GE. 1) then
    allocate(tmp_spp_wts(gis_stochy%nx,gis_stochy%ny,n_var_spp))
    call get_random_pattern_spp(rpattern_spp,nspp,gis_stochy,tmp_spp_wts)
     DO v=1,n_var_spp
       DO blk=1,nblks
         len=blksz(blk)
         DO k=1,levs
           if (spp_stddev_cutoff(v).gt.0.0) then
             spp_wts(blk,1:len,k,v)=MAX(MIN(tmp_spp_wts(1:len,blk,v)*vfact_spp(k),spp_stddev_cutoff(v)),-1.0*spp_stddev_cutoff(v))*spp_prt_list(v)
           else
             spp_wts(blk,1:len,k,v)=tmp_spp_wts(1:len,blk,v)*vfact_spp(k)*spp_prt_list(v)
           endif
         ENDDO
       ENDDO
     ENDDO
    deallocate(tmp_spp_wts)
endif
 deallocate(tmp_wts)
 deallocate(tmpu_wts)
 deallocate(tmpv_wts)


end subroutine run_stochastic_physics

subroutine run_stochastic_physics_ocn(sppt_wts,t_rp1,t_rp2)
!use MOM_forcing_type, only : mech_forcing
!use MOM_grid, only : ocean_grid_type   
use stochy_internal_state_mod
use stochy_data_mod, only : nepbl,nocnsppt,rpattern_epbl1,rpattern_epbl2,rpattern_ocnsppt, gis_stochy_ocn
use get_stochy_pattern_mod,only : get_random_pattern_scalar
use stochy_namelist_def
implicit none
!type(ocean_grid_type),       intent(in) :: G
real, intent(inout) :: sppt_wts(:,:),t_rp1(:,:),t_rp2(:,:)
real, allocatable :: tmp_wts(:,:)
if (pert_epbl .OR. do_ocnsppt) then
   allocate(tmp_wts(gis_stochy_ocn%nx,gis_stochy_ocn%ny))
   if (pert_epbl) then
      call get_random_pattern_scalar(rpattern_epbl1,nepbl,gis_stochy_ocn,tmp_wts)
      t_rp1(:,:)=2.0/(1+exp(-1*tmp_wts))
      call get_random_pattern_scalar(rpattern_epbl2,nepbl,gis_stochy_ocn,tmp_wts)
      t_rp2(:,:)=2.0/(1+exp(-1*tmp_wts))
   else
      t_rp1(:,:)=1.0
      t_rp2(:,:)=1.0
   endif
   if (do_ocnsppt) then
      call get_random_pattern_scalar(rpattern_ocnsppt,nocnsppt,gis_stochy_ocn,tmp_wts)
      sppt_wts=2.0/(1+exp(-1*tmp_wts))
   else
      sppt_wts=1.0
   endif
   deallocate(tmp_wts)
else
   sppt_wts(:,:)=1.0
   t_rp1(:,:)=1.0
   t_rp2(:,:)=1.0
endif

end subroutine run_stochastic_physics_ocn
subroutine finalize_stochastic_physics()
use stochy_data_mod, only : nshum,rpattern_shum,rpattern_sppt,nsppt,rpattern_skeb,nskeb,&
                            vfact_sppt,vfact_shum,vfact_skeb, skeb_vwts,skeb_vpts, &
                            rpattern_spp, vfact_spp, nspp, &
                            rpattern_sfc, nlndp,gg_lats,gg_lons,sl,skebu_save,skebv_save,gis_stochy
use spectral_transforms, only : lat1s_a ,lon_dims_a,wgt_a,sinlat_a,coslat_a,colrad_a,rcs2_a
implicit none

   if (allocated(gg_lats)) deallocate (gg_lats)
   if (allocated(gg_lons)) deallocate (gg_lons)
   if (allocated(sl)) deallocate (sl)
   if (nsppt > 0) then 
      if (allocated(rpattern_sppt)) deallocate(rpattern_sppt)
      if (allocated(vfact_sppt)) deallocate(vfact_sppt)
   endif
   if (nshum > 0) then
      if (allocated(rpattern_shum)) deallocate(rpattern_shum)
      if (allocated(vfact_shum)) deallocate(vfact_shum)
   endif
   if (nskeb > 0) then
      if (allocated(rpattern_skeb)) deallocate(rpattern_skeb)
      if (allocated(skeb_vwts)) deallocate (skeb_vwts)
      if (allocated(skeb_vpts)) deallocate (skeb_vpts)
      if (allocated(skebu_save)) deallocate (skebu_save)
      if (allocated(skebv_save)) deallocate (skebv_save)
      if (allocated(vfact_skeb)) deallocate(vfact_skeb)
   endif
   if (nlndp > 0) then
      if (allocated(rpattern_sfc)) deallocate(rpattern_sfc)
   endif
   if (nspp > 0) then 
      if (allocated(rpattern_spp)) deallocate(rpattern_spp)
      if (allocated(vfact_spp)) deallocate(vfact_spp)
   endif

deallocate(lat1s_a)
deallocate(lon_dims_a)
deallocate(wgt_a)
deallocate(rcs2_a)
deallocate(colrad_a)
deallocate(sinlat_a)
deallocate(coslat_a)
deallocate(gis_stochy%ls_node)
deallocate(gis_stochy%ls_nodes)
deallocate(gis_stochy%max_ls_nodes)
deallocate(gis_stochy%lats_nodes_a)
deallocate(gis_stochy%global_lats_a)
deallocate(gis_stochy%epse)
deallocate(gis_stochy%epso)
deallocate(gis_stochy%epsedn)
deallocate(gis_stochy%epsodn)
deallocate(gis_stochy%kenorm_e)
deallocate(gis_stochy%kenorm_o)
deallocate(gis_stochy%snnp1ev)
deallocate(gis_stochy%snnp1od)
deallocate(gis_stochy%plnev_a)
deallocate(gis_stochy%plnod_a)
deallocate(gis_stochy%plnew_a)
deallocate(gis_stochy%plnow_a)

end subroutine finalize_stochastic_physics

end module stochastic_physics
