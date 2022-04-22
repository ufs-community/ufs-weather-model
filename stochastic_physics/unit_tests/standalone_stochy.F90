program  standalone_stochy

use stochastic_physics,  only : init_stochastic_physics,run_stochastic_physics
use get_stochy_pattern_mod,  only : write_stoch_restart_atm

use atmosphere_stub_mod, only: Atm,atmosphere_init_stub
!use mpp_domains
use mpp_mod,             only: mpp_set_current_pelist,mpp_get_current_pelist,mpp_init,mpp_pe,mpp_npes ,mpp_declare_pelist,mpp_root_pe
use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
use fms_mod,             only:  fms_init
use xgrid_mod,           only: grid_box_type
use netcdf
use kinddef,             only : kind_dbl_prec,kind_phys
use stochy_namelist_def, only : stochini

implicit none
integer, parameter      :: nlevs=3
integer, parameter :: max_n_var_lndp = 6
integer                 :: ntasks,fid
integer                 :: nthreads
integer                 :: ncid,xt_dim_id,yt_dim_id,time_dim_id,xt_var_id,yt_var_id,time_var_id,var_id_lat,var_id_lon,var_id_tile
integer                 :: varid1,varid2,varid3,varid4,varid_lon,varid_lat,varid_tile
integer                 :: varidl(max_n_var_lndp)
integer                 :: zt_dim_id,zt_var_id
character*2             :: strid

character(len=3), dimension(max_n_var_lndp)         ::  lndp_var_list
real(kind=kind_dbl_prec), dimension(max_n_var_lndp) ::  lndp_prt_list
include 'mpif.h'
include 'netcdf.inc'
real :: ak(nlevs+1),bk(nlevs+1)
real(kind=4) :: ts,undef

data ak(1:4) /0.0, 306.1489, 13687.72    , 0.99/
data bk(1:4) /1.0,   0.9284,     0.013348, 0.0/
integer     :: nb,blksz_1,nblks,ierr,my_id,i,j,k,l,nx,ny,id
integer     :: isc,iec,jsc,jec,isd,ied,jsd,jed
integer :: halo_update_type = 1
real        :: dx,dy,pi,rd,cp
logical   :: write_this_tile
integer  :: nargs,ntile_out,nlunit,pe,npes,stackmax=4000000
integer  :: i1,i2,j1,npts,istart,tpt
character*80 :: fname
character*1  :: ntile_out_str
integer :: comm

real(kind=4),allocatable,dimension(:,:) :: workg,tile_number
real(kind=4),allocatable,dimension(:,:,:) :: workg3d
real(kind=4),allocatable,dimension(:) :: grid_xt,grid_yt
real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlat
real(kind=kind_phys), dimension(:,:),   allocatable, save :: xlon
real(kind=kind_dbl_prec)                    :: ex3d(nlevs+1),pressi(nlevs+1),pressl(nlevs),p1000,exn

type(grid_box_type)           :: grid_box
real (kind=kind_phys),allocatable :: shum_wts  (:,:,:)
real (kind=kind_phys),allocatable :: sppt_wts  (:,:,:)
real (kind=kind_phys),allocatable :: sppt_pattern(:,:)
real (kind=kind_phys),allocatable :: skebu_wts (:,:,:)
real (kind=kind_phys),allocatable :: skebv_wts (:,:,:)
real (kind=kind_phys),allocatable :: sfc_wts   (:,:,:)
integer,allocatable :: blksz(:)
integer              :: me              !< MPI rank designator
integer              :: root_pe         !< MPI rank of root atmosphere processor
real(kind=kind_phys) :: dtp             !< physics timestep in seconds
real(kind=kind_phys) :: fhour           !< previous forecast hour
real(kind=kind_phys) :: sppt_amp        !< amplitude of sppt (to go to cld scheme)
logical  :: do_sppt,do_shum,do_skeb,use_zmtnblck
integer  ::  skeb_npass,n_var_lndp, lndp_type
character(len=65) :: fn_nml                   !< namelist filename
character(len=256),allocatable :: input_nml_file(:) !< character string containing full namelist

      namelist /gfs_physics_nml/do_sppt,do_skeb,do_shum,lndp_type,n_var_lndp
write_this_tile=.false.
ntile_out_str='0'
nlunit=23
nargs=iargc()
if (nargs.EQ.1) then
   call getarg(1,ntile_out_str)
endif
read(ntile_out_str,'(I1.1)') ntile_out
open (unit=nlunit, file='input.nml', status='OLD')
n_var_lndp=0
lndp_type=0
do_sppt=.false.
do_shum=.false.
do_skeb=.false.
read(nlunit,gfs_physics_nml)
close(nlunit)
! define stuff
pi=3.14159265359
undef=9.99e+20
p1000=100000.0
!define mid-layer pressure
rd=287.0
cp=1004.0
DO k=1,nlevs
   pressi(k)=ak(k)+p1000*bk(k)
ENDDO
ex3d=cp*(pressi/p1000)**(rd/cp)
DO k=1,3 !nlevs
   exn = (ex3d(k)*pressi(k)-ex3d(k+1)*pressi(k+1))/((cp+rd)*(pressi(k)-pressi(k+1)))
   pressl(k)=p1000*exn**(cp/rd)
ENDDO
pressl(4:)=0.01

call fms_init()
call mpp_init()
call fms_init
my_id=mpp_pe()
ntasks=mpp_npes()

call atmosphere_init_stub (grid_box)
isd=Atm(1)%bd%isd
ied=Atm(1)%bd%ied
jsd=Atm(1)%bd%jsd
jed=Atm(1)%bd%jed
isc=Atm(1)%bd%isc
iec=Atm(1)%bd%iec
jsc=Atm(1)%bd%jsc
jec=Atm(1)%bd%jec
nx=iec-isc+1
ny=jec-jsc+1
allocate(workg(nx,ny))
allocate(tile_number(nx,ny))
allocate(workg3d(nx,ny,nlevs))
print*,'nx,ny=',nx,ny
blksz_1=nx
nblks=nx*ny/blksz_1
allocate(blksz(nblks))
do i=1,nblks
  blksz(i)=blksz_1
enddo
nthreads = 1
me=my_id
fhour=0
dtp=600
fn_nml='input.nml'
nlunit=21

!define model grid
dx=360.0/nx
dy=180.0/ny
allocate(xlat(nblks,blksz_1))
allocate(xlon(nblks,blksz_1))
i1=isc
j1=jsc
do nb=1,nblks
    i2=i1+blksz_1-1
    if (i2 .le. iec) then 
       xlon(nb,1:blksz_1) = Atm(1)%gridstruct%agrid_64(i1:i2,j1,1)
       xlat(nb,1:blksz_1) = Atm(1)%gridstruct%agrid_64(i1:i2,j1,2)
       i1=i1+blksz_1
    else
       npts=iec-i1+1
       xlon(nb,1:npts) = Atm(1)%gridstruct%agrid_64(i1:iec,j1,1)
       xlat(nb,1:npts) = Atm(1)%gridstruct%agrid_64(i1:iec,j1,2)
       if (j1.LT. jec) then
          xlon(nb,npts+1:blksz_1) = Atm(1)%gridstruct%agrid_64(isc:isc+(blksz_1-npts+1),j1+1,1)
          xlat(nb,npts+1:blksz_1) = Atm(1)%gridstruct%agrid_64(isc:isc+(blksz_1-npts+1),j1+1,2)
       endif
       i1=npts+1
       j1=j1+1
    endif
    if (i2.EQ.iec) then
       i1=isc
       j1=j1+1
    endif
end do

allocate(grid_xt(nx),grid_yt(ny))
do i=1,nx
  grid_xt(i)=i
enddo
do j=1,ny
  grid_yt(j)=j
enddo
print*,'calling init_stochastic_physics',nlevs
root_pe=mpp_root_pe()
allocate(input_nml_file(1))
input_nml_file='input.nml'
comm=MPI_COMM_WORLD
call init_stochastic_physics(nlevs, blksz, dtp, sppt_amp,                         &
     input_nml_file, fn_nml, nlunit, xlon, xlat, do_sppt, do_shum,                &
     do_skeb, lndp_type, n_var_lndp, use_zmtnblck, skeb_npass, &
     lndp_var_list, lndp_prt_list,    &
     ak, bk, nthreads, root_pe, comm, ierr)
if (ierr .ne. 0) print *, 'ERROR init_stochastic_physics call' ! Draper - need proper error trapping here
call get_outfile(fname)
write(strid,'(I2.2)') my_id+1
if (ntile_out.EQ.0) write_this_tile=.true.
if ((my_id+1).EQ.ntile_out) write_this_tile=.true.
print*,trim(fname)//'.tile'//strid//'.nc',write_this_tile
if (write_this_tile) then
   fid=30+my_id
   ierr=nf90_create(trim(fname)//'.tile'//strid//'.nc',cmode=NF90_CLOBBER,ncid=ncid)
   ierr=NF90_DEF_DIM(ncid,"grid_xt",nx,xt_dim_id)
   ierr=NF90_DEF_DIM(ncid,"grid_yt",ny,yt_dim_id)
   if (do_skeb)ierr=NF90_DEF_DIM(ncid,"p_ref",nlevs,zt_dim_id)
   ierr=NF90_DEF_DIM(ncid,"time",NF90_UNLIMITED,time_dim_id)
  !> - Define the dimension variables.
   ierr=NF90_DEF_VAR(ncid,"grid_xt",NF90_FLOAT,(/ xt_dim_id /), xt_var_id)
   ierr=NF90_PUT_ATT(ncid,xt_var_id,"long_name","T-cell longitude")
   ierr=NF90_PUT_ATT(ncid,xt_var_id,"cartesian_axis","X")
   ierr=NF90_PUT_ATT(ncid,xt_var_id,"units","degrees_E")
   ierr=NF90_DEF_VAR(ncid,"grid_yt",NF90_FLOAT,(/ yt_dim_id /), yt_var_id)
   ierr=NF90_PUT_ATT(ncid,yt_var_id,"long_name","T-cell latitude")
   ierr=NF90_PUT_ATT(ncid,yt_var_id,"cartesian_axis","Y")
   ierr=NF90_PUT_ATT(ncid,yt_var_id,"units","degrees_N")
   ierr=NF90_DEF_VAR(ncid,"grid_lat",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_lat)
   ierr=NF90_PUT_ATT(ncid,var_id_lat,"long_name","T-cell latitudes")
   ierr=NF90_PUT_ATT(ncid,var_id_lat,"units","degrees_N")
   ierr=NF90_PUT_ATT(ncid,var_id_lat,"missing_value",undef)
   ierr=NF90_PUT_ATT(ncid,var_id_lat,"_FillValue",undef)
   ierr=NF90_DEF_VAR(ncid,"grid_lon",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_lon)
   ierr=NF90_PUT_ATT(ncid,var_id_lon,"long_name","T-cell longitudes")
   ierr=NF90_PUT_ATT(ncid,var_id_lon,"units","degrees_N")
   ierr=NF90_PUT_ATT(ncid,var_id_lon,"missing_value",undef)
   ierr=NF90_PUT_ATT(ncid,var_id_lon,"_FillValue",undef)
   ierr=NF90_DEF_VAR(ncid,"tile_num",NF90_FLOAT,(/ xt_dim_id, yt_dim_id, time_dim_id /), var_id_tile)
   ierr=NF90_PUT_ATT(ncid,var_id_tile,"long_name","tile number")
   ierr=NF90_PUT_ATT(ncid,var_id_tile,"missing_value",undef)
   ierr=NF90_PUT_ATT(ncid,var_id_tile,"_FillValue",undef)
   if (do_skeb)then
      ierr=NF90_DEF_VAR(ncid,"p_ref",NF90_FLOAT,(/ zt_dim_id /), zt_var_id)
      ierr=NF90_PUT_ATT(ncid,zt_var_id,"long_name","reference pressure")
      ierr=NF90_PUT_ATT(ncid,zt_var_id,"cartesian_axis","Z")
      ierr=NF90_PUT_ATT(ncid,zt_var_id,"units","Pa")
   endif
   ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
   ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
   ierr=NF90_PUT_ATT(ncid,time_var_id,"long_name","time")
   ierr=NF90_PUT_ATT(ncid,time_var_id,"units","hours since 2014-08-01 00:00:00")
   ierr=NF90_PUT_ATT(ncid,time_var_id,"cartesian_axis","T")
   ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar_type","JULIAN")
   ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar","JULIAN")
   if (do_sppt)then
      ierr=NF90_DEF_VAR(ncid,"sppt_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varid1)
      ierr=NF90_PUT_ATT(ncid,varid1,"long_name","sppt pattern")
      ierr=NF90_PUT_ATT(ncid,varid1,"units","None")
      ierr=NF90_PUT_ATT(ncid,varid1,"missing_value",undef)
      ierr=NF90_PUT_ATT(ncid,varid1,"_FillValue",undef)
      ierr=NF90_PUT_ATT(ncid,varid1,"cell_methods","time: point")
   endif
   if (do_shum)then
      ierr=NF90_DEF_VAR(ncid,"shum_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varid2)
      ierr=NF90_PUT_ATT(ncid,varid2,"long_name","shum pattern")
      ierr=NF90_PUT_ATT(ncid,varid2,"units","None")
      ierr=NF90_PUT_ATT(ncid,varid2,"missing_value",undef)
      ierr=NF90_PUT_ATT(ncid,varid2,"_FillValue",undef)
      ierr=NF90_PUT_ATT(ncid,varid2,"cell_methods","time: point")
   endif
   if (do_skeb)then
      ierr=NF90_DEF_VAR(ncid,"skebu_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,zt_dim_id,time_dim_id/), varid3)
      ierr=NF90_DEF_VAR(ncid,"skebv_wts",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,zt_dim_id,time_dim_id/), varid4)
      ierr=NF90_PUT_ATT(ncid,varid3,"long_name","skeb u pattern")
      ierr=NF90_PUT_ATT(ncid,varid3,"units","None")
      ierr=NF90_PUT_ATT(ncid,varid3,"missing_value",undef)
      ierr=NF90_PUT_ATT(ncid,varid3,"_FillValue",undef)
      ierr=NF90_PUT_ATT(ncid,varid3,"cell_methods","time: point")
      ierr=NF90_PUT_ATT(ncid,varid4,"long_name","skeb v pattern")
      ierr=NF90_PUT_ATT(ncid,varid4,"units","None")
      ierr=NF90_PUT_ATT(ncid,varid4,"missing_value",undef)
      ierr=NF90_PUT_ATT(ncid,varid4,"_FillValue",undef)
      ierr=NF90_PUT_ATT(ncid,varid4,"cell_methods","time: point")
   endif
   if (lndp_type > 0)then
      do l=1,n_var_lndp
         ierr=NF90_DEF_VAR(ncid,lndp_var_list(l),NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), varidl(l))
         ierr=NF90_PUT_ATT(ncid,varidl(l),"long_name",lndp_var_list(l)//" pattern")
         ierr=NF90_PUT_ATT(ncid,varidl(l),"units","None")
         ierr=NF90_PUT_ATT(ncid,varidl(l),"missing_value",undef)
         ierr=NF90_PUT_ATT(ncid,varidl(l),"_FillValue",undef)
         ierr=NF90_PUT_ATT(ncid,varidl(l),"cell_methods","time: point")
      enddo
   endif
   ierr=NF90_ENDDEF(ncid)
   ierr=NF90_PUT_VAR(ncid,xt_var_id,grid_xt)
   ierr=NF90_PUT_VAR(ncid,yt_var_id,grid_yt)
   if (do_skeb)then
      ierr=NF90_PUT_VAR(ncid,zt_var_id,pressl)
   endif
endif
! put lat lon and tile number
!ierr=NF90_PUT_VAR(ncid,var_id_lon,transpose(xlon(isc:iec,jsc:iec)),(/1,1,1/))
!ierr=NF90_PUT_VAR(ncid,var_id_lat,transpose(xlat(isc:iec,jsc:iec)),(/1,1,1/))
ierr=NF90_PUT_VAR(ncid,var_id_lon,transpose(xlon(:,:)),(/1,1,1/))
ierr=NF90_PUT_VAR(ncid,var_id_lat,transpose(xlat(:,:)),(/1,1,1/))
tile_number=my_id+1
ierr=NF90_PUT_VAR(ncid,var_id_tile,tile_number,(/1,1,1/))
if (do_sppt)allocate(sppt_wts(nblks,blksz_1,nlevs))
if (do_shum)allocate(shum_wts(nblks,blksz_1,nlevs))
if (do_skeb)allocate(skebu_wts(nblks,blksz_1,nlevs))
if (do_skeb)allocate(skebv_wts(nblks,blksz_1,nlevs))
if (lndp_type > 0)allocate(sfc_wts(nblks,blksz_1,n_var_lndp))
if (stochini) then
   istart=11
else
   istart=1
endif
tpt=1
do i=istart,21
   ts=i/4.0
   call run_stochastic_physics(nlevs, i-1, fhour, blksz, &
                               sppt_wts=sppt_wts, shum_wts=shum_wts, skebu_wts=skebu_wts, skebv_wts=skebv_wts, sfc_wts=sfc_wts, &
                               nthreads=nthreads)
  
   if (me.EQ.0 .and. do_sppt) print*,'SPPT_WTS=',i,sppt_wts(1,1,2)
   if (i.EQ. 10) call write_stoch_restart_atm('stochy_middle.nc')
   if (i.eq.1 .OR. i.eq.20) then
      if (me.EQ.0 .and. do_sppt) print*,'writing sppt_wts=',i,sppt_wts(1,1,2)
      if (write_this_tile) then
      if (do_sppt)then
         do j=1,ny
            workg(:,j)=sppt_wts(j,:,2)   
         enddo
         ierr=NF90_PUT_VAR(ncid,varid1,workg,(/1,1,tpt/))
      endif
      if (do_shum)then
         do j=1,ny
            workg(:,j)=shum_wts(j,:,1)
         enddo
         ierr=NF90_PUT_VAR(ncid,varid2,workg,(/1,1,tpt/))
      endif
      if (do_skeb)then
         do k=1,nlevs
            do j=1,ny
               workg3d(:,j,k)=skebu_wts(j,:,k)
            enddo
         enddo
         ierr=NF90_PUT_VAR(ncid,varid3,workg3d,(/1,1,1,tpt/))
         do k=1,nlevs
            do j=1,ny
               workg3d(:,j,k)=skebv_wts(j,:,k)
            enddo
         enddo
         ierr=NF90_PUT_VAR(ncid,varid4,workg3d,(/1,1,1,tpt/))
      endif
      if (lndp_type > 0)then
         do l=1,n_var_lndp
            do j=1,ny
               workg(:,j)=sfc_wts(j,:,l)
            enddo
            ierr=NF90_PUT_VAR(ncid,varidl(l),workg,(/1,1,tpt/))
         enddo
      endif
      ierr=NF90_PUT_VAR(ncid,time_var_id,ts,(/tpt/))
      endif
      tpt=tpt+1
   endif
enddo
if (write_this_tile) ierr=NF90_CLOSE(ncid)
if (stochini) then
   call write_stoch_restart_atm('stochy_final_2.nc')
else
   call write_stoch_restart_atm('stochy_final.nc')
endif
end
subroutine get_outfile(fname)
use stochy_namelist_def
character*80,intent(out) :: fname
character*4   :: s_ntrunc,s_lat,s_lon
   write(s_ntrunc,'(I4)') ntrunc
   write(s_lat,'(I4)') lat_s 
   write(s_lon,'(I4)') lon_s  
   fname=trim('workg_T'//trim(adjustl(s_ntrunc))//'_'//trim(adjustl(s_lon))//'x'//trim(adjustl(s_lat)))
   return
end
