program  standalone_ca_global

use cellular_automata_global_mod, only : cellular_automata_global
use cellular_automata_sgs_mod, only : cellular_automata_sgs
use update_ca, only : write_ca_restart,read_ca_restart
use atmosphere_stub_mod, only: Atm,atmosphere_init_stub
!use mpp_domains
use mpp_mod,             only: mpp_set_current_pelist,mpp_get_current_pelist,mpp_init,mpp_pe,mpp_npes ,mpp_declare_pelist,mpp_root_pe
use mpp_domains_mod,     only: mpp_broadcast_domain,MPP_DOMAIN_TIME,mpp_domains_init ,mpp_domains_set_stack_size
use fms_mod,             only:  fms_init
!use time_manager_mod,    only: time_type
use xgrid_mod,           only: grid_box_type
use netcdf
use kinddef,             only : kind_dbl_prec,kind_phys


implicit none
integer                 :: ntasks,fid,ct,levs,ntiles
integer                 :: ncid_in,varid,ncid,xt_dim_id,yt_dim_id,time_dim_id,xt_var_id,yt_var_id,time_var_id,ca_out_id
integer                 :: ca1_id,ca2_id,ca3_id,ca_deep_id!,ca_turb_id,ca_shal_id
integer                 :: root_pe,comm,dump_time
real(kind=kind_phys)    :: dtf, nthresh
character*4             :: strid
character*1             :: tileid
character*4             :: CRES
!type(GFS_statein_type),allocatable :: Statein(:)
include 'mpif.h'
include 'netcdf.inc'
real(kind=4) :: ts,undef

integer     :: nblks,blksz,ierr,my_id,i,j,nx,ny,id,i1,i2
integer     :: isc,iec,jsc,jec,nb,npts
logical  :: first_time_step
integer  :: istart

real(kind=4),allocatable,dimension(:,:) :: workg
real(kind=4),allocatable,dimension(:) :: grid_xt,grid_yt
type(grid_box_type)           :: grid_box
!---cellular automata control parameters
integer              :: nca             !< number of independent cellular automata
integer              :: nlives          !< cellular automata lifetime
integer              :: ncells          !< cellular automata finer grid
integer              :: nca_g           !< number of independent cellular automata
integer              :: nlives_g        !< cellular automata lifetime
integer              :: ncells_g        !< cellular automata finer grid
real(kind=kind_phys) :: nfracseed       !< cellular automata seed probability
integer              :: nseed           !< cellular automata seed frequency
integer              :: nseed_g         !< cellular automata seed frequency
logical              :: do_ca           !< cellular automata main switch
logical              :: ca_sgs          !< switch for sgs ca
logical              :: ca_global       !< switch for global ca
logical              :: ca_smooth       !< switch for gaussian spatial filter
integer*8            :: iseed_ca        !< seed for random number generation in ca scheme
integer              :: nspinup         !< number of iterations to spin up the ca
real(kind=kind_phys) :: rcell           !< threshold used for CA scheme
real                 :: ca_amplitude    !< amplitude of ca trigger perturbation
integer              :: nsmooth         !< number of passes through smoother
logical              :: ca_closure      !< logical switch for ca on closure
logical              :: ca_entr         !< logical switch for ca on entrainment
logical              :: ca_trigger      !< logical switch for ca on trigger
logical              :: warm_start      !< logical switch for ca on trigger

real(kind=kind_phys), dimension(:,:),   allocatable :: cond_in,condition, sst,lmsk,lake
real(kind=kind_phys), dimension(:,:),   allocatable :: ca_deep, ca_turb, ca_shal

real(kind=kind_phys), dimension(:,:),   allocatable :: ca1, ca2, ca3

NAMELIST /gfs_physics_nml/ do_ca, ca_sgs, ca_global, nca, ncells, nlives, nseed,       &
                          nfracseed, rcell, ca_trigger, ca_entr, ca_closure, nca_g,    &
                          ncells_g, nlives_g, nseed_g, ca_smooth, nspinup, iseed_ca,   &
                          nsmooth, ca_amplitude, warm_start
! get mpi info,

first_time_step=.true.
warm_start=.false.
dtf=720/12.0
! default values
levs=63
nca            = 0
nca_g          = 0
ncells_g       = 1
nlives_g       = 1
nfracseed      = 0.5
nseed          = 100000
iseed_ca       = 0
nspinup        = 1
do_ca          = .false.
ca_sgs         = .false.
ca_global      = .false.
ca_smooth      = .false.
ca_amplitude   = 500.
rcell          = 0.0

! open namelist file
open (unit=565, file='input.nml', status='OLD', iostat=ierr)
read(565,gfs_physics_nml)
close(565)
! define stuff
undef=9.99e+20
print*,'ca_sgs,ca_global',ca_sgs,ca_global
if (.not. ca_sgs) then
   nca=0
endif
if (.not. ca_global) then
   nca_g=0
endif

! initialize fms
!call fms_init()
call mpp_init()
call fms_init
root_pe=mpp_root_pe()
comm=MPI_COMM_WORLD
my_id=mpp_pe()
ntasks=mpp_npes()
ntiles=6

call atmosphere_init_stub (grid_box)
!define domain
isc=Atm(1)%bd%isc
iec=Atm(1)%bd%iec
jsc=Atm(1)%bd%jsc
jec=Atm(1)%bd%jec
write(CRES,'(I4)') Atm(1)%npx-1
print*,'ATM npx,npy=',Atm(1)%npx,Atm(1)%npy

nx=iec-isc+1
ny=jec-jsc+1
allocate(workg(nx,ny))
print*,'after init',my_id,Atm(1)%tile_of_mosaic,isc,jec

! for this simple test, nblocks = ny, blksz=ny
blksz=nx
nblks=ny

! setup GFS_init parameters

!define model grid

allocate(grid_xt(nx),grid_yt(ny))
do i=1,nx
  grid_xt(i)=i
enddo
do i=1,ny
  grid_yt(i)=i
enddo

!setup GFS_coupling
if ( ntasks .GT. 1000) then
   write(strid,'(I4.4)') my_id+1
else if ( ntasks .GT. 100) then
   write(strid,'(I3.3)') my_id+1
else if ( ntasks .GT. 10) then
   write(strid,'(I2.2)') my_id+1
else
   write(strid,'(I1.1)') my_id+1
endif
fid=30+my_id
ierr=nf90_create('ca_out.tile'//trim(strid)//'.nc',cmode=NF90_CLOBBER,ncid=ncid)
ierr=NF90_DEF_DIM(ncid,"grid_xt",nx,xt_dim_id)
ierr=NF90_DEF_DIM(ncid,"grid_yt",ny,yt_dim_id)
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
ierr=NF90_DEF_VAR(ncid,"time",NF90_FLOAT,(/ time_dim_id /), time_var_id)
ierr=NF90_PUT_ATT(ncid,time_var_id,"long_name","time")
ierr=NF90_PUT_ATT(ncid,time_var_id,"units","hours since 2014-08-01 00:00:00")
ierr=NF90_PUT_ATT(ncid,time_var_id,"cartesian_axis","T")
ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar_type","JULIAN")
ierr=NF90_PUT_ATT(ncid,time_var_id,"calendar","JULIAN")
!ierr=NF90_DEF_VAR(ncid,"ca_out",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca_out_id)
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"long_name","random pattern")
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"units","None")
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"missing_value",undef)
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"_FillValue",undef)
!ierr=NF90_PUT_ATT(ncid,ca_out_id,"cell_methods","time: point")
if (ca_global) then
   ierr=NF90_DEF_VAR(ncid,"ca1",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca1_id)
   print*,'nc 1',ierr
   ierr=NF90_PUT_ATT(ncid,ca1_id,"long_name","random pattern")
   print*,'nc 2',ierr
   ierr=NF90_PUT_ATT(ncid,ca1_id,"units","None")
   print*,'nc 3',ierr
   ierr=NF90_PUT_ATT(ncid,ca1_id,"missing_value",undef)
   print*,'nc 4',ierr
   ierr=NF90_PUT_ATT(ncid,ca1_id,"_FillValue",undef)
   print*,'nc 5',ierr
   ierr=NF90_PUT_ATT(ncid,ca1_id,"cell_methods","time: point")
   print*,'nc 6',ierr
   ierr=NF90_DEF_VAR(ncid,"ca2",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca2_id)
   print*,'nc 7',ierr
   ierr=NF90_PUT_ATT(ncid,ca2_id,"long_name","random pattern")
   print*,'nc 8',ierr
   ierr=NF90_PUT_ATT(ncid,ca2_id,"units","None")
   print*,'nc 9',ierr
   ierr=NF90_PUT_ATT(ncid,ca2_id,"missing_value",undef)
   print*,'nc10',ierr
   ierr=NF90_PUT_ATT(ncid,ca2_id,"_FillValue",undef)
   print*,'nc11',ierr
   ierr=NF90_PUT_ATT(ncid,ca2_id,"cell_methods","time: point")
   print*,'nc12',ierr
   ierr=NF90_DEF_VAR(ncid,"ca3",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca3_id)
   print*,'nc13',ierr
   ierr=NF90_PUT_ATT(ncid,ca3_id,"long_name","random pattern")
   print*,'nc14',ierr
   ierr=NF90_PUT_ATT(ncid,ca3_id,"units","None")
   print*,'nc15',ierr
   ierr=NF90_PUT_ATT(ncid,ca3_id,"missing_value",undef)
   print*,'nc16',ierr
   ierr=NF90_PUT_ATT(ncid,ca3_id,"_FillValue",undef)
   print*,'nc17',ierr
   ierr=NF90_PUT_ATT(ncid,ca3_id,"cell_methods","time: point")
   print*,'nc18',ierr
endif
if (ca_sgs) then
   ierr=NF90_DEF_VAR(ncid,"ca_deep",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca_deep_id)
   print*,'ca_deep',ierr
   ierr=NF90_PUT_ATT(ncid,ca_deep_id,"long_name","CA field for deep convection")
   print*,'nc18',ierr
   ierr=NF90_PUT_ATT(ncid,ca_deep_id,"units","None")
   print*,'nc19',ierr
   ierr=NF90_PUT_ATT(ncid,ca_deep_id,"missing_value",undef)
   print*,'nc20',ierr
   ierr=NF90_PUT_ATT(ncid,ca_deep_id,"_FillValue",undef)
   print*,'nc21',ierr
   ierr=NF90_PUT_ATT(ncid,ca_deep_id,"cell_methods","time: point")
   print*,'nc22',ierr
   !ierr=NF90_DEF_VAR(ncid,"ca_turb",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca_turb_id)
   !ierr=NF90_PUT_ATT(ncid,ca_turb_id,"long_name","CA field for PBL")
   !ierr=NF90_PUT_ATT(ncid,ca_turb_id,"long_name","random pattern")
   !ierr=NF90_PUT_ATT(ncid,ca_turb_id,"units","None")
   !ierr=NF90_PUT_ATT(ncid,ca_turb_id,"missing_value",undef)
   !ierr=NF90_PUT_ATT(ncid,ca_turb_id,"_FillValue",undef)
   !ierr=NF90_PUT_ATT(ncid,ca_turb_id,"cell_methods","time: point")
   !ierr=NF90_DEF_VAR(ncid,"ca_shal",NF90_FLOAT,(/xt_dim_id, yt_dim_id ,time_dim_id/), ca_shal_id)
   !ierr=NF90_PUT_ATT(ncid,ca_shal_id,"long_name","CA field for shallow convection")
   !ierr=NF90_PUT_ATT(ncid,ca_shal_id,"units","None")
   !ierr=NF90_PUT_ATT(ncid,ca_shal_id,"missing_value",undef)
   !ierr=NF90_PUT_ATT(ncid,ca_shal_id,"_FillValue",undef)
   !ierr=NF90_PUT_ATT(ncid,ca_shal_id,"cell_methods","time: point")
endif
ierr=NF90_ENDDEF(ncid)
   print*,'nc23',ierr
ierr=NF90_PUT_VAR(ncid,xt_var_id,grid_xt)
   print*,'nc24',ierr
ierr=NF90_PUT_VAR(ncid,yt_var_id,grid_yt)
   print*,'nc25',ierr
! allocate diagnostics
if(ca_global)then
   allocate(ca1 (nblks,blksz))
   allocate(ca2 (nblks,blksz))
   allocate(ca3 (nblks,blksz))
endif

if(ca_sgs)then
   allocate(ca_deep (nblks,blksz))
   allocate(ca_turb (nblks,blksz))
   allocate(ca_shal (nblks,blksz))
   allocate(condition   (nblks,blksz))
   allocate(cond_in     (isc:iec,jsc:jec))
   allocate(sst         (nblks,blksz))
   allocate(lmsk        (nblks,blksz))
   allocate(lake        (nblks,blksz))
   sst(:,:)=303.
   lmsk(:,:)=0.
   lake(:,:)=0
! read in condtion
   write(tileid,'(I1)') Atm(1)%tile_of_mosaic
   ierr=NF90_OPEN('INPUT/C'//trim(adjustl(CRES))//'_ca_condition.tile'//tileid//'.nc',NF90_NOWRITE,ncid_in)
   if (ierr.NE.0) then
       print*,'error INPUT/C'//trim(adjustl(CRES))//'_ca_condition.tile'//tileid//'.nc'
       call MPI_ABORT(ierr)
   endif
   ierr=NF90_INQ_VARID(ncid_in,'ca_condition',varid)
   if (ierr.NE.0) then
       print*,'error gettinv varid for ca_condition'
       call MPI_ABORT(ierr)
   endif
   ierr=NF90_GET_VAR(ncid_in,varid,cond_in,start=(/isc,jsc,1/),count=(/nx,ny,1/))
   if (ierr.NE.0) then
       print*,'error getting var',isc,jsc,nx,ny
       call MPI_ABORT(ierr)
   endif
   ierr=NF90_CLOSE(ncid_in)
   
   i1=isc
   j=jsc
   do nb=1,nblks
       i2=i1+blksz-1
       if (i2 .le. iec) then  
          condition(nb,1:blksz) = cond_in(i1:i2,j)
          i1=i1+blksz
       else
          npts=iec-i1+1
          condition(nb,1:npts) = cond_in(i1:iec,j)
          if (j.LT. jec) then
             condition(nb,npts+1:blksz) = cond_in(isc:isc+(blksz-npts+1),j+1) 
          endif
          i1=npts+1
          j=j+1
       endif
       if (i2.EQ.iec) then
          i1=isc
          j=j+1
       endif
   end do
endif

dump_time=50
if (warm_start) then
   istart=dump_time+1
   call read_ca_restart(Atm(1)%domain,ncells,nca,ncells_g,nca_g)
else
   istart=1
endif
ct=1
do i=istart,101
   ts=i/4.0  ! hard coded to write out hourly based on a 900 second time-step
   if (ca_sgs) then
       call cellular_automata_sgs(i,dtf,warm_start,first_time_step,                            &
            sst,lmsk,lake,condition,ca_deep,ca_turb,ca_shal, &
            Atm(1)%domain_for_coupler,nblks,                                          &
            isc,iec,jsc,jec,Atm(1)%npx,Atm(1)%npy, levs,                                           &
            nthresh,Atm(1)%tile_of_mosaic,nca,ncells,nlives,nfracseed,                       & ! for new random number
            nseed,iseed_ca ,nspinup,ca_trigger,blksz,root_pe,comm)
   endif
   if (ca_global) then
      call cellular_automata_global(i,warm_start,first_time_step,ca1,ca2,ca3,Atm(1)%domain_for_coupler, &
           nblks,isc,iec,jsc,jec,Atm(1)%npx,Atm(1)%npy,levs,      &
           nca_g,ncells_g,nlives_g,nfracseed,nseed_g,                         &
           iseed_ca,Atm(1)%tile_of_mosaic, ca_smooth,nspinup,blksz,    &
           nsmooth,ca_amplitude,root_pe,comm)
   endif
   if (i.EQ. dump_time) call write_ca_restart('mid_run')
   first_time_step=.false.
   if (mod(i-1,5).eq.0) then
      if (ca_global) then
         workg(:,:)=TRANSPOSE(ca1(:,:))
         ierr=NF90_PUT_VAR(ncid,ca1_id,workg,(/1,1,ct/))
         print*,'put ca 1',ierr
         workg(:,:)=TRANSPOSE(ca2(:,:))
         ierr=NF90_PUT_VAR(ncid,ca2_id,workg,(/1,1,ct/))
         workg(:,:)=TRANSPOSE(ca3(:,:))
         ierr=NF90_PUT_VAR(ncid,ca3_id,workg,(/1,1,ct/))
      endif
      if (ca_sgs) then
         workg(:,:)=TRANSPOSE(ca_deep(:,:))
         ierr=NF90_PUT_VAR(ncid,ca_deep_id,workg,(/1,1,ct/))
         print*,'put ca_deep',ierr
         !workg(:,:)=TRANSPOSE(ca_turb(:,:))
         !ierr=NF90_PUT_VAR(ncid,ca_turb_id,workg,(/1,1,ct/))
         !workg(:,:)=ca_shal(:,:)   
         !workg(:,:)=cond_in(:,:)   
         !ierr=NF90_PUT_VAR(ncid,ca_shal_id,workg,(/1,1,ct/))
      endif
      ierr=NF90_PUT_VAR(ncid,time_var_id,ts,(/ct/))
      ct=ct+1
   endif
   if (ca_global) then
      if (my_id.EQ.0) write(6,fmt='(a,i7,f8.3)') 'ca glob =',i,maxval(ca1)
   endif
   if (ca_sgs) then
      if (my_id.EQ.0) write(6,fmt='(a,i7,f8.3)') 'ca sgs=',i,maxval(ca_deep)
   endif
enddo
call write_ca_restart()
!close(fid)
ierr=NF90_CLOSE(ncid)
end
