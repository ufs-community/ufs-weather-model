program compare_output


use netcdf

implicit none

integer, parameter :: CRES=96
integer, parameter :: ntiles=6
integer, parameter :: npz=3
integer, parameter :: nexpt=3
integer, allocatable :: layout_x(:),layout_y(:)

real, allocatable :: truth(:,:,:),expt(:,:,:),diff(:,:,:),accum_error(:)
real, allocatable :: truth_3d(:,:,:,:),expt_3d(:,:,:,:),diff_3d(:,:,:,:)

integer :: i1,i2,j1,j2,nx,ny,ierr1,ierr2,ierr3,i,j,k,ncid,varid,ncid2,t,t2,var

character*240 :: baseline_path
character*2   :: tile1,tile2
character*1   :: lx,ly
character*10  :: var_list_2d(4)
character*10  :: var_list_3d(2)

baseline_path='/scratch2/BMC/gsienkf/Philip.Pegion/stochastic_physics_unit_tests/baseline_20210806/'

var_list_2d(1)='sppt_wts'
var_list_2d(2)='shum_wts'
var_list_2d(3)='vgf'
var_list_2d(4)='smc'
var_list_3d(1)='skebu_wts'
var_list_3d(2)='skebv_wts'

allocate(layout_x(nexpt))
allocate(layout_y(nexpt))
allocate(accum_error(nexpt))
allocate(truth(CRES,CRES,2))
allocate(truth_3d(CRES,CRES,npz,2))
layout_x(1)=1
layout_x(2)=2
layout_x(3)=4
layout_y(1)=4
layout_y(2)=2
layout_y(3)=1

accum_error(:)=0.0
! loop through 2-d vars (sppt,shum,landp)
do k=1,nexpt
   nx=CRES/layout_x(k)
   ny=CRES/layout_y(k)
   allocate(expt(nx,ny,2))
   allocate(diff(nx,ny,2))
   do var=1,4 ! 2d-vars
      t2=1
      do t=1,ntiles
         write(tile1,fmt='(I2.2)') t
         ierr1=nf90_open(trim(baseline_path)//'workg_T162_984x488.tile'//tile1//'.nc',mode=nf90_nowrite,ncid=ncid)
         ierr2=nf90_inq_varid(ncid,trim(var_list_2d(var)),varid)
         ierr3=nf90_get_var(ncid,varid,truth,count=(/CRES,CRES,2/))
         if (ierr1+ierr2+ierr3.NE.0) then
            print*,'error reading in truth files'
            call exit(1)
         endif
         i1=1
         i2=i1+nx-1
         j1=1
         j2=j1+ny-1
         do j=1,layout_y(k)
            do i=1,layout_x(k)
               write(tile2,fmt='(I2.2)') t2
               write(lx,fmt='(I1)') layout_x(k)
               write(ly,fmt='(I1)') layout_y(k)
               ierr1=nf90_open('layout_'//lx//'x'//ly//'/workg_T162_984x488.tile'//tile2//'.nc',mode=nf90_nowrite,ncid=ncid2)
               ierr2=nf90_inq_varid(ncid2,trim(var_list_2d(var)),varid)
               ierr3=nf90_get_var(ncid2,varid,expt,count=(/nx,ny,2/))
               if (ierr1+ierr2+ierr3.NE.0) then
                  print*,'error reading in expt files',ierr1,ierr2,ierr3
                  call exit(2)
               endif
               diff(:,:,:)=expt(:,:,:)-truth(i1:i2,j1:j2,:)
               accum_error(k)=accum_error(k)+sum(abs(diff))
               i1=i1+nx
               i2=i2+nx
               if (i2.GT.CRES) then
                  i1=1
                  i2=i1+nx-1
                  j1=j1+ny
                  j2=j2+ny
               endif
               ierr1=nf90_close(ncid2)
               t2=t2+1
            enddo
         enddo
         ierr1=nf90_close(ncid)
      enddo
   enddo
   deallocate(expt)
   deallocate(diff)
enddo

! loop through 3-d vars (sppt,shum,landp)
do k=1,nexpt
   nx=CRES/layout_x(k)
   ny=CRES/layout_y(k)
   allocate(expt_3d(nx,ny,npz,2))
   allocate(diff_3d(nx,ny,npz,2))
   do var=1,2
      t2=1
      do t=1,ntiles
         write(tile1,fmt='(I2.2)') t
         ierr1=nf90_open(trim(baseline_path)//'workg_T162_984x488.tile'//tile1//'.nc',mode=nf90_nowrite,ncid=ncid)
         ierr2=nf90_inq_varid(ncid,trim(var_list_3d(var)),varid)
         ierr3=nf90_get_var(ncid,varid,truth_3d,count=(/CRES,CRES,npz,2/))
         if (ierr1+ierr2+ierr3.NE.0) then
            print*,'error reading in truth files'
            call exit(1)
         endif
         i1=1
         i2=i1+nx-1
         j1=1
         j2=j1+ny-1
         do j=1,layout_y(k)
            do i=1,layout_x(k)
               write(tile2,fmt='(I2.2)') t2
               write(lx,fmt='(I1)') layout_x(k)
               write(ly,fmt='(I1)') layout_y(k)
               ierr1=nf90_open('layout_'//lx//'x'//ly//'/workg_T162_984x488.tile'//tile2//'.nc',mode=nf90_nowrite,ncid=ncid2)
               ierr2=nf90_inq_varid(ncid2,trim(var_list_3d(var)),varid)
               ierr3=nf90_get_var(ncid2,varid,expt_3d,count=(/nx,ny,npz,2/))
               if (ierr1+ierr2+ierr3.NE.0) then
                  print*,'error reading in expt files',ierr1,ierr2,ierr3
                  call exit(2)
               endif
               diff_3d(:,:,:,:)=expt_3d(:,:,:,:)-truth_3d(i1:i2,j1:j2,:,:)
               accum_error(k)=accum_error(k)+sum(abs(diff_3d))
               i1=i1+nx
               i2=i2+nx
               if (i2.GT.CRES) then
                  i1=1
                  i2=i1+nx-1
                  j1=j1+ny
                  j2=j2+ny
               endif
               ierr1=nf90_close(ncid2)
               t2=t2+1
            enddo
         enddo
         ierr1=nf90_close(ncid)
      enddo
   enddo
   deallocate(expt_3d)
   deallocate(diff_3d)
enddo

if (sum(accum_error(:)) .EQ. 0.0 )then
   print*,'all tests pass'
else
   print*,'decomposition test fail',accum_error
   call exit(3)
endif

end
