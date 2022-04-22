program compare_output


use netcdf

implicit none

integer, parameter :: CRES=96
integer, parameter :: ntiles=6
integer, parameter :: npz=3
integer, parameter :: nexpt=3
integer, allocatable :: layout_x(:),layout_y(:)

real, allocatable :: truth(:,:,:),expt(:,:,:),diff(:,:,:),accum_error(:)

integer :: i1,i2,j1,j2,nx,ny,ierr1,ierr2,ierr3,i,j,k,ncid,varid,ncid2,t,t2,var

character*240 :: baseline_path
character*1   :: tile1
character*2   :: tile2
character*1   :: lx,ly
character*10  :: var_list_2d(2)

baseline_path='/scratch2/BMC/gsienkf/Philip.Pegion/stochastic_physics_unit_tests/baseline_20210806/'

var_list_2d(1)='ca_deep'
var_list_2d(2)='ca1'

allocate(layout_x(nexpt))
allocate(layout_y(nexpt))
allocate(accum_error(nexpt))
allocate(truth(CRES,CRES,20))
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
   allocate(expt(nx,ny,20))
   allocate(diff(nx,ny,20))
   do var=1,2 ! 2d-vars
      t2=1
      do t=1,ntiles
         write(tile1,fmt='(I1.1)') t
         ierr1=nf90_open(trim(baseline_path)//'ca_out.tile'//tile1//'.nc',mode=nf90_nowrite,ncid=ncid)
         print*,k,var,t,ierr1
         ierr2=nf90_inq_varid(ncid,trim(var_list_2d(var)),varid)
         print*,'ierr2=',ierr2
         ierr3=nf90_get_var(ncid,varid,truth,count=(/CRES,CRES,20/))
         print*,'ierr3=',ierr3
         if (ierr1+ierr2+ierr3.NE.0) then
            print*,'error reading in truth files',ierr1,ierr2,ierr3
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
               ierr1=nf90_open('ca_layout_'//lx//'x'//ly//'/ca_out.tile'//tile2//'.nc',mode=nf90_nowrite,ncid=ncid2)
               print*,'opened','ca_layout_'//lx//'x'//ly//'/ca_out.tile'//tile2//'.nc'
               ierr2=nf90_inq_varid(ncid2,trim(var_list_2d(var)),varid)
               ierr3=nf90_get_var(ncid2,varid,expt,count=(/nx,ny,20/))
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

if (sum(accum_error(:)) .EQ. 0.0 )then
   print*,'all tests pass'
else
   print*,'decomposition test fail',accum_error
   call exit(3)
endif

end
