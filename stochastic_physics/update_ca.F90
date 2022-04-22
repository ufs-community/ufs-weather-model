module update_ca
!Main module for evolving CA in time, also includes
!read and write restart routines, to restart fields 
!on the ncellsxncells CA grid

use kinddef,           only: kind_dbl_prec
use halo_exchange,    only: atmosphere_scalar_field_halo
use random_numbers,   only: random_01_CB
use mpi_wrapper,      only: mype,mp_reduce_min,mp_reduce_max
use mpp_domains_mod,  only: domain2D,mpp_get_global_domain,CENTER, mpp_get_data_domain, mpp_get_compute_domain,mpp_get_ntile_count,&
                            mpp_define_mosaic,mpp_get_layout,mpp_define_io_domain,mpp_get_io_domain_layout
use mpp_mod,          only: mpp_error,  mpp_pe, mpp_root_pe, &
                            NOTE,   FATAL
use fms2_io_mod,        only: FmsNetcdfDomainFile_t, unlimited,      &
                                open_file, close_file,                 &
                                register_axis, register_restart_field, &
                                register_variable_attribute, register_field, &
                                read_restart, write_restart, write_data,     &
                                get_global_io_domain_indices, variable_exists


implicit none

public  write_ca_restart
public  read_ca_restart
public  update_cells_sgs
public  update_cells_global

integer,allocatable :: board(:,:,:), lives(:,:,:)
integer,allocatable :: board_g(:,:,:), lives_g(:,:,:)
integer,public :: isdnx,iednx,jsdnx,jednx,nxncells,nyncells
integer,public :: iscnx,iecnx,jscnx,jecnx,nxncells_g,nyncells_g
integer,public :: isdnx_g,iednx_g,jsdnx_g,jednx_g
integer,public :: iscnx_g,iecnx_g,jscnx_g,jecnx_g
integer*8, public :: csum
type(domain2D),public :: domain_sgs,domain_global
logical, public  :: cold_start_ca_sgs=.true.,cold_start_ca_global=.true.


contains

!Compute CA domain:--------------------------------------------------------------------------
subroutine define_ca_domain(domain_in,domain_out,ncells,nxncells_local,nyncells_local)
implicit none

type(domain2D),intent(inout) :: domain_in
type(domain2D),intent(inout) :: domain_out
integer,intent(in)           :: ncells
integer,intent(out) :: nxncells_local, nyncells_local
integer :: halo1 = 1
integer :: layout(2)
integer :: ntiles
integer, allocatable :: pe_start(:), pe_end(:)

integer :: i, j, k, n
integer :: nx, ny
integer :: isc,iec,jsc,jec

!--- get params from fv domain mosaic for building domain_out
  call mpp_get_global_domain(domain_in,xsize=nx,ysize=ny,position=CENTER)
  call mpp_get_layout(domain_in,layout)
  ntiles = mpp_get_ntile_count(domain_in)
  !write(1000+mpp_pe(),*) "nx,ny: ",nx,ny
  !write(1000+mpp_pe(),*) "layout: ",layout

!--- define mosaic for domain_out refined by 'ncells' from domain_in
  nxncells_local=nx*ncells+1
  nyncells_local=ny*ncells+1

  allocate(pe_start(ntiles))
  allocate(pe_end(ntiles))
  do n = 1, ntiles
    pe_start(n) = mpp_root_pe() + (n-1)*layout(1)*layout(2)
    pe_end(n)   = mpp_root_pe() +     n*layout(1)*layout(2)-1
  enddo
  call define_cubic_mosaic(domain_out, nxncells_local-1, nyncells_local-1, layout, pe_start, pe_end, halo1 )
  deallocate(pe_start)
  deallocate(pe_end)

end subroutine define_ca_domain
!---------------------------------------------------------------------------------------------

subroutine write_ca_restart(timestamp)
!Write restart files 


implicit none
character(len=*), optional, intent(in)    :: timestamp
character(len=32)  :: fn_ca = 'ca_data.nc'

type(FmsNetcdfDomainFile_t) :: CA_restart
integer :: id_restart,nx,ny,i
integer :: is,ie,js,je,nca,nca_g

integer, allocatable, dimension(:) :: buffer
character(7) :: indir='RESTART'
character(72) :: infile
logical :: amiopen
amiopen=.false.

!Return if not allocated:
if(.not. allocated(board) .and. .not. allocated(lives) .and. .not.  allocated(board_g) .and. .not. allocated(lives_g))return

infile=trim(indir)//'/'//trim(fn_ca)
if( present(timestamp) ) infile=trim(indir)//'/'//trim(timestamp)//'.'//trim(fn_ca)
   !--- register axis
if (allocated(board)) then
   amiopen=open_file(CA_restart, trim(infile), 'overwrite', domain=domain_sgs, is_restart=.true., dont_add_res_to_filename=.true.)
   if( amiopen ) then
      nca=SIZE(board,3)
      call mpp_get_compute_domain (domain_sgs,is,ie,js,je)
      call register_axis(CA_restart, 'xaxis_1', 'X')
      call register_field(CA_restart, 'xaxis_1', 'double', (/'xaxis_1'/))
      call register_variable_attribute(CA_restart, 'xaxis_1', 'cartesian_axis', 'X', str_len=1)
      call get_global_io_domain_indices(CA_restart, 'xaxis_1', is, ie, indices=buffer)
      call write_data(CA_restart, "xaxis_1", buffer)
      deallocate(buffer)

      call register_axis(CA_restart, 'yaxis_1', 'Y')
      call register_field(CA_restart, 'yaxis_1', 'double', (/'yaxis_1'/))
      call register_variable_attribute(CA_restart, 'yaxis_1', 'cartesian_axis', 'Y', str_len=1)
      call get_global_io_domain_indices(CA_restart, 'yaxis_1', js, je, indices=buffer)
      call write_data(CA_restart, "yaxis_1", buffer)
      deallocate(buffer)

      call register_axis(CA_restart, 'zaxis_1', nca )
      call register_field(CA_restart, 'zaxis_1', 'double', (/'zaxis_1'/))
      call register_variable_attribute(CA_restart, 'zaxis_1', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(nca) )
      do i=1, nca
         buffer(i)=i
      end do
      call write_data(CA_restart, "zaxis_1", buffer)
      deallocate(buffer)
      call register_restart_field(CA_restart, "board", board(:,:,:), dimensions=(/'xaxis_1','yaxis_1','zaxis_1'/),is_optional=.false.)
      call register_restart_field(CA_restart, "lives", lives(:,:,:), dimensions=(/'xaxis_1','yaxis_1','zaxis_1'/),is_optional=.false.)
      call write_restart(CA_restart)
      call close_file(CA_restart)
   else
      call mpp_error(FATAL, 'Error opening file '//trim(infile))
   endif
endif
if (allocated(board_g)) then
   if ( amiopen) then 
      amiopen=open_file(CA_restart, trim(infile), 'append', domain=domain_global, is_restart=.true., dont_add_res_to_filename=.true.)
   else
      amiopen=open_file(CA_restart, trim(infile), 'overwrite', domain=domain_global, is_restart=.true., dont_add_res_to_filename=.true.)
   endif
   if( amiopen ) then
      nca_g=SIZE(board_g,3)
      call mpp_get_compute_domain (domain_global,is,ie,js,je)
      call register_axis(CA_restart, 'xaxis_2', 'X')
      call register_field(CA_restart, 'xaxis_2', 'double', (/'xaxis_2'/))
      call register_variable_attribute(CA_restart, 'xaxis_2', 'cartesian_axis', 'X', str_len=1)
      call get_global_io_domain_indices(CA_restart, 'xaxis_2', is, ie, indices=buffer)
      call write_data(CA_restart, "xaxis_2", buffer)
      deallocate(buffer)

      call register_axis(CA_restart, 'yaxis_2', 'Y')
      call register_field(CA_restart, 'yaxis_2', 'double', (/'yaxis_2'/))
      call register_variable_attribute(CA_restart, 'yaxis_2', 'cartesian_axis', 'Y', str_len=1)
      call get_global_io_domain_indices(CA_restart, 'yaxis_2', js, je, indices=buffer)
      call write_data(CA_restart, "yaxis_2", buffer)
      deallocate(buffer)

      call register_axis(CA_restart, 'zaxis_2', nca_g)
      call register_field(CA_restart, 'zaxis_2', 'double', (/'zaxis_2'/))
      call register_variable_attribute(CA_restart, 'zaxis_2', 'cartesian_axis', 'Z', str_len=1)
      allocate( buffer(nca_g) )
      do i=1, nca_g
         buffer(i)=i
      end do
      call write_data(CA_restart, "zaxis_2", buffer)
      deallocate(buffer)
      call register_restart_field(CA_restart, "board_g", board_g(:,:,:), dimensions=(/'xaxis_2','yaxis_2','zaxis_2'/),is_optional=.false.)
      call register_restart_field(CA_restart, "lives_g", lives_g(:,:,:), dimensions=(/'xaxis_2','yaxis_2','zaxis_2'/),is_optional=.false.)
      call write_restart(CA_restart)
      call close_file(CA_restart)
   else
      call mpp_error(FATAL, 'Error opening file '//trim(infile))
   endif
endif


end subroutine write_ca_restart

subroutine read_ca_restart(domain_in,ncells,nca,ncells_g,nca_g)
!Read restart files
implicit none
type(FmsNetcdfDomainFile_t) :: CA_restart
type(domain2D), intent(inout) :: domain_in
integer,intent(in) :: ncells,nca,nca_g,ncells_g
character(len=32)  :: fn_ca = 'ca_data.nc'

character(len=64) :: fname
integer :: id_restart
integer :: nxc,nyc,i
real    :: pi,re,dx
integer :: nx,ny
character(5)  :: indir='INPUT'
logical :: amiopen
integer, allocatable, dimension(:) :: io_layout(:)


call mpp_get_global_domain(domain_in,xsize=nx,ysize=ny,position=CENTER)

 fname = trim(indir)//'/'//trim(fn_ca)
 if (nca .gt. 0 ) then
    allocate(io_layout(2))
    io_layout=mpp_get_io_domain_layout(domain_in)
    call define_ca_domain(domain_in,domain_sgs,ncells,nxncells,nyncells)
    call mpp_define_io_domain(domain_sgs, io_layout)
    call mpp_get_compute_domain (domain_sgs,iscnx,iecnx,jscnx,jecnx)
    amiopen=open_file(CA_restart, trim(fname), 'read', domain=domain_sgs, is_restart=.true., dont_add_res_to_filename=.true.)
    if( amiopen ) then
       call register_axis(CA_restart, 'xaxis_1', 'X')
       call register_axis(CA_restart, 'yaxis_1', 'Y')
       call register_axis(CA_restart, 'nca', nca)
       !Get CA SGS domain

       nxc = iecnx-iscnx+1
       nyc = jecnx-jscnx+1
       if (.not. allocated(board))then
          allocate(board(nxc,nyc,nca))
       endif
       if (.not. allocated(lives))then
          allocate(lives(nxc,nyc,nca))
       endif
   
      !Read restart
      call register_restart_field(CA_restart, "board", board(:,:,:), dimensions=(/'xaxis_1','yaxis_1','zaxis_1'/),is_optional=.false.)
      call register_restart_field(CA_restart, "lives", lives(:,:,:), dimensions=(/'xaxis_1','yaxis_1','zaxis_1'/),is_optional=.false.)
      !--- read the CA restart data
      call mpp_error(NOTE,'reading CA_sgs restart data from INPUT/ca_data.tile*.nc')
      call read_restart(CA_restart)
      call close_file(CA_restart)
      cold_start_ca_sgs=.false.
    else
      call mpp_error(NOTE,'No CA_sgs restarts - cold starting CA')
      cold_start_ca_sgs=.true.
    endif
endif

if (nca_g .gt. 0 ) then
   domain_global=domain_in
   amiopen=open_file(CA_restart, trim(fname), 'read', domain=domain_global, is_restart=.true., dont_add_res_to_filename=.true.)
   if( amiopen ) then
      call register_axis(CA_restart, 'xaxis_2', 'X')
      call register_axis(CA_restart, 'yaxis_2', 'Y')
      call register_axis(CA_restart, 'nca_g', nca_g)
      !call define_ca_domain(domain_in,domain_global,ncells_g,nxncells_g,nyncells_g)
      call mpp_get_compute_domain (domain_global,iscnx_g,iecnx_g,jscnx_g,jecnx_g)
      nxc = iecnx_g-iscnx_g+1
      nyc = jecnx_g-jscnx_g+1
      if (.not. allocated(board_g))then
         allocate(board_g(nxc,nyc,nca_g))
      endif
      if (.not. allocated(lives_g))then
         allocate(lives_g(nxc,nyc,nca_g))
      endif

      !Read restart
      call register_restart_field(CA_restart, "board_g", board_g(:,:,:), dimensions=(/'xaxis_2','yaxis_2','zaxis_2'/),is_optional=.false.)
      call register_restart_field(CA_restart, "lives_g", lives_g(:,:,:), dimensions=(/'xaxis_2','yaxis_2','zaxis_2'/),is_optional=.false.)
      call mpp_error(NOTE,'reading CA_global restart data from INPUT/ca_data.tile*.nc')
      call read_restart(CA_restart)
      call close_file(CA_restart)
      cold_start_ca_global=.false.
      
   else
      call mpp_error(NOTE,'No CA_global restarts - cold starting CA')
      cold_start_ca_global=.true.
   endif
endif


end subroutine read_ca_restart

subroutine update_cells_sgs(kstep,initialize_ca,iseed_ca,first_flag,restart,first_time_step,nca,nxc,nyc,nxch,nych,nlon,&
                            nlat,isc,iec,jsc,jec, npx,npy,  &
                            CA,ca_plumes,iini,ilives_in,nlives,     &
                            nfracseed,nseed,nspinup,nf,nca_plumes,ncells,mytile)

implicit none

integer, intent(in)  :: kstep,nxc,nyc,nlon,nlat,nxch,nych,nca,isc,iec,jsc,jec,npx,npy
integer(kind=kind_dbl_prec), intent(in) :: iseed_ca
integer, intent(in)  :: iini(nxc,nyc,nca),initialize_ca,ilives_in(nxc,nyc,nca)
integer, intent(in)  :: mytile
real,    intent(out) :: CA(nlon,nlat)
integer, intent(out) :: ca_plumes(nlon,nlat)
integer, intent(in)  :: nlives,nseed, nspinup, nf,ncells
real,    intent(in)  :: nfracseed
logical, intent(in)  :: nca_plumes,restart,first_flag,first_time_step
integer, allocatable  :: V(:),L(:),B(:)
integer, allocatable  :: AG(:,:)
integer              :: inci, incj, i, j, k,sub,spinup,it,halo,k_in,isize,jsize
integer              :: ih, jh,kend, boardmax,livemax
real,    allocatable :: board_halo(:,:,:)
integer, dimension(nxc,nyc) :: neighbours, birth, thresh
integer, dimension(nxc,nyc) :: newcell, temp,newseed
integer, dimension(ncells,ncells) :: onegrid
integer(8)           :: nx_full,ny_full
integer(8)           :: iscale = 10000000000
logical, save        :: start_from_restart

real, dimension(nxc,nyc) :: noise_b
integer(8) :: count, count_rate, count_max, count_trunc
integer    :: count4
integer*8            :: i1,j1
real                 :: ncells2inv


!------------------------------------------------------------------------------------------------

if(first_time_step)then
start_from_restart = .False.
endif

!-------------------------------------------------------------------------------------------------
halo=1
isize=nlon+2*halo
jsize=nlat+2*halo
k_in=1
 
  if (.not. allocated(board))then
     allocate(board(nxc,nyc,nca))
     board=0.0
  endif
  if (.not. allocated(lives))then
     allocate(lives(nxc,nyc,nca))
     lives=0.0
  endif
  if(.not. allocated(board_halo))then
     allocate(board_halo(nxch,nych,1))
  endif
 

 !Step 2: Initialize CA, if restart data exist (board,lives > 0) initialize from restart file, otherwise initialize at time-
 !step initialize_ca.
 boardmax=maxval(board)
 call mp_reduce_max(boardmax)
 livemax=maxval(lives)
 call mp_reduce_max(livemax)

 if(restart .and. first_time_step .and. boardmax > 0 .and. livemax > 0)then
    !restart
    start_from_restart = .true.
    spinup = 1
 else

   if(kstep < initialize_ca .and. .not. start_from_restart)then
    do j=1,nyc
     do i=1,nxc
      board(i,j,nf) = 0
      lives(i,j,nf) = 0
     enddo
    enddo
   endif

  if(kstep == initialize_ca .and. .not. start_from_restart)then 
   do j=1,nyc
    do i=1,nxc
    board(i,j,nf) = iini(i,j,nf)
    lives(i,j,nf) = ilives_in(i,j,nf)*iini(i,j,nf)
   enddo
   enddo
   spinup=nspinup
  else
   spinup=1
  endif

 endif

  newseed = 0

!seed with new active cells each nseed time-step regardless of restart/cold start

nx_full=int(ncells,kind=8)*int(npx-1,kind=8)
ny_full=int(ncells,kind=8)*int(npy-1,kind=8)
if(mod(kstep,nseed)==0. .and. (kstep >= initialize_ca .or. start_from_restart))then
   do j=1,nyc
      j1=j+(jsc-1)*ncells
      do i=1,nxc
         i1=i+(isc-1)*ncells
         if (iseed_ca <= 0) then
            !call system_clock(count, count_rate, count_max)
            count_trunc = iscale*(count/iscale)
            count4 = count - count_trunc + mytile *( i1+nx_full*(j1-1)) ! no need to multply by 7 since time will be different in sgs
         else
            count4 = mod((iseed_ca*nf+mytile)*(i1+nx_full*(j1-1))+ 2147483648, 4294967296) - 2147483648
         endif
         noise_b(i,j)=real(random_01_CB(kstep,count4),kind=8)
      enddo
   enddo
   do j=1,nyc
      do i=1,nxc
         if(board(i,j,nf) == 0 .and. noise_b(i,j)>0.90 )then
            newseed(i,j) = 1
         endif
         board(i,j,nf) = board(i,j,nf) + newseed(i,j)
      enddo
   enddo
endif

 
 !Step 3: Evolve CA
 do it = 1,spinup
 
 CA=0
 neighbours=0
 birth=0
 newcell=0

 
 !--- copy board into the halo-augmented board_halo                                                         
 board_halo(1+halo:nxc+halo,1+halo:nyc+halo,1) = real(board(1:nxc,1:nyc,1),kind=8)
! write(1000+mpp_pe(),*) "board_halo pre: ",board_halo(20,1:50,1)

 !--- perform halo update
 call atmosphere_scalar_field_halo (board_halo, halo, nxch, nych, 1, &
                                     iscnx, iecnx, jscnx, jecnx, &
                                     nxncells, nyncells, domain_sgs)

 !--- output data to ensure proper update                                                                            
 !write(1000+mpp_pe(),*) "board_halo post: ",board_halo(20,1:50,1)

 !--- Count the neighbours
  do j=1,nyc
     do i=1,nxc
        ih=i+halo
        jh=j+halo
        neighbours(i,j)=board_halo(ih-1,jh-1,1)+board_halo(ih-1,jh,1)+ &
                        board_halo(ih-1,jh+1,1)+board_halo(ih,jh+1,1)+board_halo(ih+1,jh+1,1)+&
                        board_halo(ih+1,jh,1)+board_halo(ih+1,jh-1,1)+board_halo(ih,jh-1,1)
     enddo
  enddo

 !--- Check rules; 

 !birth
  do j=1,nyc
   do i=1,nxc
     if((neighbours(i,j) == 3 .or. neighbours(i,j) == 2))then 
     birth(i,j)=1
     endif
   enddo
  enddo
 
 !death                                                                                                                        
  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j) < 2 .or. neighbours(i,j) > 3)then  
     lives(i,j,nf)=lives(i,j,nf) - 1
     endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
   if(lives(i,j,nf) < 0)then
     lives(i,j,nf)=0
   endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(birth(i,j)==1 .and. lives(i,j,nf)==0)then
    newcell(i,j)=1
    endif
   enddo
  enddo


  do j=1,nyc
     do i=1,nxc
        lives(i,j,nf)=lives(i,j,nf)+newcell(i,j)*ilives_in(i,j,nf)
     enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(neighbours(i,j)==3 .or. (board(i,j,nf)==1 .and. neighbours(i,j)==2))then
    board(i,j,nf)=1
    else
    board(i,j,nf)=0
    endif
   enddo
  enddo


 enddo !spinup


!COARSE-GRAIN BACK TO NWP GRID
 
  inci=ncells
  incj=ncells
  sub=ncells-1
  ncells2inv=real(1.0/(ncells*ncells))
  DO j=1,nlat
     DO i=1,nlon
        CA(i,j)=(SUM(lives(inci-sub:inci,incj-sub:incj,nf)))*ncells2inv
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

if(nca_plumes) then
!COMPUTE NUMBER OF CLUSTERS (CONVECTIVE PLUMES) IN EACH CA-CELL
!Note, at the moment we only use the count of the plumes found in a grid-cell
!In the future the routine "plumes" can also be used to give the size of 
!each individual plume for better coupling to the convection scheme.

  temp=0
  do j=1,nyc
   do i=1,nxc
     if(lives(i,j,1) > 0)then
      temp(i,j)=1  
     endif
   enddo
  enddo

  kend=ceiling((ncells*ncells)/2.)
  if (.not. allocated(V))then
  allocate(V(kend))
  endif
  if (.not. allocated(L))then
  allocate(L(kend))
  endif
  if (.not. allocated(B))then
  allocate(B(kend))
  endif
  if (.not. allocated(AG))then
  allocate(AG(ncells,ncells))
  endif
  
  ca_plumes(:,:)=0
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        B(:)=0
        L(:)=0
        V(:)=0
        onegrid(1:ncells,1:ncells)=temp(inci-sub:inci,incj-sub:incj)
        call plumes(V,L,AG,onegrid,ncells,ncells,kend)
        do k=1,kend
           if(V(k)==1)then
              B(k)=L(k) !to avoid considering clusters of 0
           endif
        enddo
        ca_plumes(i,j)=MAXVAL(B(1:kend))
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

else

ca_plumes(:,:)=0.

endif ! nca_plumes

end subroutine update_cells_sgs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine update_cells_global(kstep,first_time_step,iseed_ca,restart,nca,nxc,nyc,nxch,nych,nlon,nlat,isc,iec,jsc,jec, &
                        npx,npy,CA,iini_g,ilives_g,                &
                        nlives,ncells,nfracseed,nseed,nspinup,nf,mytile)

implicit none

integer, intent(in) :: kstep,nxc,nyc,nlon,nlat,nxch,nych,nca,isc,iec,jsc,jec,npx,npy
integer, intent(in) :: iini_g(nxc,nyc,nca), ilives_g(nxc,nyc)
integer(kind=kind_dbl_prec), intent(in) :: iseed_ca
real, intent(out) :: CA(nlon,nlat)
logical, intent(in) :: first_time_step
logical, intent(in) :: restart
integer, intent(in) :: nlives, ncells, nseed, nspinup, nf
real, intent(in) :: nfracseed
integer, intent(in) :: mytile
integer,allocatable :: V(:),L(:)
integer :: inci, incj, i, j, k ,sub,spinup,it,halo,k_in,isize,jsize
integer :: ih, jh,kend
real, allocatable :: board_halo(:,:,:)
integer, dimension(nxc,nyc) :: neighbours, birth, thresh
integer, dimension(nxc,nyc) :: newcell, temp,newseed
real, dimension(nxc,nyc) :: noise_b
integer(8) :: count, count_rate, count_max, count_trunc
integer    :: count4
integer(8) :: nx_full,ny_full
integer(8) :: iscale = 10000000000
integer*8            :: i1,j1

!-------------------------------------------------------------------------------------------------

halo=1
isize=nlon+2*halo
jsize=nlat+2*halo
k_in=1

 if (.not. allocated(board_g)) allocate(board_g(nxc,nyc,nca))
 if (.not. allocated(lives_g)) allocate(lives_g(nxc,nyc,nca))
 if (.not. allocated(board_halo)) allocate(board_halo(nxch,nych,1))   

  if(first_time_step .and. cold_start_ca_global)then
   do j=1,nyc
    do i=1,nxc
     board_g(i,j,nf) = iini_g(i,j,nf)
     lives_g(i,j,nf) = ilives_g(i,j)*iini_g(i,j,nf)
    enddo
   enddo

  endif

!Seed with new CA cells at each nseed step
newseed=0
if(mod(kstep,nseed) == 0)then
   nx_full=int(npx-1,kind=8)
   ny_full=int(npy-1,kind=8)
   !random numbers:
   do j=1,nyc
      j1=j+(jsc-1)*ncells
      do i=1,nxc
         i1=i+(isc-1)*ncells
         if (iseed_ca <= 0) then
            !call system_clock(count, count_rate, count_max)
            count_trunc = iscale*(count/iscale)
            count4 = count - count_trunc + mytile *( i1+nx_full*(j1-1)) ! no need to multply by 7 since time will be different in sgs
         else
            count4 = mod(iseed_ca*nf+(7*mytile)*(i1+nx_full*(j1-1))+ 2147483648, 4294967296) - 2147483648
         endif
         noise_b(i,j)=real(random_01_CB(kstep,count4),kind=8)
      enddo
   enddo

   do j=1,nyc
      do i=1,nxc
         if(board_g(i,j,nf) == 0 .and. noise_b(i,j)>0.75 )then
            newseed(i,j)=1
         endif
         board_g(i,j,nf) = board_g(i,j,nf) + newseed(i,j)
      enddo
   enddo
endif

  if(first_time_step .and. cold_start_ca_global)then
  spinup=nspinup
  else
  spinup = 1
  endif
 

do it=1,spinup
!Step 2 - Initialize variables to 0 and extract the halo
 
 neighbours=0
 birth=0
 newcell=0
 CA=0
 board_halo=0

!The input to scalar_field_halo needs to be 1D.
!take the updated board_g fields and extract the halo
! in order to have updated values in the halo region. 

 !--- copy board into the halo-augmented board_halo                                                                                                    
 board_halo(1+halo:nxc+halo,1+halo:nyc+halo,1) = real(board_g(1:nxc,1:nyc,nf),kind=8)
 !write(1000+mpp_pe(),*) "board_halo pre: ",board_halo(:,:,1)                                                                                          

 !--- perform halo update                                                                                                                              
 call atmosphere_scalar_field_halo (board_halo, halo, nxch, nych, 1, &
                                     iscnx_g, iecnx_g, jscnx_g, jecnx_g, &
                                     nxncells_g, nyncells_g, domain_global)

  do j=1,nyc
     do i=1,nxc
        ih=i+halo
        jh=j+halo
        neighbours(i,j)=board_halo(ih-1,jh-1,1)+board_halo(ih-1,jh,1)+ &
                        board_halo(ih-1,jh+1,1)+board_halo(ih,jh+1,1)+board_halo(ih+1,jh+1,1)+&
                        board_halo(ih+1,jh,1)+board_halo(ih+1,jh-1,1)+board_halo(ih,jh-1,1)
     enddo
  enddo



  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j)==2 .or. neighbours(i,j)==3)then
     birth(i,j)=1
     endif
   enddo
  enddo


  do j=1,nyc
   do i=1,nxc
     if(neighbours(i,j)<2 .or. neighbours(i,j)>3)then
     lives_g(i,j,nf)=lives_g(i,j,nf) - 1
     endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
   if(lives_g(i,j,nf)<0)then
     lives_g(i,j,nf)=0
   endif
   enddo
  enddo

  do j=1,nyc
   do i=1,nxc
    if(birth(i,j)==1 .and. lives_g(i,j,nf)==0)then
    newcell(i,j)=1
    endif
   enddo
  enddo


  do j=1,nyc
   do i=1,nxc
    lives_g(i,j,nf)=lives_g(i,j,nf)+newcell(i,j)*ilives_g(i,j)
   enddo
  enddo


   do j=1,nyc
   do i=1,nxc
    if( (board_g(i,j,nf) ==1 .and. (neighbours(i,j)==3 .or. neighbours(i,j)==2) ).or. (board_g(i,j,nf)==0 .and. neighbours(i,j)==3) )then
    board_g(i,j,nf)=1
    else
    board_g(i,j,nf)=0
    endif
   enddo
  enddo

enddo !spinup

!COARSE-GRAIN BACK TO NWP GRID
 
  inci=ncells
  incj=ncells
  sub=ncells-1
  DO j=1,nlat
     DO i=1,nlon
        CA(i,j)=(SUM(lives_g(inci-sub:inci,incj-sub:incj,nf)))/real(ncells*ncells)
        inci=inci+ncells
     ENDDO
     inci=ncells
     incj=incj+ncells
  ENDDO

end subroutine update_cells_global

!================================
 ! This subroutine is copied from FMS/test_fms/test_mpp_domains.F90
  ! and modified to make it simpler to use.
  ! domain_decomp in fv_mp_mod.F90 does something similar, but it does a
  ! few other unnecessary things (and requires more arguments).
  subroutine define_cubic_mosaic(domain, ni, nj, layout, pe_start, pe_end, halo)
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: ni, nj
    integer,        intent(in)    :: layout(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    integer,        intent(in)    :: halo
    !--- local variables
    integer                       :: global_indices(4,6), layout2D(2,6)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact
    integer                       :: i

    ntiles = 6
    num_contact = 12
    if(size(pe_start(:)) .NE. 6 .OR. size(pe_end(:)) .NE. 6 ) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of pe_start and pe_end should be 6")
    if(size(layout) .NE. 2) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of layout should be 2")
    do i = 1, 6
       layout2D(:,i) = layout(:)
       global_indices(1,i) = 1
       global_indices(2,i) = ni
       global_indices(3,i) = 1
       global_indices(4,i) = nj
    enddo
!--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1;     tile2(1) = 2
    istart1(1) = ni;  iend1(1) = ni;  jstart1(1) = 1;      jend1(1) = nj
    istart2(1) = 1;   iend2(1) = 1;   jstart2(1) = 1;      jend2(1) = nj

    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni;  jstart1(2) = nj;  jend1(2) = nj
    istart2(2) = 1;      iend2(2) = 1;   jstart2(2) = nj;  jend2(2) = 1

    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1;     tile2(3) = 5
    istart1(3) = 1;   iend1(3) = 1;      jstart1(3) = 1;   jend1(3) = nj
    istart2(3) = ni;  iend2(3) = 1;      jstart2(3) = nj;  jend2(3) = nj

    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni;  jstart1(4) = 1;   jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni;  jstart2(4) = nj;  jend2(4) = nj

    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2;        tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni;  jstart1(5) = nj;  jend1(5) = nj
    istart2(5) = 1;      iend2(5) = ni;  jstart2(5) = 1;   jend2(5) = 1

    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni;  iend1(6) = ni;  jstart1(6) = 1;      jend1(6) = nj
    istart2(6) = ni;  iend2(6) = 1;   jstart2(6) = 1;      jend2(6) = 1

    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;   iend1(7) = ni;  jstart1(7) = 1;   jend1(7) = 1
    istart2(7) = ni;  iend2(7) = ni;  jstart2(7) = nj;  jend2(7) = 1

    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni;  iend1(8) = ni;  jstart1(8) = 1;      jend1(8) = nj
    istart2(8) = 1;   iend2(8) = 1;   jstart2(8) = 1;      jend2(8) = nj

    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni;  jstart1(9) = nj;  jend1(9) = nj
    istart2(9) = 1;      iend2(9) = 1;   jstart2(9) = nj;  jend2(9) = 1

    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni; jstart1(10) = nj; jend1(10) = nj
    istart2(10) = 1;     iend2(10) = ni; jstart2(10) = 1;  jend2(10) = 1

    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni; iend1(11) = ni; jstart1(11) = 1;     jend1(11) = nj
    istart2(11) = ni; iend2(11) = 1;  jstart2(11) = 1;     jend2(11) = 1

    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni; iend1(12) = ni; jstart1(12) = 1;     jend1(12) = nj
    istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;     jend2(12) = nj

    call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, &
         num_contact, tile1, tile2, istart1, iend1, jstart1, jend1, &
         istart2, iend2, jstart2, jend2, pe_start, pe_end, symmetry=.true., &
         whalo=halo, ehalo=halo, shalo=halo, nhalo=halo, &
         name='CA cubic mosaic')

  end subroutine define_cubic_mosaic

end module update_ca
