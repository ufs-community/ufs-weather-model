
!***********************************************************************
!*                   GNU Lesser General Public License                 
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it 
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or 
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be 
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty 
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.  
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'FV3_control' is for initialization and termination
!! of the model, and controls namelist parameters in FV3.
!----------------
! FV control panel
!----------------

module fv_control_stub_mod
   use constants_mod,       only: pi=>pi_8, kappa, radius, grav, rdgas
   use fms_mod,             only: write_version_number, open_namelist_file, &
                                  check_nml_error, close_file, file_exist,  &
                                  get_mosaic_tile_grid
   use fms_io_mod,          only: set_domain
   use fms_io_mod,          only: field_exist, read_data, &
                               get_global_att_value, get_var_att_value
   use mpp_mod,             only: FATAL, mpp_error, mpp_pe, stdlog, &
                                  mpp_npes, mpp_get_current_pelist, &
                                  input_nml_file, get_unit, WARNING, &
                                  read_ascii_file
   use mpp_domains_mod,     only: mpp_get_data_domain, mpp_get_compute_domain, mpp_get_tile_id
   use tracer_manager_mod,  only: tm_get_number_tracers => get_number_tracers, &
                                  tm_get_tracer_index   => get_tracer_index,   &
                                  tm_get_tracer_indices => get_tracer_indices, &
                                  tm_set_tracer_profile => set_tracer_profile, &
                                  tm_get_tracer_names   => get_tracer_names,   &
                                  tm_check_if_prognostic=> check_if_prognostic,&
                                  tm_register_tracers   => register_tracers

   use fv_mp_stub_mod,      only: mp_start, domain_decomp, mp_assign_gid
   use fv_mp_stub_mod,      only: broadcast_domains, setup_master, grids_master_procs
   use fv_mp_stub_mod,      only: MAX_NNEST, MAX_NTILE,fill_corners,XDir,YDir,ng
   use mpp_domains_mod,     only: domain2D
   use mpp_domains_mod,     only: mpp_get_global_domain
   use mpp_domains_mod,     only: mpp_get_C2F_index, mpp_get_F2C_index
   use mpp_domains_mod,     only: CENTER
   use mpp_domains_mod,     only: mpp_update_domains
   use mpp_mod,             only: mpp_send, mpp_sync, mpp_transmit, mpp_set_current_pelist, &
                                  mpp_declare_pelist, mpp_root_pe, mpp_recv, mpp_sync_self, read_input_nml, &
                                  mpp_max
   use mosaic_mod,       only : get_mosaic_ntiles
   use fv_arrays_stub_mod,       only: fv_atmos_type, allocate_fv_atmos_type,R_GRID

   implicit none
   private

#ifdef OVERLOAD_R4
   real    :: too_big  = 1.E8
#else
   real    :: too_big  = 1.E35
#endif
   public :: fv_control_init

   integer, public :: ngrids = 1
   integer :: commID, global_commID

   integer :: halo_update_type = 1 ! 1 for two-interfaces non-block
                                   ! 2 for block
                                   ! 3 for four-interfaces non-block
#ifdef NO_QUAD_PRECISION
! 64-bit precision (kind=8)
 integer, parameter:: f_p = selected_real_kind(15)
#else
! Higher precision (kind=16) for grid geometrical factors:
 integer, parameter:: f_p = selected_real_kind(20)
#endif

! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

 contains

!-------------------------------------------------------------------------------
         
   subroutine fv_control_init(Atm, dt_atmos, this_grid, grids_on_this_pe)

     type(fv_atmos_type), allocatable, intent(inout), target :: Atm(:)
     real,                intent(in)    :: dt_atmos
     integer,             intent(OUT)   :: this_grid
     logical, allocatable, intent(OUT) :: grids_on_this_pe(:)

     character(100) :: pe_list_name, errstring
     integer :: n, npes, pecounter, i, num_family
     integer, allocatable :: global_pelist(:)
     integer, dimension(MAX_NNEST) :: grid_pes = 0
     integer, dimension(MAX_NNEST) :: all_npx = 0
     integer, dimension(MAX_NNEST) :: all_npy = 0
     integer, dimension(MAX_NNEST) :: all_ntiles = 0

     real :: sdt
     integer :: unit, ens_root_pe, tile_id(1)
     integer :: ngrids

     !!!!!!!!!! POINTERS FOR READING NAMELISTS !!!!!!!!!!

     !------------------------------------------
     ! Model Domain parameters
     ! See fv_arrays.F90 for descriptions
     !------------------------------------------
     !CLEANUP module pointers
     character(len=80) , pointer :: grid_name
     character(len=120), pointer :: grid_file
     integer, pointer :: grid_type

     integer , pointer :: npx           
     integer , pointer :: npy           

     integer , pointer :: ntiles        
     integer , pointer :: ndims        
     real(kind=R_GRID), pointer :: deglon_start, deglon_stop, &  ! boundaries of latlon patch
          deglat_start, deglat_stop
     real(kind=R_GRID), pointer :: deglat

     integer, pointer :: layout(:), io_layout(:)
     logical :: nested

     !!!!!!!!!! END POINTERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!

     this_grid = -1 ! default
     call mp_assign_gid
     ens_root_pe = mpp_root_pe()

     ! 2. Set up Atm and PElists

     ngrids = 1
     allocate(Atm(ngrids))
      Atm(1)%gridstruct%bounded_domain=.false.
     npes = mpp_npes() ! now on global pelist

     allocate(global_pelist(npes))
     call mpp_get_current_pelist(global_pelist, commID=global_commID) ! for commID


     allocate(grids_master_procs(ngrids))
     pecounter = 0
     allocate(grids_on_this_pe(ngrids))
     grids_on_this_pe(:) = .false.

     do n=1,ngrids

        if (ngrids == 1 .or. grid_pes(n) == 0) then
           grid_pes(n) = npes - sum(grid_pes)
           if (grid_pes(n) == 0) then
              if ( n > 1 ) then
                 call mpp_error(FATAL, 'Only one zero entry in grid_pes permitted.') 
              else
                 grid_pes(n) = npes
              endif
           endif
        endif

        allocate(Atm(n)%pelist(grid_pes(n)))
        grids_master_procs(n) = pecounter
        do i=1,grid_pes(n)
           if (pecounter >= npes) then
              if (mpp_pe() == 0) then
                 print*, 'ngrids = ', ngrids, ', grid_pes = ', grid_pes(1:ngrids)
              endif
              call mpp_error(FATAL, 'grid_pes assigns more PEs than are available.')
           endif
           Atm(n)%pelist(i) = pecounter + ens_root_pe
           pecounter = pecounter + 1
           Atm(n)%npes_this_grid = grid_pes(n)
        enddo
        Atm(n)%grid_number = n

        !TODO: we are required to use PE name for reading INTERNAL namelist
        ! and the actual file name for EXTERNAL namelists. Need to clean up this code
        if (n == 1) then
           pe_list_name = ''
        else
           write(pe_list_name,'(A4, I2.2)') 'nest', n
        endif
        call mpp_declare_pelist(Atm(n)%pelist, pe_list_name)
        !If nest need to re-initialize internal NML
        if (n > 1) then
           Atm(n)%nml_filename = 'input_'//trim(pe_list_name)//'.nml'
        else
           Atm(n)%nml_filename = 'input.nml'
        endif
     enddo

     do n=1,ngrids        
        !ONE grid per pe
        if (ANY(mpp_pe() == Atm(n)%pelist)) then
           if (this_grid > 0) then
              print*, mpp_pe(), this_grid, n
              call mpp_error(FATAL, " Grid assigned to multiple pes")
           endif
           call mpp_set_current_pelist(Atm(n)%pelist)
           call setup_master(Atm(n)%pelist)
           this_grid = n
           grids_on_this_pe(n) = .true.
        endif

     enddo

     if (pecounter /= npes) then
        if (mpp_pe() == 0) then
           print*, 'npes = ', npes, ', grid_pes = ', grid_pes(1:ngrids)
           call mpp_error(FATAL, 'grid_pes in fv_nest_Nml does not assign all of the available PEs')
        endif
     endif
 
     ! 3pre. 
     
     ! 3. Read namelists, do option processing and I/O

     call set_namelist_pointers(Atm(this_grid))
     call read_namelist_fv_grid_nml
     call read_namelist_fv_core_nml(Atm(this_grid)) ! do options processing here too?

     call mpp_get_current_pelist(Atm(this_grid)%pelist, commID=commID) ! for commID
     call mp_start(commID,halo_update_type)

     ! 4. Set up domains
     all_ntiles(this_grid) = ntiles
     call mpp_max(all_ntiles, ngrids, global_pelist)

     all_npx(this_grid) = npx
     call mpp_max(all_npx, ngrids, global_pelist)

     all_npy(this_grid) = npy
     call mpp_max(all_npy, ngrids, global_pelist)

     do n=1,ngrids
        if (n/=this_grid) then
           Atm(n)%flagstruct%npx = all_npx(n)
           Atm(n)%flagstruct%npy = all_npy(n)
           Atm(n)%flagstruct%ntiles = all_ntiles(n)
        endif

     enddo

         
     ! 5. domain_decomp()
     nested=.false.
     call domain_decomp(Atm(this_grid)%flagstruct%npx,Atm(this_grid)%flagstruct%npy,Atm(this_grid)%flagstruct%ntiles,&
          Atm(this_grid)%flagstruct%grid_type,nested, &
          Atm(this_grid)%layout,Atm(this_grid)%io_layout,Atm(this_grid)%bd,Atm(this_grid)%tile_of_mosaic, &
          Atm(this_grid)%gridstruct%square_domain,Atm(this_grid)%npes_per_tile,Atm(this_grid)%domain, &
          Atm(this_grid)%domain_for_coupler,Atm(this_grid)%num_contact,Atm(this_grid)%pelist)
     call set_domain(Atm(this_grid)%domain)
     call broadcast_domains(Atm,Atm(this_grid)%pelist,size(Atm(this_grid)%pelist))
     do n=1,ngrids
        tile_id = mpp_get_tile_id(Atm(n)%domain)       
        Atm(n)%global_tile = tile_id(1) ! only meaningful locally
        Atm(n)%npes_per_tile = size(Atm(n)%pelist)/Atm(n)%flagstruct%ntiles ! domain decomp doesn't set this globally
     enddo

     ! 6. Set up domain and Atm structure
     do n=1,ngrids
        call allocate_fv_atmos_type(Atm(n), &
             Atm(n)%bd%isd, Atm(n)%bd%ied, &
             Atm(n)%bd%jsd, Atm(n)%bd%jed, &
             Atm(n)%bd%isc, Atm(n)%bd%iec, &
             Atm(n)%bd%jsc, Atm(n)%bd%jec, &
             Atm(n)%flagstruct%npx,    Atm(n)%flagstruct%npy,   Atm(n)%flagstruct%npz, &
             Atm(n)%flagstruct%ndims,   &
             n/=this_grid, n==this_grid, ngrids) !TODO don't need both of the last arguments
     enddo
     call init_grid(Atm(this_grid), Atm(this_grid)%flagstruct%grid_name, Atm(this_grid)%flagstruct%grid_file, &
          Atm(this_grid)%flagstruct%npx, Atm(this_grid)%flagstruct%npy, Atm(this_grid)%flagstruct%ndims, Atm(this_grid)%flagstruct%ntiles, ng)


   contains
!>@brief The subroutine 'setup_namelist_pointers' associates the MODULE flag pointers
!! with the ARRAY flag variables for the grid active on THIS pe so the flags
!! can be read in from the namelist.
     subroutine set_namelist_pointers(Atm)
       type(fv_atmos_type), intent(INOUT), target :: Atm

       !This routine associates the MODULE flag pointers with the ARRAY flag variables for the grid active on THIS pe so the flags can be read in from the namelist.

       grid_type                     => Atm%flagstruct%grid_type
       grid_name                     => Atm%flagstruct%grid_name
       grid_file                     => Atm%flagstruct%grid_file
       npx                           => Atm%flagstruct%npx
       npy                           => Atm%flagstruct%npy
       ntiles                        => Atm%flagstruct%ntiles
       ndims                         => Atm%flagstruct%ndims
       layout                        => Atm%layout
       io_layout                     => Atm%io_layout

     end subroutine set_namelist_pointers


     subroutine read_namelist_fv_grid_nml

       integer :: f_unit, ios, ierr
       !  local version of these variables to allow PGI compiler to compile
       character(len=80)  :: grid_name = ''
       character(len=120) :: grid_file = ''
       namelist /fv_grid_nml/ grid_name, grid_file

       f_unit=open_namelist_file()
       rewind (f_unit)
       ! Read Main namelist
       read (f_unit,fv_grid_nml,iostat=ios)
       ierr = check_nml_error(ios,'fv_grid_nml')
       call close_file (f_unit)
       call write_version_number ( 'FV_CONTROL_MOD', version )
       unit = stdlog()
       write(unit, nml=fv_grid_nml)

       !Basic option processing
       if (len_trim(grid_file) /= 0) Atm(this_grid)%flagstruct%grid_file = grid_file
       if (len_trim(grid_name) /= 0) Atm(this_grid)%flagstruct%grid_name = grid_name


     end subroutine read_namelist_fv_grid_nml

     subroutine read_namelist_fv_core_nml(Atm)

       type(fv_atmos_type), intent(inout) :: Atm
       integer :: f_unit, ios, ierr
       real :: dim0 = 180.           ! base dimension
       real :: dt0  = 1800.          ! base time step
       real :: dimx, dl, dp, dxmin, dymin, d_fac
       real :: umax = 350.           ! max wave speed for grid_type>3

       integer :: n0split

       !  local version of these variables to allow PGI compiler to compile

       namelist /fv_core_nml/npx, npy, ntiles, layout, io_layout, grid_type


       f_unit = open_namelist_file(Atm%nml_filename)
       ! Read FVCORE namelist
       read (f_unit,fv_core_nml,iostat=ios)
       ierr = check_nml_error(ios,'fv_core_nml')
       call close_file(f_unit)
       call write_version_number ( 'FV_CONTROL_MOD', version )
       unit = stdlog()
       write(unit, nml=fv_core_nml)

       !*** single tile for Cartesian grids
       if (grid_type>3) then
          ntiles=1
       else
          ntiles=6
       endif


197    format(A,l7)
198    format(A,i2.2,A,i4.4,'x',i4.4,'x',i1.1,'-',f9.3)
199    format(A,i3.3)

     end subroutine read_namelist_fv_core_nml

   end subroutine fv_control_init

!>@brief The subroutine 'init_grid' reads the grid from the input file
!! and sets up grid descriptors.
  subroutine init_grid(Atm, grid_name, grid_file, npx, npy, ndims, nregions, ng)
!--------------------------------------------------------
    type(fv_atmos_type), intent(inout), target :: Atm
    character(len=80), intent(IN) :: grid_name
    character(len=120),intent(IN) :: grid_file
    integer,      intent(IN) :: npx, npy
    integer,      intent(IN) :: ndims
    integer,      intent(IN) :: nregions
    integer,      intent(IN) :: ng
!--------------------------------------------------------
    real(kind=R_GRID)   ::  ys(npx,npy)

    real(kind=R_GRID)  :: dp, dl
    real(kind=R_GRID)  :: x1,x2,y1,y2,z1,z2
    integer :: i,j,k,n,nreg
    integer :: fileLun

    real(kind=R_GRID)  :: p1(3), p2(3), p3(3), p4(3)
    real(kind=R_GRID)  :: dist,dist1,dist2, pa(2), pa1(2), pa2(2), pb(2)
    real(kind=R_GRID)  :: pt(3), pt1(3), pt2(3), pt3(3)

    real(kind=R_GRID)  :: angN,angM,angAV,ang
    real(kind=R_GRID)  :: aspN,aspM,aspAV,asp

    real(kind=R_GRID)  :: vec1(3), vec2(3), vec3(3), vec4(3)
    real(kind=R_GRID)  :: vecAvg(3), vec3a(3), vec3b(3), vec4a(3), vec4b(3)
    real(kind=R_GRID)  :: xyz1(3), xyz2(3)

    integer :: ios, ip, jp
    
    integer :: igrid
    
    integer :: tmplun
    character(len=80) :: tmpFile   

    real(kind=R_GRID), dimension(Atm%bd%is:Atm%bd%ie) :: sbuffer, nbuffer
    real(kind=R_GRID), dimension(Atm%bd%js:Atm%bd%je) :: wbuffer, ebuffer

    real(kind=R_GRID), pointer, dimension(:,:,:) :: agrid, grid

    real(kind=R_GRID), pointer, dimension(:,:) :: sina, cosa

    integer, pointer, dimension(:,:,:) ::  iinta, jinta, iintb, jintb

    integer, pointer :: ntiles_g, tile
    logical, pointer :: sw_corner, se_corner, ne_corner, nw_corner
    logical, pointer :: latlon, cubed_sphere, have_south_pole, have_north_pole

    type(domain2d), pointer :: domain
    integer :: is,  ie,  js,  je
    integer :: isd, ied, jsd, jed

    is  = Atm%bd%is
    ie  = Atm%bd%ie
    js  = Atm%bd%js
    je  = Atm%bd%je
    isd = Atm%bd%isd
    ied = Atm%bd%ied
    jsd = Atm%bd%jsd
    jed = Atm%bd%jed

    !!! Associate pointers
    agrid => Atm%gridstruct%agrid_64
    grid  => Atm%gridstruct%grid_64

    iinta                         => Atm%gridstruct%iinta
    jinta                         => Atm%gridstruct%jinta
    iintb                         => Atm%gridstruct%iintb
    jintb                         => Atm%gridstruct%jintb
    ntiles_g                      => Atm%gridstruct%ntiles_g
    sw_corner                     => Atm%gridstruct%sw_corner
    se_corner                     => Atm%gridstruct%se_corner
    ne_corner                     => Atm%gridstruct%ne_corner
    nw_corner                     => Atm%gridstruct%nw_corner
    latlon                        => Atm%gridstruct%latlon
    cubed_sphere                  => Atm%gridstruct%cubed_sphere
    have_south_pole               => Atm%gridstruct%have_south_pole
    have_north_pole               => Atm%gridstruct%have_north_pole

    tile                          => Atm%tile_of_mosaic

    domain                        => Atm%domain

    ntiles_g = nregions
    latlon = .false.
    cubed_sphere = .true.
        
    call read_grid(Atm, grid_file, ndims, nregions, ng)

    call sorted_inta(isd, ied, jsd, jed, cubed_sphere, grid, iinta, jinta)

    agrid(:,:,:) = -1.e25
 
    do j=js,je
       do i=is,ie
            call cell_center2(grid(iinta(1,i,j),jinta(1,i,j),1:2),  &
                              grid(iinta(2,i,j),jinta(2,i,j),1:2),  &
                              grid(iinta(3,i,j),jinta(3,i,j),1:2),  &
                              grid(iinta(4,i,j),jinta(4,i,j),1:2),  &
                              agrid(i,j,1:2) )
       enddo
    enddo
    call sorted_intb(isd, ied, jsd, jed, is, ie, js, je, npx, npy, &
                         cubed_sphere, agrid, iintb, jintb)
    end subroutine init_grid


!-------------------------------------------------------------------------------
 subroutine cell_center2(q1, q2, q3, q4, e2)
      real(kind=R_GRID) , intent(in ) :: q1(2), q2(2), q3(2), q4(2)
      real(kind=R_GRID) , intent(out) :: e2(2)
! Local
      real(kind=R_GRID) p1(3), p2(3), p3(3), p4(3)
      real(kind=R_GRID) ec(3)
      real(kind=R_GRID) dd
      integer k

      call latlon2xyz(q1, p1)
      call latlon2xyz(q2, p2)
      call latlon2xyz(q3, p3)
      call latlon2xyz(q4, p4)

      do k=1,3
         ec(k) = p1(k) + p2(k) + p3(k) + p4(k)
      enddo
      dd = sqrt( ec(1)**2 + ec(2)**2 + ec(3)**2 )

      do k=1,3
         ec(k) = ec(k) / dd
      enddo
      call cart_to_latlon(1, ec, e2(1), e2(2))

 end subroutine cell_center2
!>@brief The subroutine 'read_grid' reads the grid from the mosaic grid file.
  subroutine read_grid(Atm, grid_file, ndims, nregions, ng)
    type(fv_atmos_type), intent(inout), target :: Atm
    character(len=*),    intent(IN)    :: grid_file
    integer,             intent(IN)    :: ndims
    integer,             intent(IN)    :: nregions
    integer,             intent(IN)    :: ng

    real, allocatable, dimension(:,:)  :: tmpx, tmpy
    real(kind=R_GRID), pointer, dimension(:,:,:)    :: grid
    character(len=128)                 :: units = ""
    character(len=256)                 :: atm_mosaic, atm_hgrid, grid_form
    character(len=1024)                :: attvalue
    integer                            :: ntiles, i, j, stdunit
    integer                            :: isc2, iec2, jsc2, jec2
    integer                            :: start(4), nread(4)  
    integer                            :: is,  ie,  js,  je
    integer                            :: isd, ied, jsd, jed
    integer,save :: halo=3 ! for regional domain external tools

    is  = Atm%bd%is
    ie  = Atm%bd%ie
    js  = Atm%bd%js
    je  = Atm%bd%je
    isd = Atm%bd%isd
    ied = Atm%bd%ied
    jsd = Atm%bd%jsd
    jed = Atm%bd%jed
    grid  => Atm%gridstruct%grid_64

    if(.not. file_exist(grid_file)) call mpp_error(FATAL, 'fv_grid_tools(read_grid): file '// &
         trim(grid_file)//' does not exist')

    !--- make sure the grid file is mosaic file.
    if(field_exist(grid_file, 'atm_mosaic_file')) then
       call read_data(grid_file, "atm_mosaic_file", atm_mosaic)
       atm_mosaic = "INPUT/"//trim(atm_mosaic)
    else 
       atm_mosaic = trim(grid_file)
    endif

    call get_mosaic_tile_grid(atm_hgrid, atm_mosaic, Atm%domain)

    grid_form = "none"    
    if( get_global_att_value(atm_hgrid, "history", attvalue) ) then
       if( index(attvalue, "gnomonic_ed") > 0) grid_form = "gnomonic_ed"
    endif
    if(grid_form .NE. "gnomonic_ed") call mpp_error(FATAL, &
         "fv_grid_tools(read_grid): the grid should be 'gnomonic_ed' when reading from grid file, contact developer")

    ntiles = get_mosaic_ntiles(atm_mosaic)
    if( .not. Atm%gridstruct%bounded_domain) then  !<-- The regional setup has only 1 tile so do not shutdown in that case.
       if(ntiles .NE. 6) call mpp_error(FATAL, &
            'fv_grid_tools(read_grid): ntiles should be 6 in mosaic file '//trim(atm_mosaic) )
       if(nregions .NE. 6) call mpp_error(FATAL, &
            'fv_grid_tools(read_grid): nregions should be 6 when reading from mosaic file '//trim(grid_file) )
    endif

    call get_var_att_value(atm_hgrid, 'x', 'units', units)

    !--- get the geographical coordinates of super-grid.
    isc2 = 2*is-1; iec2 = 2*ie+1
    jsc2 = 2*js-1; jec2 = 2*je+1  
    if( Atm%gridstruct%bounded_domain ) then
      isc2 = 2*(isd+halo)-1; iec2 = 2*(ied+1+halo)-1   ! For the regional domain the cell corner locations must be transferred
      jsc2 = 2*(jsd+halo)-1; jec2 = 2*(jed+1+halo)-1   ! from the entire supergrid to the compute grid, including the halo region. 
    endif
    allocate(tmpx(isc2:iec2, jsc2:jec2) )
    allocate(tmpy(isc2:iec2, jsc2:jec2) )
    start = 1; nread = 1
    start(1) = isc2; nread(1) = iec2 - isc2 + 1
    start(2) = jsc2; nread(2) = jec2 - jsc2 + 1
    call read_data(atm_hgrid, 'x', tmpx, start, nread, no_domain=.TRUE.)  !<-- tmpx (lon, deg east) is on the supergrid
    call read_data(atm_hgrid, 'y', tmpy, start, nread, no_domain=.TRUE.)  !<-- tmpy (lat, deg) is on the supergrid

    !--- geographic grid at cell corner
    grid(isd: is-1, jsd:js-1,1:ndims)=0.
    grid(isd: is-1, je+2:jed+1,1:ndims)=0.
    grid(ie+2:ied+1,jsd:js-1,1:ndims)=0.
    grid(ie+2:ied+1,je+2:jed+1,1:ndims)=0.
    if(len_trim(units) < 6) call mpp_error(FATAL, &
          "fv_grid_tools_mod(read_grid): the length of units must be no less than 6")
    if(units(1:6) == 'degree') then
    if( .not. Atm%gridstruct%bounded_domain) then
       do j = js, je+1
          do i = is, ie+1
             grid(i,j,1) = tmpx(2*i-1,2*j-1)*pi/180.
             grid(i,j,2) = tmpy(2*i-1,2*j-1)*pi/180.
          enddo
       enddo
    else
!
!***  In the regional case the halo surrounding the domain was included in the read.
!***  Transfer the compute and halo regions to the compute grid.
!
          do j = jsd, jed+1
          do i = isd, ied+1
             grid(i,j,1) = tmpx(2*i+halo+2,2*j+halo+2)*pi/180.
             grid(i,j,2) = tmpy(2*i+halo+2,2*j+halo+2)*pi/180.
          enddo
          enddo
       endif

    else if(units(1:6) == 'radian') then
       do j = js, je+1
          do i = is, ie+1
             grid(i,j,1) = tmpx(2*i-1,2*j-1)
             grid(i,j,2) = tmpy(2*i-1,2*j-1)
          enddo
       enddo
    else
       print*, 'units is ' , trim(units), len_trim(units), mpp_pe()
       call mpp_error(FATAL, 'fv_grid_tools_mod(read_grid): units must start with degree or radian')
    endif

    deallocate(tmpx, tmpy)
    nullify(grid)
  end subroutine read_grid
 subroutine latlon2xyz(p, e, id)

 real(kind=R_GRID), intent(in) :: p(2)
 real(kind=R_GRID), intent(out):: e(3)
 integer, optional, intent(in):: id !< id=0 do nothing; id=1, right_hand

 integer n
 real (f_p):: q(2)
 real (f_p):: e1, e2, e3

    do n=1,2
       q(n) = p(n)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3

 end subroutine latlon2xyz
 subroutine cart_to_latlon(np, q, xs, ys)
! vector version of cart_to_latlon1
  integer, intent(in):: np
  real(kind=R_GRID), intent(inout):: q(3,np)
  real(kind=R_GRID), intent(inout):: xs(np), ys(np)
! local
  real(kind=R_GRID), parameter:: esl=1.d-10
  real (f_p):: p(3)
  real (f_p):: dist, lat, lon
  integer i,k

  do i=1,np
     do k=1,3
        p(k) = q(k,i)
     enddo
     dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
     do k=1,3
        p(k) = p(k) / dist
     enddo

     if ( (abs(p(1))+abs(p(2)))  < esl ) then
          lon = real(0.,kind=f_p)
     else
          lon = atan2( p(2), p(1) )   ! range [-pi,pi]
     endif

     if ( lon < 0.) lon = real(2.,kind=f_p)*pi + lon
! RIGHT_HAND system:
     lat = asin(p(3))

     xs(i) = lon
     ys(i) = lat
! q Normalized:
     do k=1,3
        q(k,i) = p(k)
     enddo
  enddo

 end  subroutine cart_to_latlon
  
subroutine sorted_inta(isd, ied, jsd, jed, cubed_sphere, bgrid, iinta, jinta)
    integer, intent(in) :: isd, ied, jsd, jed
    real(kind=R_GRID),    intent(in), dimension(isd:ied+1,jsd:jed+1,2) :: bgrid
    logical, intent(in) :: cubed_sphere

    integer, intent(out), dimension(4,isd:ied,jsd:jed) :: iinta, jinta
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real,    dimension(4) :: xsort, ysort
    integer, dimension(4) :: isort, jsort
    integer :: i, j
    !------------------------------------------------------------------!
    ! special treatment for cubed sphere                               !
    !------------------------------------------------------------------!
    if (cubed_sphere) then
       !---------------------------------------------------------------!
       ! get order of indices for line integral around a-grid cell     ! 
       !---------------------------------------------------------------!
       do j=jsd,jed
          do i=isd,ied
             xsort(1)=bgrid(i  ,j  ,1); ysort(1)=bgrid(i  ,j  ,2); isort(1)=i  ; jsort(1)=j
             xsort(2)=bgrid(i  ,j+1,1); ysort(2)=bgrid(i  ,j+1,2); isort(2)=i  ; jsort(2)=j+1
             xsort(3)=bgrid(i+1,j+1,1); ysort(3)=bgrid(i+1,j+1,2); isort(3)=i+1; jsort(3)=j+1
             xsort(4)=bgrid(i+1,j  ,1); ysort(4)=bgrid(i+1,j  ,2); isort(4)=i+1; jsort(4)=j
             call sort_rectangle(iinta(1,i,j), jinta(1,i,j))
          enddo
       enddo
    else
       !---------------------------------------------------------------!
       ! default behavior for other grids                              !
       !---------------------------------------------------------------!
       do j=jsd,jed
          do i=isd,ied
             iinta(i,j,1)=i  ; jinta(i,j,1)=j
             iinta(i,j,2)=i  ; jinta(i,j,2)=j+1
             iinta(i,j,3)=i+1; jinta(i,j,3)=j+1
             iinta(i,j,4)=i+1; jinta(i,j,4)=j  
          enddo
       enddo
    endif

  contains
    !------------------------------------------------------------------!
    subroutine sort_rectangle(iind, jind)
      integer, dimension(4), intent(inout) :: iind, jind
      !----------------------------------------------------------------!
      ! local variables                                                !
      !----------------------------------------------------------------!
      real,    dimension(4) :: xsorted, ysorted
      integer, dimension(4) :: isorted, jsorted
      integer :: l, ll, lll
      !----------------------------------------------------------------!
      ! sort in east west                                              !
      !----------------------------------------------------------------!
      xsorted(:)=10.
      ysorted(:)=10.
      isorted(:)=0
      jsorted(:)=0
             
      do l=1,4
         do ll=1,4
            if (xsort(l)<xsorted(ll)) then
               do lll=3,ll,-1
                  xsorted(lll+1)=xsorted(lll)
                  ysorted(lll+1)=ysorted(lll)
                  isorted(lll+1)=isorted(lll)
                  jsorted(lll+1)=jsorted(lll)
               enddo
               xsorted(ll)=xsort(l)
               ysorted(ll)=ysort(l)
               isorted(ll)=isort(l)
               jsorted(ll)=jsort(l)
               exit
            endif
         enddo
      enddo
      !----------------------------------------------------------------!
      ! sort in north south                                            !
      !----------------------------------------------------------------!
      do l=1,4
         xsort(l)=xsorted(l); ysort(l)=ysorted(l)
         isort(l)=isorted(l); jsort(l)=jsorted(l)
      enddo
      xsorted(:)=10.
      ysorted(:)=10.
      isorted(:)=0
      jsorted(:)=0
      
      do l=1,4
         do ll=1,4
            if (ysort(l)<ysorted(ll)) then
               do lll=3,ll,-1
                  xsorted(lll+1)=xsorted(lll)
                  ysorted(lll+1)=ysorted(lll)
                  isorted(lll+1)=isorted(lll)
                  jsorted(lll+1)=jsorted(lll)
               enddo
               xsorted(ll)=xsort(l)
               ysorted(ll)=ysort(l)
               isorted(ll)=isort(l)
               jsorted(ll)=jsort(l)
               exit
            endif
         enddo
      enddo
      !----------------------------------------------------------------!
      ! use first two grid point for start and orientation             !
      !----------------------------------------------------------------!
      if ( isorted(1)==i .and. jsorted(1)==j ) then
         if ( isorted(2)==i+1 .and. jsorted(2)==j+1 ) then
            isorted(2)=isorted(3); jsorted(2)=jsorted(3)
         endif
         if ( isorted(2)==i   .and. jsorted(2)==j+1 ) then
            iind(1)=i  ; jind(1)=j
            iind(2)=i  ; jind(2)=j+1
            iind(3)=i+1; jind(3)=j+1
            iind(4)=i+1; jind(4)=j  
         elseif ( isorted(2)==i+1 .and. jsorted(2)==j ) then
            iind(1)=i  ; jind(1)=j
            iind(2)=i+1; jind(2)=j
            iind(3)=i+1; jind(3)=j+1
            iind(4)=i  ; jind(4)=j+1
         endif
         
      elseif ( isorted(1)==i .and. jsorted(1)==j+1 ) then
         if ( isorted(2)==i+1 .and. jsorted(2)==j ) then
            isorted(2)=isorted(3); jsorted(2)=jsorted(3)
         endif
         if ( isorted(2)==i+1 .and. jsorted(2)==j+1 ) then
            iind(1)=i  ; jind(1)=j+1
            iind(2)=i+1; jind(2)=j+1
            iind(3)=i+1; jind(3)=j
            iind(4)=i  ; jind(4)=j  
         elseif ( isorted(2)==i   .and. jsorted(2)==j ) then
            iind(1)=i  ; jind(1)=j+1
            iind(2)=i  ; jind(2)=j
            iind(3)=i+1; jind(3)=j
            iind(4)=i+1; jind(4)=j+1
         endif
         
      elseif ( isorted(1)==i+1 .and. jsorted(1)==j+1 ) then
         if ( isorted(2)==i .and. jsorted(2)==j ) then
            isorted(2)=isorted(3); jsorted(2)=jsorted(3)
         endif
         if ( isorted(2)==i+1 .and. jsorted(2)==j ) then
            iind(1)=i+1; jind(1)=j+1
            iind(2)=i+1; jind(2)=j
            iind(3)=i  ; jind(3)=j
            iind(4)=i  ; jind(4)=j+1  
         elseif ( isorted(2)==i   .and. jsorted(2)==j+1 ) then
            iind(1)=i+1; jind(1)=j+1
            iind(2)=i  ; jind(2)=j+1
            iind(3)=i  ; jind(3)=j
            iind(4)=i+1; jind(4)=j
         endif
         
      elseif ( isorted(1)==i+1 .and. jsorted(1)==j ) then
         if ( isorted(2)==i .and. jsorted(2)==j+1 ) then
            isorted(2)=isorted(3); jsorted(2)=jsorted(3)
         endif
         if ( isorted(2)==i   .and. jsorted(2)==j ) then
            iind(1)=i+1; jind(1)=j
            iind(2)=i  ; jind(2)=j
            iind(3)=i  ; jind(3)=j+1
            iind(4)=i+1; jind(4)=j+1
         elseif ( isorted(2)==i+1 .and. jsorted(2)==j+1 ) then
            iind(1)=i+1; jind(1)=j
            iind(2)=i+1; jind(2)=j+1
            iind(3)=i  ; jind(3)=j+1
            iind(4)=i  ; jind(4)=j  
         endif
         
      endif

    end subroutine sort_rectangle
    !------------------------------------------------------------------!
  end subroutine sorted_inta

!>@brief The subroutine 'sorted_intb' sorts cell corner indices in latlon space
!!  based on grid locations in index space.
!>@details If not the grid is notcubed_sphere, it assumes that
!! the orientations in index  and latlon space are identical.
!! i/jinta are indices of b-grid locations needed for line integrals 
!! around an a-grid cell including ghosting.
!! i/jintb are indices of a-grid locations needed for line integrals
!! around a b-grid cell, no ghosting.
  subroutine sorted_intb(isd, ied, jsd, jed, is, ie, js, je, npx, npy, &
                          cubed_sphere, agrid, iintb, jintb)
    integer, intent(in) :: isd, ied, jsd, jed, is, ie, js, je, npx, npy
    real(kind=R_GRID),    intent(in), dimension(isd:ied,jsd:jed,2) :: agrid
    logical, intent(in) :: cubed_sphere

    integer, dimension(4,is:ie+1,js:je+1), intent(out) :: iintb, jintb
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real,    dimension(4) :: xsort, ysort, xsorted, ysorted 
    integer, dimension(4) :: isort, jsort, isorted, jsorted
    integer :: i, j, l, ll, lll
    !------------------------------------------------------------------!
    ! special treatment for cubed sphere                               !
    !------------------------------------------------------------------!
    if (cubed_sphere) then
       !---------------------------------------------------------------!
       ! get order of indices for line integral around b-grid cell     ! 
       !---------------------------------------------------------------!
       do j=js,je+1
          do i=is,ie+1
             xsort(1)=agrid(i  ,j  ,1); ysort(1)=agrid(i  ,j  ,2); isort(1)=i  ; jsort(1)=j
             xsort(2)=agrid(i  ,j-1,1); ysort(2)=agrid(i  ,j-1,2); isort(2)=i  ; jsort(2)=j-1
             xsort(3)=agrid(i-1,j-1,1); ysort(3)=agrid(i-1,j-1,2); isort(3)=i-1; jsort(3)=j-1
             xsort(4)=agrid(i-1,j  ,1); ysort(4)=agrid(i-1,j  ,2); isort(4)=i-1; jsort(4)=j
             call sort_rectangle(iintb(1,i,j), jintb(1,i,j))
          enddo
       enddo
       !---------------------------------------------------------------!
       ! take care of corner points                                    !
       !---------------------------------------------------------------!
       if ( (is==1) .and. (js==1) ) then
          i=1
          j=1
          xsort(1)=agrid(i  ,j  ,1); ysort(1)=agrid(i  ,j  ,2); isort(1)=i  ; jsort(1)=j  
          xsort(2)=agrid(i  ,j-1,1); ysort(2)=agrid(i  ,j-1,2); isort(2)=i  ; jsort(2)=j-1
          xsort(3)=agrid(i-1,j  ,1); ysort(3)=agrid(i-1,j  ,2); isort(3)=i-1; jsort(3)=j
          call sort_triangle()
          iintb(4,i,j)=i-1; jintb(4,i,j)=j-1
       endif

       if ( (ie+1==npx) .and. (js==1) ) then
          i=npx
          j=1
          xsort(1)=agrid(i  ,j  ,1); ysort(1)=agrid(i  ,j  ,2); isort(1)=i  ; jsort(1)=j
          xsort(2)=agrid(i-1,j  ,1); ysort(2)=agrid(i-1,j  ,2); isort(2)=i-1; jsort(2)=j
          xsort(3)=agrid(i-1,j-1,1); ysort(3)=agrid(i-1,j-1,2); isort(3)=i-1; jsort(3)=j-1
          call sort_triangle()
          iintb(4,i,j)=i; jintb(4,i,j)=j-1
       endif

       if ( (ie+1==npx) .and. (je+1==npy) ) then
          i=npx
          j=npy
          xsort(1)=agrid(i-1,j-1,1); ysort(1)=agrid(i-1,j-1,2); isort(1)=i-1; jsort(1)=j-1
          xsort(2)=agrid(i  ,j-1,1); ysort(2)=agrid(i  ,j-1,2); isort(2)=i  ; jsort(2)=j-1
          xsort(3)=agrid(i-1,j  ,1); ysort(3)=agrid(i-1,j  ,2); isort(3)=i-1; jsort(3)=j
          call sort_triangle()
          iintb(4,i,j)=i; jintb(4,i,j)=j
       endif
       
       if ( (is==1) .and. (je+1==npy) ) then
          i=1
          j=npy
          xsort(1)=agrid(i  ,j  ,1); ysort(1)=agrid(i  ,j  ,2); isort(1)=i  ; jsort(1)=j
          xsort(2)=agrid(i-1,j-1,1); ysort(2)=agrid(i-1,j-1,2); isort(2)=i-1; jsort(2)=j-1
          xsort(3)=agrid(i  ,j-1,1); ysort(3)=agrid(i  ,j-1,2); isort(3)=i  ; jsort(3)=j-1
          call sort_triangle()
          iintb(4,i,j)=i-1; jintb(4,i,j)=j
       endif
    else
       !---------------------------------------------------------------!
       ! default behavior for other grids                              !
       !---------------------------------------------------------------!
       do j=js,je+1
          do i=is,ie+1
             iintb(1,i,j)=i  ; jintb(1,i,j)=j
             iintb(2,i,j)=i  ; jintb(2,i,j)=j-1
             iintb(3,i,j)=i-1; jintb(3,i,j)=j-1
             iintb(4,i,j)=i-1; jintb(4,i,j)=j  
          enddo
       enddo
    endif

  contains
    !------------------------------------------------------------------!
    subroutine sort_rectangle(iind, jind)

      integer, dimension(4), intent(inout) :: iind, jind
      !----------------------------------------------------------------!
      ! local variables                                                !
      !----------------------------------------------------------------!
      real,    dimension(4) :: xsorted, ysorted 
      integer, dimension(4) :: isorted, jsorted
      !----------------------------------------------------------------!
      ! sort in east west                                              !
      !----------------------------------------------------------------!
      xsorted(:)=10.
      ysorted(:)=10.
      isorted(:)=0
      jsorted(:)=0
             
      do l=1,4
         do ll=1,4
            if (xsort(l)<xsorted(ll)) then
               do lll=3,ll,-1
                  xsorted(lll+1)=xsorted(lll)
                  ysorted(lll+1)=ysorted(lll)
                  isorted(lll+1)=isorted(lll)
                  jsorted(lll+1)=jsorted(lll)
               enddo
               xsorted(ll)=xsort(l)
               ysorted(ll)=ysort(l)
               isorted(ll)=isort(l)
               jsorted(ll)=jsort(l)
               exit
            endif
         enddo
      enddo
      !----------------------------------------------------------------!
      ! sort in north south                                            !
      !----------------------------------------------------------------!
      do l=1,4
         xsort(l)=xsorted(l); ysort(l)=ysorted(l)
         isort(l)=isorted(l); jsort(l)=jsorted(l)
      enddo
      xsorted(:)=10.
      ysorted(:)=10.
      isorted(:)=0
      jsorted(:)=0
      
      do l=1,4
         do ll=1,4
            if (ysort(l)<ysorted(ll)) then
               do lll=3,ll,-1
                  xsorted(lll+1)=xsorted(lll)
                  ysorted(lll+1)=ysorted(lll)
                  isorted(lll+1)=isorted(lll)
                  jsorted(lll+1)=jsorted(lll)
               enddo
               xsorted(ll)=xsort(l)
               ysorted(ll)=ysort(l)
               isorted(ll)=isort(l)
               jsorted(ll)=jsort(l)
               exit
            endif
         enddo
      enddo
      !----------------------------------------------------------------!
      ! use first two grid point for start and orientation             !
      !----------------------------------------------------------------!
      if ( isorted(1)==i .and. jsorted(1)==j ) then
         if ( isorted(2)==i-1 .and. jsorted(2)==j-1 ) then
            isorted(2)=isorted(3); jsorted(2)=jsorted(3)
         endif
         if ( isorted(2)==i   .and. jsorted(2)==j-1 ) then
            iind(1)=i  ; jind(1)=j
            iind(2)=i  ; jind(2)=j-1
            iind(3)=i-1; jind(3)=j-1
            iind(4)=i-1; jind(4)=j  
         elseif ( isorted(2)==i-1 .and. jsorted(2)==j ) then
            iind(1)=i  ; jind(1)=j
            iind(2)=i-1; jind(2)=j
            iind(3)=i-1; jind(3)=j-1
            iind(4)=i  ; jind(4)=j-1
         endif
         
      elseif ( isorted(1)==i .and. jsorted(1)==j-1 ) then
         if ( isorted(2)==i-1 .and. jsorted(2)==j ) then
            isorted(2)=isorted(3); jsorted(2)=jsorted(3)
         endif
         if ( isorted(2)==i-1 .and. jsorted(2)==j-1 ) then
            iind(1)=i  ; jind(1)=j-1
            iind(2)=i-1; jind(2)=j-1
            iind(3)=i-1; jind(3)=j
            iind(4)=i  ; jind(4)=j  
         elseif ( isorted(2)==i   .and. jsorted(2)==j ) then
            iind(1)=i  ; jind(1)=j-1
            iind(2)=i  ; jind(2)=j
            iind(3)=i-1; jind(3)=j
            iind(4)=i-1; jind(4)=j-1
         endif
         
      elseif ( isorted(1)==i-1 .and. jsorted(1)==j-1 ) then
         if ( isorted(2)==i .and. jsorted(2)==j ) then
            isorted(2)=isorted(3); jsorted(2)=jsorted(3)
         endif
         if ( isorted(2)==i-1 .and. jsorted(2)==j ) then
            iind(1)=i-1; jind(1)=j-1
            iind(2)=i-1; jind(2)=j
            iind(3)=i  ; jind(3)=j
            iind(4)=i  ; jind(4)=j-1  
         elseif ( isorted(2)==i   .and. jsorted(2)==j-1 ) then
            iind(1)=i-1; jind(1)=j-1
            iind(2)=i  ; jind(2)=j-1
            iind(3)=i  ; jind(3)=j
            iind(4)=i-1; jind(4)=j
         endif
         
      elseif ( isorted(1)==i-1 .and. jsorted(1)==j ) then
         if ( isorted(2)==i .and. jsorted(2)==j-1 ) then
            isorted(2)=isorted(3); jsorted(2)=jsorted(3)
         endif
         if ( isorted(2)==i   .and. jsorted(2)==j ) then
            iind(1)=i-1; jind(1)=j
            iind(2)=i  ; jind(2)=j
            iind(3)=i  ; jind(3)=j-1
            iind(4)=i-1; jind(4)=j-1
         elseif ( isorted(2)==i-1 .and. jsorted(2)==j-1 ) then
            iind(1)=i-1; jind(1)=j
            iind(2)=i-1; jind(2)=j-1
            iind(3)=i  ; jind(3)=j-1
            iind(4)=i  ; jind(4)=j  
         endif
         
      endif

    end subroutine sort_rectangle
    !------------------------------------------------------------------!
    subroutine sort_triangle()

      xsorted(1:3)=10.
      ysorted(1:3)=10.
      isorted(1:3)=0
      jsorted(1:3)=0
      !----------------------------------------------------------------!
      ! sort in east west                                              !
      !----------------------------------------------------------------!
      do l=1,3
         do ll=1,3
            if (xsort(l)<xsorted(ll)) then
               do lll=2,ll,-1
                  xsorted(lll+1)=xsorted(lll)
                  ysorted(lll+1)=ysorted(lll)
                  isorted(lll+1)=isorted(lll)
                  jsorted(lll+1)=jsorted(lll)
               enddo
               xsorted(ll)=xsort(l)
               ysorted(ll)=ysort(l)
               isorted(ll)=isort(l)
               jsorted(ll)=jsort(l)
               exit
            endif
         enddo
      enddo
      !----------------------------------------------------------------!
      ! sort in north south                                            !
      !----------------------------------------------------------------!
      do l=1,3
         xsort(l)=xsorted(l); ysort(l)=ysorted(l)
         isort(l)=isorted(l); jsort(l)=jsorted(l)
      enddo
      xsorted(1:3)=10.
      ysorted(1:3)=10.
      isorted(1:3)=0
      jsorted(1:3)=0
      
      do l=1,3
         do ll=1,3
            if (ysort(l)<ysorted(ll)) then
               do lll=2,ll,-1
                  xsorted(lll+1)=xsorted(lll)
                  ysorted(lll+1)=ysorted(lll)
                  isorted(lll+1)=isorted(lll)
                  jsorted(lll+1)=jsorted(lll)
               enddo
               xsorted(ll)=xsort(l)
               ysorted(ll)=ysort(l)
               isorted(ll)=isort(l)
               jsorted(ll)=jsort(l)
               exit
            endif
         enddo
      enddo
      !----------------------------------------------------------------!
      ! use first two grid point for start and orientation             !
      !----------------------------------------------------------------!
      iintb(1,i,j)=isorted(1) ; jintb(1,i,j)=jsorted(1)
      iintb(2,i,j)=isorted(2) ; jintb(2,i,j)=jsorted(2)
      iintb(3,i,j)=isorted(3) ; jintb(3,i,j)=jsorted(3)
   
    end subroutine sort_triangle
    !------------------------------------------------------------------!
  end subroutine sorted_intb

end module fv_control_stub_mod
