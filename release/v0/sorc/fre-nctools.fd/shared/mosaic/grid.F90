module grid_mod

use mpp_mod, only : mpp_root_pe
use constants_mod, only : PI, radius
use fms_mod, only : uppercase, lowercase, field_exist, field_size, read_data, &
     error_mesg, string, FATAL, NOTE
use fms_io_mod, only : get_great_circle_algorithm, get_global_att_value
use mosaic_mod, only : get_mosaic_ntiles, get_mosaic_xgrid_size, get_mosaic_grid_sizes, &
     get_mosaic_xgrid, calc_mosaic_grid_area, calc_mosaic_grid_great_circle_area

! the following two use statement are only needed for define_cube_mosaic
use mpp_domains_mod, only : domain2d, mpp_define_mosaic, mpp_get_compute_domain, &
                            mpp_get_global_domain
use mosaic_mod, only : get_mosaic_ncontacts, get_mosaic_contact

implicit none;private

! ==== public interfaces =====================================================
! grid dimension inquiry subroutines
public :: get_grid_ntiles ! returns number of tiles
public :: get_grid_size   ! returns horizontal sizes of the grid
! grid geometry inquiry subroutines
public :: get_grid_cell_centers 
public :: get_grid_cell_vertices
! grid area inquiry subroutines
public :: get_grid_cell_area
public :: get_grid_comp_area
! decompose cubed sphere domains -- probably does not belong here, but it should 
! be in some place available for component models
public :: define_cube_mosaic
! ==== end of public interfaces ==============================================

interface get_grid_size
   module procedure get_grid_size_for_all_tiles
   module procedure get_grid_size_for_one_tile
end interface

interface get_grid_cell_vertices
   module procedure get_grid_cell_vertices_1D
   module procedure get_grid_cell_vertices_2D
end interface

interface get_grid_cell_centers
   module procedure get_grid_cell_centers_1D
   module procedure get_grid_cell_centers_2D
end interface

! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'grid_mod', &
     version     = '$Id$', &
     tagname     = '$Name$'

character(len=*), parameter :: &
     grid_dir  = 'INPUT/',     &      ! root directory for all grid files
     grid_file = 'INPUT/grid_spec.nc' ! name of the grid spec file

integer, parameter :: &
     MAX_NAME = 256,  & ! max length of the variable names
     MAX_FILE = 1024, & ! max length of the file names
     VERSION_0 = 0,   &
     VERSION_1 = 1,   &
     VERSION_2 = 2

integer, parameter :: BUFSIZE = 1048576  ! This is used to control memory usage in get_grid_comp_area
                                         ! We may change this to a namelist variable is needed.

! ==== module variables ======================================================
integer :: grid_version = -1
logical :: great_circle_algorithm = .FALSE.
logical :: first_call = .TRUE.


contains 

function get_grid_version()
  integer :: get_grid_version

  if(first_call) then
     great_circle_algorithm = get_great_circle_algorithm()
     first_call = .FALSE.
  endif

  if(grid_version<0) then
    if(field_exist(grid_file, 'geolon_t')) then
       grid_version = VERSION_0 
    else if(field_exist(grid_file, 'x_T')) then
       grid_version = VERSION_1
    else if(field_exist(grid_file, 'ocn_mosaic_file') ) then
       grid_version = VERSION_2
    else
       call error_mesg(module_name//'/get_grid_version',&
            'Can''t determine the version of the grid spec: none of "x_T", "geolon_t", or "ocn_mosaic_file" exist in file "'//trim(grid_file)//'"', &
            FATAL )
    endif
  endif
  get_grid_version = grid_version
end function get_grid_version


! ============================================================================
! returns number of tiles for a given component
! ============================================================================
subroutine get_grid_ntiles(component,ntiles)
  character(len=*)     :: component
  integer, intent(out) :: ntiles

  ! local vars
  character(len=MAX_FILE) :: component_mosaic

  select case (get_grid_version())
  case(VERSION_0,VERSION_1)
     ntiles = 1
  case(VERSION_2)
     call read_data(grid_file,trim(lowercase(component))//'_mosaic_file',component_mosaic)
     ntiles = get_mosaic_ntiles(grid_dir//trim(component_mosaic))
  end select
end subroutine get_grid_ntiles


! ============================================================================
! returns size of the grid for each of the tiles
! ============================================================================
subroutine get_grid_size_for_all_tiles(component,nx,ny)
  character(len=*)     :: component
  integer, intent(inout) :: nx(:),ny(:)

  ! local vars
  integer :: siz(4) ! for the size of external fields
  character(len=MAX_NAME) :: varname1, varname2
  character(len=MAX_FILE) :: component_mosaic
  
  varname1 = 'AREA_'//trim(uppercase(component))
  varname2 = trim(lowercase(component))//'_mosaic_file'

  select case (get_grid_version())
  case(VERSION_0,VERSION_1)
     call field_size(grid_file, varname1, siz)
     nx(1) = siz(1); ny(1)=siz(2)
  case(VERSION_2) ! mosaic file
     call read_data(grid_file,varname2, component_mosaic)
     call get_mosaic_grid_sizes(grid_dir//trim(component_mosaic),nx,ny)
  end select
end subroutine get_grid_size_for_all_tiles


! ============================================================================
! returns size of the grid for one of the tiles
! ============================================================================
subroutine get_grid_size_for_one_tile(component,tile,nx,ny)
  character(len=*)       :: component
  integer, intent(in)    :: tile
  integer, intent(inout) :: nx,ny
  
  ! local vars
  integer, allocatable :: nnx(:), nny(:)
  integer :: ntiles

  call get_grid_ntiles(component, ntiles)
  if(tile>0.and.tile<=ntiles) then
     allocate(nnx(ntiles),nny(ntiles))
     call get_grid_size_for_all_tiles(component,nnx,nny)
     nx = nnx(tile); ny = nny(tile)
     deallocate(nnx,nny)
  else
     call error_mesg('get_grid_size',&
          'requested tile index '//trim(string(tile))//' is out of bounds (1:'//trim(string(ntiles))//')',&
          FATAL)
  endif
end subroutine get_grid_size_for_one_tile

! ============================================================================
! return grid cell area for the specified model component and tile
! ============================================================================
subroutine get_grid_cell_area(component, tile, cellarea, domain)
  character(len=*), intent(in)    :: component
  integer         , intent(in)    :: tile
  real            , intent(inout) :: cellarea(:,:)
  type(domain2d)  , intent(in), optional :: domain

  ! local vars
  integer :: nlon, nlat
  real, allocatable :: glonb(:,:), glatb(:,:)

  select case(get_grid_version())
  case(VERSION_0,VERSION_1)
     select case(trim(component))
     case('LND')
        call read_data(grid_file, 'AREA_LND_CELL', cellarea, &
            no_domain=.not.present(domain), domain=domain)
     case('ATM','OCN')
        call read_data(grid_file, 'AREA_'//trim(uppercase(component)),cellarea,&
            no_domain=.not.present(domain),domain=domain)
     case default
        call error_mesg(module_name//'/get_grid_cell_area',&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select
     ! convert area to m2
     cellarea = cellarea*4.*PI*radius**2
  case(VERSION_2)
     if (present(domain)) then
        call mpp_get_compute_domain(domain,xsize=nlon,ysize=nlat)
     else
        call get_grid_size(component,tile,nlon,nlat)
     endif
     allocate(glonb(nlon+1,nlat+1),glatb(nlon+1,nlat+1))
     call get_grid_cell_vertices(component, tile, glonb, glatb, domain)
     if (great_circle_algorithm) then
        call calc_mosaic_grid_great_circle_area(glonb*pi/180.0, glatb*pi/180.0, cellarea)
     else
        call calc_mosaic_grid_area(glonb*pi/180.0, glatb*pi/180.0, cellarea)
     end if
     deallocate(glonb,glatb)
  end select

end subroutine get_grid_cell_area


! ============================================================================
! get the area of the component per grid cell
! ============================================================================
subroutine get_grid_comp_area(component,tile,area,domain)
  character(len=*) :: component
  integer, intent(in) :: tile
  real, intent(inout) :: area(:,:)
  type(domain2d), intent(in), optional :: domain
  ! local vars
  integer :: n_xgrid_files ! number of exchange grid files in the mosaic
  integer :: siz(4), nxgrid
  integer :: i,j,m,n
  integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
  real, allocatable :: xgrid_area(:)
  real, allocatable :: rmask(:,:)
  character(len=MAX_NAME) :: &
     xgrid_name, & ! name of the variable holding xgrid names
     tile_name,  & ! name of the tile
     xgrid_file, & ! name of the current xgrid file
     mosaic_name,& ! name of the mosaic
     mosaic_file,&
     tilefile 
  character(len=4096)     :: attvalue
  character(len=MAX_NAME), allocatable :: nest_tile_name(:)
  character(len=MAX_NAME) :: varname1, varname2
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0, j0 ! offsets for x and y, respectively
  integer :: num_nest_tile, ntiles
  logical :: is_nest
  integer :: found_xgrid_files ! how many xgrid files we actually found in the grid spec
  integer :: ibegin, iend, bsize, l

  select case (get_grid_version())
  case(VERSION_0,VERSION_1)
     select case(component)
     case('ATM')
        call read_data(grid_file,'AREA_ATM',area, no_domain=.not.present(domain),domain=domain)
     case('OCN')
        allocate(rmask(size(area,1),size(area,2)))
        call read_data(grid_file,'AREA_OCN',area, no_domain=.not.present(domain),domain=domain)
        call read_data(grid_file,'wet',     rmask,no_domain=.not.present(domain),domain=domain)
        area = area*rmask
        deallocate(rmask)
     case('LND')
        call read_data(grid_file,'AREA_LND',area,no_domain=.not.present(domain),domain=domain)
     case default
        call error_mesg(module_name//'/get_grid_comp_area',&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select
  case(VERSION_2) ! mosaic gridspec
     select case (component)
     case ('ATM')
        ! just read the grid cell area and return
        call get_grid_cell_area(component,tile,area)
        return
     case ('LND')
        xgrid_name = 'aXl_file'
        call read_data(grid_file, 'lnd_mosaic', mosaic_name)
        tile_name  = trim(mosaic_name)//'_tile'//char(tile+ichar('0'))
     case ('OCN')
        xgrid_name = 'aXo_file'
        call read_data(grid_file, 'ocn_mosaic', mosaic_name)
        tile_name  = trim(mosaic_name)//'_tile'//char(tile+ichar('0'))
     case default
        call error_mesg(module_name//'/get_grid_comp_area',&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
             FATAL)
     end select
     ! get the boundaries of the requested domain
     if(present(domain)) then
        call mpp_get_compute_domain(domain,is,ie,js,je)
        i0 = 1-is ; j0=1-js
     else
        call get_grid_size(component,tile,ie,je)
        is = 1 ; i0 = 0
        js = 1 ; j0 = 0
     endif
     if (size(area,1)/=ie-is+1.or.size(area,2)/=je-js+1) &
        call error_mesg(module_name//'/get_grid_comp_area',&
        'size of the output argument "area" is not consistent with the domain',FATAL) 

     ! find the nest tile 
     call read_data(grid_file, 'atm_mosaic', mosaic_name)
     call read_data(grid_file,'atm_mosaic_file',mosaic_file)
     mosaic_file = grid_dir//trim(mosaic_file)
     ntiles = get_mosaic_ntiles(trim(mosaic_file))   
     allocate(nest_tile_name(ntiles))  
     num_nest_tile = 0
     do n = 1, ntiles
        call read_data(mosaic_file, 'gridfiles', tilefile, level=n)        
        tilefile = grid_dir//trim(tilefile)
        if( get_global_att_value(tilefile, "nest_grid", attvalue) ) then
           if(trim(attvalue) == "TRUE") then
              num_nest_tile = num_nest_tile + 1
              nest_tile_name(num_nest_tile) = trim(mosaic_name)//'_tile'//char(n+ichar('0'))
           else if(trim(attvalue) .NE. "FALSE") then
              call error_mesg(module_name//'/get_grid_comp_area', 'value of global attribute nest_grid in file'// &
                   trim(tilefile)//' should be TRUE of FALSE', FATAL)
           endif
        end if
     end do
     area(:,:) = 0.
     if(field_exist(grid_file,xgrid_name)) then
        ! get the number of the exchange-grid files
        call field_size(grid_file,xgrid_name,siz)
        n_xgrid_files = siz(2)
        found_xgrid_files = 0
        ! loop through all exchange grid files
        do n = 1, n_xgrid_files
           ! get the name of the current exchange grid file
           call read_data(grid_file,xgrid_name,xgrid_file,level=n)
           ! skip the rest of the loop if the name of the current tile isn't found 
           ! in the file name, but check this only if there is more than 1 tile
           if(n_xgrid_files>1) then
              if(index(xgrid_file,trim(tile_name))==0) cycle
           endif
           found_xgrid_files = found_xgrid_files + 1
           !---make sure the atmosphere grid is not a nested grid
           is_nest = .false. 
           do m = 1, num_nest_tile
              if(index(xgrid_file, trim(nest_tile_name(m))) .NE. 0) then
                 is_nest = .true.
                 exit
              end if
           end do
           if(is_nest) cycle 

           ! finally read the exchange grid
           nxgrid = get_mosaic_xgrid_size(grid_dir//xgrid_file)
           if(nxgrid < BUFSIZE) then
              allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), xgrid_area(nxgrid))
           else
              allocate(i1(BUFSIZE), j1(BUFSIZE), i2(BUFSIZE), j2(BUFSIZE), xgrid_area(BUFSIZE))
           endif
           ibegin = 1
           do l = 1,nxgrid,BUFSIZE
              bsize = min(BUFSIZE, nxgrid-l+1)
              iend = ibegin + bsize - 1
              call get_mosaic_xgrid(grid_dir//xgrid_file, i1(1:bsize), j1(1:bsize), i2(1:bsize), j2(1:bsize), &
                                    xgrid_area(1:bsize), ibegin, iend)
              ! and sum the exchange grid areas
              do m = 1, bsize
                 i = i2(m); j = j2(m)
                 if (i<is.or.i>ie) cycle
                 if (j<js.or.j>je) cycle
                 area(i+i0,j+j0) = area(i+i0,j+j0) + xgrid_area(m)
              end do
              ibegin = iend + 1
           enddo
           deallocate(i1, j1, i2, j2, xgrid_area)
        enddo
        if (found_xgrid_files == 0) &
           call error_mesg('get_grid_comp_area', 'no xgrid files were found for component '& 
                 //trim(component)//' (mosaic name is '//trim(mosaic_name)//')', FATAL)

     endif
     deallocate(nest_tile_name)
  end select ! version
  ! convert area to m2
  area = area*4.*PI*radius**2
end subroutine

! ============================================================================
! returns arrays of global grid cell boundaries for given model component and 
! mosaic tile number.
! NOTE that in case of non-lat-lon grid the returned coordinates may have be not so 
! meaningful, by the very nature of such grids. But presumably these 1D coordinate 
! arrays are good enough for diag axis and such.
! ============================================================================
subroutine get_grid_cell_vertices_1D(component, tile, glonb, glatb)
  character(len=*), intent(in) :: component
  integer,          intent(in) :: tile
  real,          intent(inout) :: glonb(:),glatb(:)

  integer                      :: nlon, nlat
  integer                      :: start(4), nread(4)
  real, allocatable            :: tmp(:,:), x_vert_t(:,:,:), y_vert_t(:,:,:)
  character(len=MAX_FILE)      :: filename1, filename2

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)
  if (size(glonb(:))/=nlon+1) &
       call error_mesg ( module_name//'/get_grid_cell_vertices_1D',&
       'Size of argument "glonb" is not consistent with the grid size',FATAL)
  if (size(glatb(:))/=nlat+1) &
       call error_mesg ( module_name//'/get_grid_cell_vertices_1D',&
       'Size of argument "glatb" is not consistent with the grid size',FATAL)
  if(trim(component) .NE. 'ATM' .AND. component .NE. 'LND' .AND. component .NE. 'OCN') then
     call error_mesg(module_name//'/get_grid_cell_vertices_1D',&
          'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
          FATAL)     
  endif

  select case(get_grid_version())
  case(VERSION_0)
     select case(trim(component))
     case('ATM','LND')
        call read_data(grid_file, 'xb'//lowercase(component(1:1)), glonb, no_domain=.true.)
        call read_data(grid_file, 'yb'//lowercase(component(1:1)), glatb, no_domain=.true.)
     case('OCN')
        call read_data(grid_file, "gridlon_vert_t", glonb, no_domain=.true.) 
        call read_data(grid_file, "gridlat_vert_t", glatb, no_domain=.true.) 
     end select
  case(VERSION_1)
     select case(trim(component))
     case('ATM','LND')
        call read_data(grid_file, 'xb'//lowercase(component(1:1)), glonb, no_domain=.true.)
        call read_data(grid_file, 'yb'//lowercase(component(1:1)), glatb, no_domain=.true.)
     case('OCN')
        allocate (x_vert_t(nlon,1,2), y_vert_t(1,nlat,2) ) 
        start = 1; nread = 1
        nread(1) = nlon; nread(2) = 1; start(3) = 1
        call read_data(grid_file, "x_vert_T", x_vert_t(:,:,1), start, nread, no_domain=.TRUE.)
        nread(1) = nlon; nread(2) = 1; start(3) = 2
        call read_data(grid_file, "x_vert_T", x_vert_t(:,:,2), start, nread, no_domain=.TRUE.)

        nread(1) = 1; nread(2) = nlat; start(3) = 1
        call read_data(grid_file, "y_vert_T", y_vert_t(:,:,1), start, nread, no_domain=.TRUE.)
        nread(1) = 1; nread(2) = nlat; start(3) = 4
        call read_data(grid_file, "y_vert_T", y_vert_t(:,:,2), start, nread, no_domain=.TRUE.)
        glonb(1:nlon) = x_vert_t(1:nlon,1,1)
        glonb(nlon+1) = x_vert_t(nlon,1,2)
        glatb(1:nlat) = y_vert_t(1,1:nlat,1)
        glatb(nlat+1) = y_vert_t(1,nlat,2)
        deallocate(x_vert_t, y_vert_t)
     end select
  case(VERSION_2)
     ! get the name of the mosaic file for the component
     call read_data(grid_file, trim(lowercase(component))//'_mosaic_file', filename1)
     filename1=grid_dir//trim(filename1)
     ! get the name of the grid file for the component and tile
     call read_data(filename1, 'gridfiles', filename2, level=tile)
     filename2 = grid_dir//trim(filename2)

     start = 1; nread = 1
     nread(1) = 2*nlon+1
     allocate( tmp(2*nlon+1,1) )
     call read_data(filename2, "x", tmp, start, nread, no_domain=.TRUE.)
     glonb(1:nlon+1) = tmp(1:2*nlon+1:2,1)
     deallocate(tmp)
     allocate(tmp(1,2*nlat+1))

     start = 1; nread = 1
     nread(2) = 2*nlat+1
     call read_data(filename2, "y", tmp, start, nread, no_domain=.TRUE.)
     glatb(1:nlat+1) = tmp(1,1:2*nlat+1:2)
     deallocate(tmp)
  end select

end subroutine get_grid_cell_vertices_1D

! ============================================================================
! returns cell vertices for the specified model component and mosaic tile number
! ============================================================================
subroutine get_grid_cell_vertices_2D(component, tile, lonb, latb, domain)
  character(len=*),         intent(in) :: component
  integer,                  intent(in) :: tile
  real,                  intent(inout) :: lonb(:,:),latb(:,:)
  type(domain2d), optional, intent(in) :: domain

  ! local vars
  character(len=MAX_FILE) :: filename1, filename2
  integer :: nlon, nlat
  integer :: i,j
  real, allocatable :: buffer(:), tmp(:,:), x_vert_t(:,:,:), y_vert_t(:,:,:)
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0,j0 ! offsets for coordinates
  integer :: isg, jsg
  integer :: start(4), nread(4)

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)
  if (present(domain)) then
    call mpp_get_compute_domain(domain,is,ie,js,je)
  else
    is = 1 ; ie = nlon
    js = 1 ; je = nlat
    !--- domain normally should be present
    call error_mesg ( module_name//'/get_grid_cell_vertices',&
       'domain is not present, global data will be read', NOTE)
  endif
  i0 = -is+1; j0 = -js+1
  
  ! verify that lonb and latb sizes are consistent with the size of domain
  if (size(lonb,1)/=ie-is+2.or.size(lonb,2)/=je-js+2) &
       call error_mesg ( module_name//'/get_grid_cell_vertices',&
       'Size of argument "lonb" is not consistent with the domain size',FATAL)
  if (size(latb,1)/=ie-is+2.or.size(latb,2)/=je-js+2) &
       call error_mesg ( module_name//'/get_grid_cell_vertices',&
       'Size of argument "latb" is not consistent with the domain size',FATAL)
  if(trim(component) .NE. 'ATM' .AND. component .NE. 'LND' .AND. component .NE. 'OCN') then
     call error_mesg(module_name//'/get_grid_cell_vertices',&
          'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
          FATAL)  
  endif

  select case(get_grid_version())
  case(VERSION_0)
     select case(component)
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)+1))
        ! read coordinates of grid cell vertices
        call read_data(grid_file, 'xb'//lowercase(component(1:1)), buffer(1:nlon+1), no_domain=.true.)
        do j = js, je+1
           do i = is, ie+1
              lonb(i+i0,j+j0) = buffer(i)
           enddo
        enddo
        call read_data(grid_file, 'yb'//lowercase(component(1:1)), buffer(1:nlat+1), no_domain=.true.)
        do j = js, je+1
           do i = is, ie+1
              latb(i+i0,j+j0) = buffer(j)
           enddo
        enddo
        deallocate(buffer)
     case('OCN')
        if (present(domain)) then
           start = 1; nread = 1
           start(1) = is; start(2) = js
           nread(1) = ie-is+2; nread(2) = je-js+2
           call read_data(grid_file, 'geolon_vert_t', lonb, start, nread, no_domain=.true. )
           call read_data(grid_file, 'geolat_vert_t', latb, start, nread, no_domain=.true. )
         else
           call read_data(grid_file, 'geolon_vert_t', lonb, no_domain=.TRUE. )
           call read_data(grid_file, 'geolat_vert_t', latb, no_domain=.TRUE. )
         endif
     end select
  case(VERSION_1)
     select case(component)
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)+1))
        ! read coordinates of grid cell vertices
        call read_data(grid_file, 'xb'//lowercase(component(1:1)), buffer(1:nlon+1), no_domain=.true.)
        do j = js, je+1
           do i = is, ie+1
              lonb(i+i0,j+j0) = buffer(i)
           enddo
        enddo
        call read_data(grid_file, 'yb'//lowercase(component(1:1)), buffer(1:nlat+1), no_domain=.true.)
        do j = js, je+1
           do i = is, ie+1
              latb(i+i0,j+j0) = buffer(j)
           enddo
        enddo
        deallocate(buffer)
     case('OCN')
        nlon=ie-is+1; nlat=je-js+1
        allocate (x_vert_t(nlon,nlat,4), y_vert_t(nlon,nlat,4) ) 
        call read_data(grid_file, 'x_vert_T', x_vert_t, no_domain=.not.present(domain), domain=domain )
        call read_data(grid_file, 'y_vert_T', y_vert_t, no_domain=.not.present(domain), domain=domain )
        lonb(1:nlon,1:nlat) = x_vert_t(1:nlon,1:nlat,1)
        lonb(nlon+1,1:nlat) = x_vert_t(nlon,1:nlat,2)
        lonb(1:nlon,nlat+1) = x_vert_t(1:nlon,nlat,4)
        lonb(nlon+1,nlat+1) = x_vert_t(nlon,nlat,3)
        latb(1:nlon,1:nlat) = y_vert_t(1:nlon,1:nlat,1)
        latb(nlon+1,1:nlat) = y_vert_t(nlon,1:nlat,2)
        latb(1:nlon,nlat+1) = y_vert_t(1:nlon,nlat,4)
        latb(nlon+1,nlat+1) = y_vert_t(nlon,nlat,3)
        deallocate(x_vert_t, y_vert_t)
     end select
  case(VERSION_2)
     ! get the name of the mosaic file for the component
     call read_data(grid_file, trim(lowercase(component))//'_mosaic_file', filename1)
     filename1=grid_dir//trim(filename1)
     ! get the name of the grid file for the component and tile
     call read_data(filename1, 'gridfiles', filename2, level=tile)
     filename2 = grid_dir//trim(filename2)
     if(PRESENT(domain)) then
        call mpp_get_global_domain(domain, xbegin=isg, ybegin=jsg)
        start = 1; nread = 1
        start(1) = 2*(is-isg+1) - 1; nread(1) = 2*(ie-is)+3
        start(2) = 2*(js-jsg+1) - 1; nread(2) = 2*(je-js)+3
        allocate(tmp(nread(1), nread(2)) )
        call read_data(filename2, 'x', tmp, start, nread, no_domain=.TRUE.)
        do j = 1, je-js+2
           do i = 1, ie-is+2
              lonb(i,j) = tmp(2*i-1,2*j-1)
           enddo
        enddo
        call read_data(filename2, 'y', tmp, start, nread, no_domain=.TRUE.)
        do j = 1, je-js+2
           do i = 1, ie-is+2
              latb(i,j) = tmp(2*i-1,2*j-1)
           enddo
        enddo        
     else
        allocate(tmp(2*nlon+1,2*nlat+1))
        call read_data(filename2, 'x', tmp, no_domain=.TRUE.)
        do j = js, je+1
           do i = is, ie+1
              lonb(i+i0,j+j0) = tmp(2*i-1,2*j-1)
           end do
        end do
        call read_data(filename2, 'y', tmp, no_domain=.TRUE.)
        do j = js, je+1
           do i = is, ie+1
              latb(i+i0,j+j0) = tmp(2*i-1,2*j-1)
           end do
        end do
     endif
     deallocate(tmp)
  end select

end subroutine get_grid_cell_vertices_2D

! ============================================================================
! returns global coordinate arrays fro given model component and mosaic tile number
! NOTE that in case of non-lat-lon grid those coordinates may have be not so 
! meaningful, by the very nature of such grids. But presumably these 1D coordinate 
! arrays are good enough for diag axis and such.
! ============================================================================
subroutine get_grid_cell_centers_1D(component, tile, glon, glat)
  character(len=*), intent(in) :: component
  integer, intent(in) :: tile
  real, intent(inout) :: glon(:),glat(:)
  integer                      :: nlon, nlat
  integer                      :: start(4), nread(4)
  real, allocatable            :: tmp(:,:)
  character(len=MAX_FILE)      :: filename1, filename2

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)
  if (size(glon(:))/=nlon) &
       call error_mesg ( module_name//'/get_grid_cell_centers_1D',&
       'Size of argument "glon" is not consistent with the grid size',FATAL)
  if (size(glat(:))/=nlat) &
       call error_mesg ( module_name//'/get_grid_cell_centers_1D',&
       'Size of argument "glat" is not consistent with the grid size',FATAL)
  if(trim(component) .NE. 'ATM' .AND. component .NE. 'LND' .AND. component .NE. 'OCN') then
     call error_mesg(module_name//'/get_grid_cell_centers_1D',&
          'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
          FATAL)     
  endif

  select case(get_grid_version())
  case(VERSION_0)
     select case(trim(component))
     case('ATM','LND')
        call read_data(grid_file, 'xt'//lowercase(component(1:1)), glon, no_domain=.true.)
        call read_data(grid_file, 'yt'//lowercase(component(1:1)), glat, no_domain=.true.)
     case('OCN')
        call read_data(grid_file, "gridlon_t", glon, no_domain=.true.) 
        call read_data(grid_file, "gridlat_t", glat, no_domain=.true.)
     end select 
  case(VERSION_1)
     select case(trim(component))
     case('ATM','LND')
        call read_data(grid_file, 'xt'//lowercase(component(1:1)), glon, no_domain=.true.)
        call read_data(grid_file, 'yt'//lowercase(component(1:1)), glat, no_domain=.true.)
     case('OCN')
        call read_data(grid_file, "grid_x_T", glon, no_domain=.true.) 
        call read_data(grid_file, "grid_y_T", glat, no_domain=.true.) 
     end select
  case(VERSION_2)
     ! get the name of the mosaic file for the component
     call read_data(grid_file, trim(lowercase(component))//'_mosaic_file', filename1)
     filename1=grid_dir//trim(filename1)
     ! get the name of the grid file for the component and tile
     call read_data(filename1, 'gridfiles', filename2, level=tile)
     filename2 = grid_dir//trim(filename2)

     start = 1; nread = 1
     nread(1) = 2*nlon+1; start(2) = 2
     allocate( tmp(2*nlon+1,1) )
     call read_data(filename2, "x", tmp, start, nread, no_domain=.TRUE.)
     glon(1:nlon) = tmp(2:2*nlon:2,1)
     deallocate(tmp)
     allocate(tmp(1, 2*nlat+1))

     start = 1; nread = 1
     nread(2) = 2*nlat+1; start(1) = 2
     call read_data(filename2, "y", tmp, start, nread, no_domain=.TRUE.)
     glat(1:nlat) = tmp(1,2:2*nlat:2)
     deallocate(tmp)
  end select

  
end subroutine get_grid_cell_centers_1D

! ============================================================================
! returns grid cell centers for specified model component and mosaic tile number
! ============================================================================
subroutine get_grid_cell_centers_2D(component, tile, lon, lat, domain)
  character(len=*), intent(in) :: component
  integer, intent(in) :: tile
  real, intent(inout) :: lon(:,:),lat(:,:)
  type(domain2d), intent(in), optional :: domain
  ! local vars
  character(len=MAX_NAME) :: varname
  character(len=MAX_FILE) :: filename1, filename2
  integer :: nlon, nlat
  integer :: i,j
  real, allocatable :: buffer(:),tmp(:,:)
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0,j0 ! offsets for coordinates
  integer :: isg, jsg
  integer :: start(4), nread(4) 

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)
  if (present(domain)) then
    call mpp_get_compute_domain(domain,is,ie,js,je)
  else
    is = 1 ; ie = nlon
    js = 1 ; je = nlat
    !--- domain normally should be present
    call error_mesg ( module_name//'/get_grid_cell_centers',&
       'domain is not present, global data will be read', NOTE)
  endif
  i0 = -is+1; j0 = -js+1

  ! verify that lon and lat sizes are consistent with the size of domain
  if (size(lon,1)/=ie-is+1.or.size(lon,2)/=je-js+1) &
       call error_mesg ( module_name//'/get_grid_cell_centers',&
       'Size of array "lon" is not consistent with the domain size',&
       FATAL )
  if (size(lat,1)/=ie-is+1.or.size(lat,2)/=je-js+1) &
       call error_mesg ( module_name//'/get_grid_cell_centers',&
       'Size of array "lat" is not consistent with the domain size',&
       FATAL )
  if(trim(component) .NE. 'ATM' .AND. component .NE. 'LND' .AND. component .NE. 'OCN') then
     call error_mesg(module_name//'/get_grid_cell_vertices',&
          'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN',&
          FATAL)  
  endif

  select case(get_grid_version())
  case(VERSION_0)
     select case (trim(component))
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)))
        ! read coordinates of grid cell vertices
        call read_data(grid_file, 'xt'//lowercase(component(1:1)), buffer(1:nlon), no_domain=.true.)
        do j = js,je
        do i = is,ie
           lon(i+i0,j+j0) = buffer(i)
        enddo
        enddo
        call read_data(grid_file, 'yt'//lowercase(component(1:1)), buffer(1:nlat), no_domain=.true.)
        do j = js,je
        do i = is,ie
           lat(i+i0,j+j0) = buffer(j)
        enddo
        enddo
        deallocate(buffer)
     case('OCN')
        call read_data(grid_file, 'geolon_t', lon, no_domain=.not.present(domain), domain=domain )
        call read_data(grid_file, 'geolat_t', lat, no_domain=.not.present(domain), domain=domain )
     end select
  case(VERSION_1)
     select case(trim(component))
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)))
        ! read coordinates of grid cell vertices
        call read_data(grid_file, 'xt'//lowercase(component(1:1)), buffer(1:nlon), no_domain=.true.)
        do j = js,je
        do i = is,ie
           lon(i+i0,j+j0) = buffer(i)
        enddo
        enddo
        call read_data(grid_file, 'yt'//lowercase(component(1:1)), buffer(1:nlat), no_domain=.true.)
        do j = js,je
        do i = is,ie
           lat(i+i0,j+j0) = buffer(j)
        enddo
        enddo
        deallocate(buffer)
     case('OCN')
        call read_data(grid_file, 'x_T', lon, no_domain=.not.present(domain), domain=domain )
        call read_data(grid_file, 'y_T', lat, no_domain=.not.present(domain), domain=domain )
     end select
  case(VERSION_2) ! mosaic grid file
     ! get the name of the mosaic file for the component
     call read_data(grid_file, trim(lowercase(component))//'_mosaic_file', filename1)
     filename1=grid_dir//trim(filename1)
     ! get the name of the grid file for the component and tile
     call read_data(filename1, 'gridfiles', filename2, level=tile)
     filename2 = grid_dir//trim(filename2)
     if(PRESENT(domain)) then
        call mpp_get_global_domain(domain, xbegin=isg, ybegin=jsg)
        start = 1; nread = 1
        start(1) = 2*(is-isg+1) - 1; nread(1) = 2*(ie-is)+3
        start(2) = 2*(js-jsg+1) - 1; nread(2) = 2*(je-js)+3
        allocate(tmp(nread(1), nread(2)))
        call read_data(filename2, 'x', tmp, start, nread, no_domain=.TRUE.)
        do j = 1, je-js+1
           do i = 1, ie-is+1
              lon(i,j) = tmp(2*i,2*j)
           enddo
        enddo
        call read_data(filename2, 'y', tmp, start, nread, no_domain=.TRUE.)
        do j = 1, je-js+1
           do i = 1, ie-is+1
              lat(i,j) = tmp(2*i,2*j)
           enddo
        enddo        
     else
        allocate(tmp(2*nlon+1,2*nlat+1))
        call read_data(filename2, 'x', tmp, no_domain=.TRUE.)
        do j = js,je
           do i = is,ie
              lon(i+i0,j+j0) = tmp(2*i,2*j)
           end do
        end do
        call read_data(filename2, 'y', tmp, no_domain=.TRUE.)
        do j = js,je
           do i = is,ie
              lat(i+i0,j+j0) = tmp(2*i,2*j)
           end do
        end do
        deallocate(tmp)
     endif
  end select

end subroutine get_grid_cell_centers_2D


! ============================================================================
! given a model component, a layout, and (optionally) a halo size, returns a 
! domain for current processor
! ============================================================================
! this subroutine probably does not belong in the grid_mod 
subroutine define_cube_mosaic ( component, domain, layout, halo, maskmap )
  character(len=*) , intent(in)    :: component
  type(domain2d)   , intent(inout) :: domain
  integer          , intent(in)    :: layout(2)
  integer, optional, intent(in)    :: halo 
  logical, optional, intent(in)    :: maskmap(:,:,:)

  ! ---- local constants
  
  ! ---- local vars
  character(len=MAX_NAME) :: varname
  character(len=MAX_FILE) :: mosaic_file
  integer :: ntiles     ! number of tiles
  integer :: ncontacts  ! number of contacts between mosaic tiles
  integer :: n
  integer :: ng, pe_pos, npes         ! halo size
  integer, allocatable :: nlon(:), nlat(:), global_indices(:,:)
  integer, allocatable :: pe_start(:), pe_end(:), layout_2d(:,:)
  integer, allocatable :: tile1(:),tile2(:)
  integer, allocatable :: is1(:),ie1(:),js1(:),je1(:)
  integer, allocatable :: is2(:),ie2(:),js2(:),je2(:)

  call get_grid_ntiles(component,ntiles)
  allocate(nlon(ntiles), nlat(ntiles))
  allocate(global_indices(4,ntiles))
  allocate(pe_start(ntiles),pe_end(ntiles))
  allocate(layout_2d(2,ntiles))
  call get_grid_size(component,nlon,nlat)

  pe_pos = mpp_root_pe()
  do n = 1, ntiles
     global_indices(:,n) = (/ 1, nlon(n), 1, nlat(n) /)
     layout_2d     (:,n) = layout
     if(present(maskmap)) then
        npes = count(maskmap(:,:,n))
     else
        npes = layout(1)*layout(2)
     endif
     pe_start(n) = pe_pos
     pe_end  (n) = pe_pos + npes - 1
     pe_pos      = pe_end(n) + 1
  enddo

  varname=trim(lowercase(component))//'_mosaic_file'
  call read_data(grid_file,varname,mosaic_file)
  mosaic_file = grid_dir//mosaic_file

  ! get the contact information from mosaic file
  ncontacts = get_mosaic_ncontacts(mosaic_file)
  allocate(tile1(ncontacts),tile2(ncontacts))
  allocate(is1(ncontacts),ie1(ncontacts),js1(ncontacts),je1(ncontacts))
  allocate(is2(ncontacts),ie2(ncontacts),js2(ncontacts),je2(ncontacts))
  call get_mosaic_contact(mosaic_file, tile1, tile2, &
       is1, ie1, js1, je1, is2, ie2, js2, je2)

  ng = 0
  if(present(halo)) ng = halo
  ! create the domain2d variable
  call mpp_define_mosaic ( global_indices, layout_2d, domain, &
       ntiles, ncontacts, tile1, tile2,                  &
       is1, ie1, js1, je1, &
       is2, ie2, js2, je2, &
       pe_start=pe_start, pe_end=pe_end, symmetry=.true.,  &
       shalo = ng, nhalo = ng, whalo = ng, ehalo = ng,     &
       maskmap = maskmap,                                  &
       name = trim(component)//'Cubic-Sphere Grid' )

  deallocate(nlon,nlat,global_indices,pe_start,pe_end,layout_2d)
  deallocate(tile1,tile2)
  deallocate(is1,ie1,js1,je1)
  deallocate(is2,ie2,js2,je2)

end subroutine define_cube_mosaic

end module grid_mod
