module mosaic_mod

! <CONTACT EMAIL="Zhi.Liang@noaa.gov">
!   Zhi Liang
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!    <TT>mosaic_mod</TT> implements some utility routines to read mosaic information.
! </OVERVIEW>

! <DESCRIPTION>
!    <TT>mosaic_mod</TT> implements some utility routines to read mosaic information.
!    The information includes number of tiles and contacts in the mosaic, 
!    mosaic grid resolution of each tile, mosaic contact information, mosaic exchange
!    grid information. Each routine will call a C-version routine to get these information.
! </DESCRIPTION>

use mpp_mod,    only : mpp_error, FATAL, mpp_pe, mpp_root_pe
use mpp_io_mod, only : MPP_MULTI
use fms_io_mod, only : dimension_size, field_exist, read_data, read_compressed
use constants_mod, only : PI, RADIUS

implicit none
private

character(len=*), parameter :: &
     grid_dir  = 'INPUT/'      ! root directory for all grid files

integer, parameter :: &
     MAX_NAME = 256,  & ! max length of the variable names
     MAX_FILE = 1024, & ! max length of the file names
     X_REFINE = 2,    & ! supergrid size/model grid size in x-direction
     Y_REFINE = 2       ! supergrid size/model grid size in y-direction

! --- public interface


public :: get_mosaic_ntiles
public :: get_mosaic_ncontacts
public :: get_mosaic_grid_sizes
public :: get_mosaic_contact
public :: get_mosaic_xgrid_size
public :: get_mosaic_xgrid
public :: calc_mosaic_grid_area
public :: calc_mosaic_grid_great_circle_area
public :: is_inside_polygon

logical :: module_is_initialized = .true.
! version information varaible
 character(len=128) :: version = '$Id$'
 character(len=128) :: tagname = '$Name$'

contains

!#######################################################################

! <SUBROUTINE NAME="mosaic_init">
!   <OVERVIEW>
!     Initialize the mosaic_mod. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initialization routine for the mosaic module. It writes the 
!     version information to the log file.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call mosaic_init ( )
!   </TEMPLATE>
subroutine mosaic_init() 

  if (module_is_initialized) return
  module_is_initialized = .TRUE.

!--------- write version number and namelist ------------------
!  call write_version_number (version, tagname)

end subroutine mosaic_init
! </SUBROUTINE>

!#######################################################################
! <FUNCTION NAME="get_mosaic_xgrid_size">
!   <OVERVIEW>
!     return exchange grid size of mosaic xgrid file.
!   </OVERVIEW>
!   <DESCRIPTION>
!     return exchange grid size of mosaic xgrid file.
!   </DESCRIPTION>
!   <TEMPLATE>
!    nxgrid = get_mosaic_xgrid_size(xgrid_file)
!   </TEMPLATE>
!   <IN NAME="xgrid_file" TYPE="character(len=*)">
!     The file that contains exchange grid information.
!   </IN>
  function get_mosaic_xgrid_size(xgrid_file)
    character(len=*), intent(in)          :: xgrid_file
    integer                               :: get_mosaic_xgrid_size

    get_mosaic_xgrid_size = dimension_size(xgrid_file, "ncells", no_domain=.TRUE.)

    return   

  end function get_mosaic_xgrid_size
! </FUNCTION>
!#######################################################################
! <SUBROUTINE NAME="get_mosaic_xgrid">
!   <OVERVIEW>
!     get exchange grid information from mosaic xgrid file.
!   </OVERVIEW>
!   <DESCRIPTION>
!     get exchange grid information from mosaic xgrid file.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_mosaic_xgrid(xgrid_file, nxgrid, i1, j1, i2, j2, area)
!   </TEMPLATE>
!   <IN NAME="xgrid_file" TYPE="character(len=*)">
!     The file that contains exchange grid information.
!   </IN>
!   <INOUT NAME="nxgrid" TYPE="integer">
!     number of exchange grid in xgrid_file
!   </INOUT>
!   <INOUT NAME="i1, j1" TYPE="integer, dimension(:)">
!     i and j-index in grid 1 of exchange grid.
!   </INOUT>
!   <INOUT NAME="i2, j2" TYPE="integer, dimension(:)">
!     i and j-index in grid 2 of exchange grid.
!   </INOUT>
!   <INOUT NAME="area" TYPE="real, dimension(:)">
!     area of the exchange grid. The area is scaled to represent unit earth area.
!   </INOUT>
  subroutine get_mosaic_xgrid(xgrid_file, i1, j1, i2, j2, area, ibegin, iend)
    character(len=*), intent(in) :: xgrid_file
    integer,       intent(inout) :: i1(:), j1(:), i2(:), j2(:)
    real,          intent(inout) :: area(:)
    integer, optional, intent(in) :: ibegin, iend

    integer                            :: start(4), nread(4), istart
    real,    dimension(2, size(i1(:))) :: tile1_cell, tile2_cell
    integer                            :: nxgrid, n
    real                               :: garea
    real                               :: get_global_area;

    garea = get_global_area();

    ! When start and nread present, make sure nread(1) is the same as the size of the data
    if(present(ibegin) .and. present(iend)) then
       istart = ibegin
       nxgrid = iend - ibegin + 1
       if(nxgrid .NE. size(i1(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(i1(:))")
       if(nxgrid .NE. size(j1(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(j1(:))")
       if(nxgrid .NE. size(i2(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(i2(:))")
       if(nxgrid .NE. size(j2(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(j2(:))")
       if(nxgrid .NE. size(area(:))) call mpp_error(FATAL, "get_mosaic_xgrid: nxgrid .NE. size(area(:))")
    else
       istart = 1
       nxgrid = size(i1(:))
    endif

    start  = 1; nread = 1
    start(1) = istart; nread(1) = nxgrid
    call read_compressed(xgrid_file, 'xgrid_area', area, start=start, nread=nread, threading=MPP_MULTI)
    start = 1; nread = 1
    nread(1) = 2
    start(2) = istart; nread(2) = nxgrid
    call read_compressed(xgrid_file, 'tile1_cell', tile1_cell, start=start, nread=nread, threading=MPP_MULTI)
    call read_compressed(xgrid_file, 'tile2_cell', tile2_cell, start=start, nread=nread, threading=MPP_MULTI)

     do n = 1, nxgrid
       i1(n) = tile1_cell(1,n) 
       j1(n) = tile1_cell(2,n)
       i2(n) = tile2_cell(1,n) 
       j2(n) = tile2_cell(2,n)
       area(n) = area(n)/garea
    end do
    
    return

  end subroutine get_mosaic_xgrid
! </SUBROUTINE>

  !###############################################################################
  ! <SUBROUTINE NAME="get_mosaic_ntiles">
  !   <OVERVIEW>
  !     get number of tiles in the mosaic_file.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     get number of tiles in the mosaic_file.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     ntiles = get_mosaic_ntiles( mosaic_file)
  !   </TEMPLATE>
  !   <IN NAME="mosaic_file" TYPE="character(len=*)">
  !     The file that contains mosaic information.
  !   </IN>
  function get_mosaic_ntiles(mosaic_file)
    character(len=*), intent(in) :: mosaic_file
    integer                      :: get_mosaic_ntiles 

    get_mosaic_ntiles = dimension_size(mosaic_file, "ntiles")

    return

  end function get_mosaic_ntiles
! </SUBROUTINE>

  !###############################################################################
  ! <SUBROUTINE NAME="get_mosaic_ncontacts">
  !   <OVERVIEW>
  !     get number of contacts in the mosaic_file.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     get number of contacts in the mosaic_file.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     ntiles = get_mosaic_ncontacts( mosaic_file)
  !   </TEMPLATE>
  !   <IN NAME="mosaic_file" TYPE="character(len=*)">
  !     The file that contains mosaic information.
  !   </IN>
  function get_mosaic_ncontacts( mosaic_file)
    character(len=*), intent(in) :: mosaic_file
    integer                      :: get_mosaic_ncontacts 

    character(len=len_trim(mosaic_file)+1) :: mfile    
    integer                                :: strlen
    integer                                :: read_mosaic_ncontacts

    if(field_exist(mosaic_file, "contacts") ) then
      get_mosaic_ncontacts = dimension_size(mosaic_file, "ncontact", no_domain=.TRUE.)
    else
      get_mosaic_ncontacts = 0
    endif

    return

  end function get_mosaic_ncontacts
! </SUBROUTINE>


  !###############################################################################
  ! <SUBROUTINE NAME="get_mosaic_grid_sizes">
  !   <OVERVIEW>
  !     get grid size of each tile from mosaic_file
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     get grid size of each tile from mosaic_file
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call get_mosaic_grid_sizes(mosaic_file, nx, ny)
  !   </TEMPLATE>
  !   <IN NAME="mosaic_file" TYPE="character(len=*)">
  !     The file that contains mosaic information.
  !   </IN>
  !   <INOUT NAME="nx" TYPE="integer, dimension(:)">
  !     List of grid size in x-direction of each tile.
  !   </INOUT>
  !   <INOUT NAME="ny" TYPE="integer, dimension(:)">
  !     List of grid size in y-direction of each tile.
  !   </INOUT>
  subroutine get_mosaic_grid_sizes( mosaic_file, nx, ny)
    character(len=*),         intent(in) :: mosaic_file
    integer, dimension(:), intent(inout) :: nx, ny

    character(len=MAX_FILE) :: gridfile
    integer                 :: ntiles, n    

    ntiles = get_mosaic_ntiles(mosaic_file)
    if(ntiles .NE. size(nx(:)) .OR. ntiles .NE. size(ny(:)) ) then
      call mpp_error(FATAL, "get_mosaic_grid_sizes: size of nx/ny does not equal to ntiles")
    endif
    do n = 1, ntiles
      call read_data(mosaic_file, 'gridfiles', gridfile, level=n)
      gridfile = grid_dir//trim(gridfile)
      nx(n) = dimension_size(gridfile, "nx")
      ny(n) = dimension_size(gridfile, "ny")
      if(mod(nx(n),x_refine) .NE. 0) call mpp_error(FATAL, "get_mosaic_grid_sizes: nx is not divided by x_refine");
      if(mod(ny(n),y_refine) .NE. 0) call mpp_error(FATAL, "get_mosaic_grid_sizes: ny is not divided by y_refine");
      nx(n) = nx(n)/x_refine;
      ny(n) = ny(n)/y_refine;
    enddo

    return

  end subroutine get_mosaic_grid_sizes
! </SUBROUTINE>

  !###############################################################################
  ! <SUBROUTINE NAME="get_mosaic_contact">
  !   <OVERVIEW>
  !     get contact information from mosaic_file
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     get contact information from mosaic_file
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call get_mosaic_contact(mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1,
  !                             istart2, iend2, jstart2, jend2)
  !   </TEMPLATE>
  !   <IN NAME="mosaic_file" TYPE="character(len=*)">
  !     The file that contains mosaic information.
  !   </IN>
  !   <INOUT NAME="tile1" TYPE="integer, dimension(:)">
  !     list tile number in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="tile1" TYPE="integer, dimension(:)">
  !     list tile number in tile 2 of each contact.
  !   </INOUT>
  !   <INOUT NAME="istart1" TYPE="integer, dimension(:)">
  !     list starting i-index in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="iend1" TYPE="integer, dimension(:)">
  !     list ending i-index in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="jstart1" TYPE="integer, dimension(:)">
  !     list starting j-index in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="jend1" TYPE="integer, dimension(:)">
  !     list ending j-index in tile 1 of each contact.
  !   </INOUT>
  !   <INOUT NAME="istart2" TYPE="integer, dimension(:)">
  !     list starting i-index in tile 2 of each contact.
  !   </INOUT>
  !   <INOUT NAME="iend2" TYPE="integer, dimension(:)">
  !     list ending i-index in tile 2 of each contact.
  !   </INOUT>
  !   <INOUT NAME="jstart2" TYPE="integer, dimension(:)">
  !     list starting j-index in tile 2 of each contact.
  !   </INOUT>
  !   <INOUT NAME="jend2" TYPE="integer, dimension(:)">
  !     list ending j-index in tile 2 of each contact.
  !   </INOUT>
  subroutine get_mosaic_contact( mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, &
                                   istart2, iend2, jstart2, jend2)
    character(len=*),         intent(in) :: mosaic_file
    integer, dimension(:), intent(inout) :: tile1, tile2
    integer, dimension(:), intent(inout) :: istart1, iend1, jstart1, jend1
    integer, dimension(:), intent(inout) :: istart2, iend2, jstart2, jend2
    character(len=MAX_NAME), allocatable :: gridtiles(:)
    character(len=MAX_NAME)              :: contacts
    character(len=MAX_NAME)              :: strlist(8)
    integer :: ntiles, n, m, ncontacts, nstr, ios
    integer :: i1_type, j1_type, i2_type, j2_type
    logical :: found

    ntiles = get_mosaic_ntiles(mosaic_file)
    allocate(gridtiles(ntiles))    
    do n = 1, ntiles
      call read_data(mosaic_file, 'gridtiles', gridtiles(n), level=n)
    enddo

    ncontacts = get_mosaic_ncontacts(mosaic_file)

    do n = 1, ncontacts
      call read_data(mosaic_file, "contacts", contacts, level=n)
      nstr = parse_string(contacts, ":", strlist)      
      if(nstr .NE. 4) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): number of elements in contact seperated by :/:: should be 4")
      found = .false.
      do m = 1, ntiles
        if(trim(gridtiles(m)) == trim(strlist(2)) ) then !found the tile name 
          found = .true.
          tile1(n) = m
          exit
        endif
      enddo
    
      if(.not.found) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact):the first tile name specified in contact is not found in tile list")

      found = .false.
      do m = 1, ntiles
        if(trim(gridtiles(m)) == trim(strlist(4)) ) then !found the tile name 
          found = .true.
          tile2(n) = m
          exit
        endif
      enddo
    
      if(.not.found) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact):the second tile name specified in contact is not found in tile list")

      call read_data(mosaic_file, "contact_index", contacts, level=n)
      nstr = parse_string(contacts, ":,", strlist)
      if(nstr .NE. 8) then
        if(mpp_pe()==mpp_root_pe()) then 
          print*, "nstr is ", nstr
          print*, "contacts is ", contacts
          do m = 1, nstr
            print*, "strlist is ", trim(strlist(m))
          enddo
        endif
        call mpp_error(FATAL, &
               "mosaic_mod(get_mosaic_contact): number of elements in contact_index seperated by :/, should be 8")
      endif
      read(strlist(1), *, iostat=ios) istart1(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading istart1")
      read(strlist(2), *, iostat=ios) iend1(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading iend1")
      read(strlist(3), *, iostat=ios) jstart1(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading jstart1")
      read(strlist(4), *, iostat=ios) jend1(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading jend1")
      read(strlist(5), *, iostat=ios) istart2(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading istart2")
      read(strlist(6), *, iostat=ios) iend2(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading iend2")
      read(strlist(7), *, iostat=ios) jstart2(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading jstart2")
      read(strlist(8), *, iostat=ios) jend2(n)
      if(ios .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): Error in reading jend2")

      i1_type = transfer_to_model_index(istart1(n), iend1(n), x_refine)
      j1_type = transfer_to_model_index(jstart1(n), jend1(n), y_refine)
      i2_type = transfer_to_model_index(istart2(n), iend2(n), x_refine)
      j2_type = transfer_to_model_index(jstart2(n), jend2(n), y_refine)

      if( i1_type == 0 .AND. j1_type == 0 ) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): istart1==iend1 and jstart1==jend1")
      if( i2_type == 0 .AND. j2_type == 0 ) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): istart2==iend2 and jstart2==jend2")
      if( i1_type + j1_type .NE. i2_type + j2_type ) call mpp_error(FATAL, &
         "mosaic_mod(get_mosaic_contact): It is not a line or overlap contact")

   enddo

      deallocate(gridtiles)

  end subroutine get_mosaic_contact
! </SUBROUTINE>


function transfer_to_model_index(istart, iend, refine_ratio)
   integer, intent(inout) :: istart, iend
   integer                :: refine_ratio
   integer                :: transfer_to_model_index
   integer                :: istart_in, iend_in

   istart_in = istart
   iend_in = iend

   if( istart_in == iend_in ) then
      transfer_to_model_index = 0
      istart = (istart_in + 1)/refine_ratio
      iend   = istart
   else 
      transfer_to_model_index = 1
      if( iend_in > istart_in ) then
        istart = istart_in + 1
        iend   = iend_in
      else 
        istart = istart_in
        iend   = iend_in + 1
      endif
      if( mod(istart, refine_ratio) .NE. 0 .OR. mod(iend,refine_ratio) .NE. 0) call mpp_error(FATAL, &
         "mosaic_mod(transfer_to_model_index): mismatch between refine_ratio and istart/iend")
      istart = istart/refine_ratio
      iend = iend/refine_ratio

   endif

   return

end function transfer_to_model_index

  !###############################################################################
  ! <SUBROUTINE NAME="calc_mosaic_grid_area">
  !   <OVERVIEW>
  !     calculate grid cell area.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     calculate the grid cell area. The purpose of this routine is to make 
  !     sure the consistency between model grid area and exchange grid area.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call calc_mosaic_grid_area(lon, lat, area)
  !   </TEMPLATE>
  !   <IN NAME="lon" TYPE="real, dimension(:,:)">
  !     geographical longitude of grid cell vertices.
  !   </IN>
  !   <IN NAME="lat" TYPE="real, dimension(:,:)">
  !     geographical latitude of grid cell vertices.
  !   </IN>
  !   <INOUT NAME="area" TYPE="real, dimension(:,:)">
  !     grid cell area.
  !   </INOUT>
  subroutine calc_mosaic_grid_area(lon, lat, area)
     real, dimension(:,:), intent(in)    :: lon
     real, dimension(:,:), intent(in)    :: lat
     real, dimension(:,:), intent(inout) :: area
     integer                             :: nlon, nlat

     nlon = size(area,1)
     nlat = size(area,2)
     ! make sure size of lon, lat and area are consitency
     if( size(lon,1) .NE. nlon+1 .OR. size(lat,1) .NE. nlon+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,1) and size(lat,1) should equal to size(area,1)+1")
     if( size(lon,2) .NE. nlat+1 .OR. size(lat,2) .NE. nlat+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,2) and size(lat,2) should equal to size(area,2)+1")

     call get_grid_area( nlon, nlat, lon, lat, area)

  end subroutine calc_mosaic_grid_area
  ! </SUBROUTINE>

  !###############################################################################
  ! <SUBROUTINE NAME="calc_mosaic_grid_great_circle_area">
  !   <OVERVIEW>
  !     calculate grid cell area using great cirlce algorithm
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !     calculate the grid cell area. The purpose of this routine is to make 
  !     sure the consistency between model grid area and exchange grid area.
  !   </DESCRIPTION>
  !   <TEMPLATE>
  !     call calc_mosaic_grid_great_circle_area(lon, lat, area)
  !   </TEMPLATE>
  !   <IN NAME="lon" TYPE="real, dimension(:,:)">
  !     geographical longitude of grid cell vertices.
  !   </IN>
  !   <IN NAME="lat" TYPE="real, dimension(:,:)">
  !     geographical latitude of grid cell vertices.
  !   </IN>
  !   <INOUT NAME="area" TYPE="real, dimension(:,:)">
  !     grid cell area.
  !   </INOUT>
  subroutine calc_mosaic_grid_great_circle_area(lon, lat, area)
     real, dimension(:,:), intent(in)    :: lon
     real, dimension(:,:), intent(in)    :: lat
     real, dimension(:,:), intent(inout) :: area
     integer                             :: nlon, nlat
     

     nlon = size(area,1)
     nlat = size(area,2)
     ! make sure size of lon, lat and area are consitency
     if( size(lon,1) .NE. nlon+1 .OR. size(lat,1) .NE. nlon+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,1) and size(lat,1) should equal to size(area,1)+1")
     if( size(lon,2) .NE. nlat+1 .OR. size(lat,2) .NE. nlat+1 ) &
        call mpp_error(FATAL, "mosaic_mod: size(lon,2) and size(lat,2) should equal to size(area,2)+1")

     call get_grid_great_circle_area( nlon, nlat, lon, lat, area)

  end subroutine calc_mosaic_grid_great_circle_area
  ! </SUBROUTINE>

  !#####################################################################
  ! This function check if a point (lon1,lat1) is inside a polygon (lon2(:), lat2(:))
  ! lon1, lat1, lon2, lat2 are in radians.
  function is_inside_polygon(lon1, lat1, lon2, lat2 )
     real, intent(in) :: lon1, lat1
     real, intent(in) :: lon2(:), lat2(:)
     logical          :: is_inside_polygon
     real, dimension(size(lon2(:))) :: x2, y2, z2
     integer                        :: npts, isinside
     integer                        :: inside_a_polygon

     npts = size(lon2(:))

     isinside = inside_a_polygon(lon1, lat1, npts, lon2, lat2)
     if(isinside == 1) then
        is_inside_polygon = .TRUE.
     else
        is_inside_polygon = .FALSE.
     endif

     return

  end function is_inside_polygon

  function parse_string(string, set, value)
  character(len=*),  intent(in) :: string
  character(len=*),  intent(in) :: set
  character(len=*), intent(out) :: value(:)
  integer                       :: parse_string
  integer :: nelem, length, first, last

     nelem = size(value(:))
     length = len_trim(string)

     first = 1; last = 0
     parse_string = 0

     do while(first .LE. length)
       parse_string = parse_string + 1
       if(parse_string>nelem) then
         call mpp_error(FATAL, "mosaic_mod(parse_string) : number of element is greater than size(value(:))")
       endif
       last = first - 1 + scan(string(first:length), set)
       if(last == first-1 ) then  ! not found, end of string
         value(parse_string) = string(first:length)
         exit
       else
         if(last <= first) then
           call mpp_error(FATAL, "mosaic_mod(parse_string) : last <= first")
         endif
         value(parse_string) = string(first:(last-1))
         first = last + 1
         ! scan to make sure the next is not the character in the set
         do while (first == last+1)
           last = first - 1 + scan(string(first:length), set)
           if(last == first) then
             first = first+1
           else
             exit
           endif
         end do
       endif
     enddo

     return

  end function parse_string


     
end module mosaic_mod


#ifdef TEST_MOSAIC
program test_mosaic

use mosaic_mod, only : get_mosaic_ntiles, get_mosaic_ncontacts
use mosaic_mod, only : get_mosaic_grid_sizes, get_mosaic_contact

implicit none

integer              :: ntiles, ncontacts, n
integer, allocatable :: tile1(:), tile2(:), nx(:), ny(:)
integer, allocatable :: istart1(:), iend1(:), jstart1(:), jend1(:)
integer, allocatable :: istart2(:), iend2(:), jstart2(:), jend2(:)
character(len=128)   :: mosaic_file = "INPUT/mosaic.nc"

ntiles = get_mosaic_ntiles(mosaic_file)
ncontacts = get_mosaic_ncontacts(mosaic_file)
allocate(nx(ntiles), ny(ntiles))
allocate(tile1(ncontacts), tile2(ncontacts) )
allocate(istart1(ncontacts), iend1(ncontacts), jstart1(ncontacts), jend1(ncontacts) )
allocate(istart2(ncontacts), iend2(ncontacts), jstart2(ncontacts), jend2(ncontacts) )

call get_mosaic_grid_sizes(mosaic_file, nx, ny )
call get_mosaic_contact(mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2)

! print out information

print '(a,i3,a,a)', "****** There is ", ntiles, " tiles in ", trim(mosaic_file)
do n = 1, ntiles
   print '(a,i3,a,i3,a,i3)', " tile = ", n, ", nx = ", nx(n), ", ny = ", ny(n)
end do

print '(a,i3,a,a)', "****** There is ", ncontacts, " contacts in ", trim(mosaic_file)
do n = 1, ncontacts
   print '(a,i3,a,i3,a,i3,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)', &
           "contact=", n, ": tile1=", tile1(n), " tile2=", tile2(n),   &
           " is1=", istart1(n), " ie1=", iend1(n),                   &
           " js1=", jstart1(n), " je1=", jend1(n),                   &
           " is2=", istart2(n), " ie2=", iend2(n),                   &
           " js2=", jstart2(n), " je2=", jend2(n)
end do

deallocate(tile1, tile2, nx, ny)
deallocate(istart1, iend1, jstart1, jend1)
deallocate(istart2, iend2, jstart2, jend2)


end program test_mosaic
#endif
