module gradient_mod
! <CONTACT EMAIL="Zhi.Liang@noaa.gov">
!   Zhi Liang
! </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!    <TT>gradient_mod</TT> implements some utility routines to calculate gradient.
! </OVERVIEW>

! <DESCRIPTION>
!    <TT>gradient_mod</TT> implements some utility routines to calculate gradient.
!    Currently only gradient on cubic grid is implemented. Also a public interface 
!    is provided to calculate grid information needed to calculate gradient.

use mpp_mod,       only : mpp_error, FATAL
use constants_mod, only : RADIUS

implicit none
private


public :: gradient_cubic
public :: calc_cubic_grid_info

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

contains


!#####################################################################
!  NOTe: pin has halo size = 1.
!  the size of pin    will be (nx+2,ny+2), T-cell center, with halo = 1
!  the size of dx     will be (nx, ny+1),  N-cell center
!  the size of dy     will be (nx+1, ny),  E-cell center
!  the size of area   will be (nx, ny),    T-cell center.
!  The size of edge_w will be (ny+1),      C-cell center
!  The size of edge_e will be (ny+1),      C-cell center
!  The size of edge_s will be (nx+1),      C-cell center
!  The size of edge_n will be (nx+1),      C-cell center
!  The size of en_n   will be (3,nx,ny+1), N-cell center
!  The size of en_e   will be (3,nx+1,ny), E-cell center
!  The size of vlon   will be (3,nx, ny)   T-cell center
!  The size of vlat   will be (3,nx, ny),  T-cell center

subroutine gradient_cubic(pin, dx, dy, area, edge_w, edge_e, edge_s, edge_n,    &
                          en_n, en_e, vlon, vlat, grad_x, grad_y, on_west_edge, &
                          on_east_edge, on_south_edge, on_north_edge)

  real,    dimension(:,:  ), intent(in ) :: pin, dx, dy, area
  real,    dimension(:    ), intent(in ) :: edge_w, edge_e, edge_s, edge_n
  real,    dimension(:,:,:), intent(in ) :: en_n, en_e
  real,    dimension(:,:,:), intent(in ) :: vlon, vlat
  real,    dimension(:,:  ), intent(out) :: grad_x, grad_y
  logical,                   intent(in ) :: on_west_edge, on_east_edge, on_south_edge, on_north_edge
  integer :: nx, ny


  nx = size(grad_x,1)
  ny = size(grad_x,2)

  if(size(pin,1) .NE. nx+2 .OR. size(pin,2) .NE. ny+2)call mpp_error(FATAL, "gradient_mod:size of pin should be (nx+2, ny+2)")
  if(size(dx,1) .NE. nx .OR. size(dx,2) .NE. ny+1 ) call mpp_error(FATAL, "gradient_mod: size of dx should be (nx,ny+1)")
  if(size(dy,1) .NE. nx+1 .OR. size(dy,2) .NE. ny ) call mpp_error(FATAL, "gradient_mod: size of dy should be (nx+1,ny)")
  if(size(area,1) .NE. nx .OR. size(area,2) .NE. ny ) call mpp_error(FATAL, "gradient_mod: size of area should be (nx,ny)")
  if(size(vlon,1) .NE. 3 .OR. size(vlon,2) .NE. nx .OR. size(vlon,3) .NE. ny) &
          call mpp_error(FATAL, "gradient_mod: size of vlon should be (3,nx,ny)")
  if(size(vlat,1) .NE. 3 .OR. size(vlat,2) .NE. nx .OR. size(vlat,3) .NE. ny) &
          call mpp_error(FATAL, "gradient_mod: size of vlat should be (3,nx,ny)")
  if(size(edge_w) .NE. ny+1) call mpp_error(FATAL, "gradient_mod: size of edge_w should be (ny+1)")
  if(size(edge_e) .NE. ny+1) call mpp_error(FATAL, "gradient_mod: size of edge_e should be (ny+1)")
  if(size(edge_s) .NE. nx+1) call mpp_error(FATAL, "gradient_mod: size of edge_s should be (nx+1)")
  if(size(edge_n) .NE. nx+1) call mpp_error(FATAL, "gradient_mod: size of edge_n should be (nx+1)")
  if(size(en_n,1) .NE. 3 .OR. size(en_n,2) .NE. nx .OR.  size(en_n,3) .NE. ny+1 ) &
       call mpp_error(FATAL, "gradient_mod:size of en_n should be (3, nx, ny+1)")
  if(size(en_e,1) .NE. 3 .OR. size(en_e,2) .NE. nx+1 .OR.  size(en_e,3) .NE. ny ) &
       call mpp_error(FATAL, "gradient_mod:size of en_e should be (3, nx+1, ny)")

  call grad_c2l(nx, ny, pin, dx, dy, area, edge_w, edge_e, edge_s, edge_n, en_n, en_e, vlon, vlat, &
                grad_x, grad_y, on_west_edge, on_east_edge, on_south_edge, on_north_edge)

  return

end subroutine gradient_cubic


subroutine calc_cubic_grid_info(xt, yt, xc, yc, dx, dy, area, edge_w, edge_e, edge_s, edge_n, &
                           en_n, en_e, vlon, vlat, on_west_edge, on_east_edge, on_south_edge, on_north_edge )
  real,    dimension(:,:  ), intent(in ) :: xt, yt, xc, yc
  real,    dimension(:,:  ), intent(out) :: dx, dy, area
  real,    dimension(:    ), intent(out) :: edge_w, edge_e, edge_s, edge_n
  real,    dimension(:,:,:), intent(out) :: en_n, en_e
  real,    dimension(:,:,:), intent(out) :: vlon, vlat
  logical,                   intent(in ) :: on_west_edge, on_east_edge, on_south_edge, on_north_edge
  integer :: nx, ny, nxp, nyp


  nx  = size(area,1)
  ny  = size(area,2)
  nxp = nx+1
  nyp = ny+1

  if(size(xt,1) .NE. nx+2 .OR. size(xt,2) .NE. ny+2 ) call mpp_error(FATAL, "gradient_mod: size of xt should be (nx+2,ny+2)")
  if(size(yt,1) .NE. nx+2 .OR. size(yt,2) .NE. ny+2 ) call mpp_error(FATAL, "gradient_mod: size of yt should be (nx+2,ny+2)")
  if(size(xc,1) .NE. nxp .OR. size(xc,2) .NE. nyp ) call mpp_error(FATAL, "gradient_mod: size of xc should be (nx+1,ny+1)")
  if(size(yc,1) .NE. nxp .OR. size(yc,2) .NE. nyp ) call mpp_error(FATAL, "gradient_mod: size of yc should be (nx+1,ny+1)")
  if(size(dx,1) .NE. nx .OR. size(dx,2) .NE. nyp ) call mpp_error(FATAL, "gradient_mod: size of dx should be (nx,ny+1)")
  if(size(dy,1) .NE. nxp .OR. size(dy,2) .NE. ny ) call mpp_error(FATAL, "gradient_mod: size of dy should be (nx+1,ny)")
  if(size(area,1) .NE. nx .OR. size(area,2) .NE. ny ) call mpp_error(FATAL, "gradient_mod: size of area should be (nx,ny)")
  if(size(vlon,1) .NE. 3 .OR. size(vlon,2) .NE. nx .OR. size(vlon,3) .NE. ny) &
          call mpp_error(FATAL, "gradient_mod: size of vlon should be (3,nx,ny)")
  if(size(vlat,1) .NE. 3 .OR. size(vlat,2) .NE. nx .OR. size(vlat,3) .NE. ny) &
          call mpp_error(FATAL, "gradient_mod: size of vlat should be (3,nx,ny)")
  if(size(edge_w) .NE. ny+1) call mpp_error(FATAL, "gradient_mod: size of edge_w should be (ny-1)")
  if(size(edge_e) .NE. ny+1) call mpp_error(FATAL, "gradient_mod: size of edge_e should be (ny-1)")
  if(size(edge_s) .NE. nx+1) call mpp_error(FATAL, "gradient_mod: size of edge_s should be (nx-1)")
  if(size(edge_n) .NE. nx+1) call mpp_error(FATAL, "gradient_mod: size of edge_n should be (nx-1)")
  if(size(en_n,1) .NE. 3 .OR. size(en_n,2) .NE. nx .OR.  size(en_n,3) .NE. nyp ) &
       call mpp_error(FATAL, "gradient_mod:size of en_n should be (3, nx, ny+1)")
  if(size(en_e,1) .NE. 3 .OR. size(en_e,2) .NE. nxp .OR.  size(en_e,3) .NE. ny ) &
       call mpp_error(FATAL, "gradient_mod:size of en_e should be (3, nx+1, ny)")


  call calc_c2l_grid_info(nx, ny, xt, yt, xc, yc, dx, dy, area, edge_w, edge_e, edge_s, edge_n, &
                          en_n, en_e, vlon, vlat, on_west_edge, on_east_edge, on_south_edge, on_north_edge )


  return

end subroutine calc_cubic_grid_info

end module gradient_mod
