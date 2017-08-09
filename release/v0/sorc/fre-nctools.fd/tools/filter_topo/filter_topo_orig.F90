program filter_topo


  implicit none

#include <netcdf.inc>

#ifdef NO_QUAD_PRECISION
  ! 64-bit precision (kind=8)
  integer, parameter:: f_p = selected_real_kind(15)
#else
  ! Higher precision (kind=16) for grid geometrical factors:
  integer, parameter:: f_p = selected_real_kind(20)
#endif

  integer, parameter :: XDir=1
  integer, parameter :: YDir=2
  real, parameter :: pi = 3.14159265358979323846d0
  real, parameter :: radius = 6371.d03
  real, parameter ::  big_number=1.d8
  real, parameter :: tiny_number=1.d-8


  real:: cd4 = 0.16        ! Dimensionless coeff for del-4 difussion (with FCT)
  real:: peak_fac = 1.05   ! overshoot factor for the mountain peak
  real:: max_slope = 0.15       ! max allowable terrain slope: 1 --> 45 deg
                                ! 0.15 for C768 or lower; 0.25 C1536; 0.3 for C3072
  integer :: n_del2_weak = 12

  logical :: zs_filter = .true. 
  logical :: zero_ocean = .true.          ! if true, no diffusive flux into water/ocean area 
  real    :: stretch_fac = 1.0
  logical :: nested = .false.
  integer :: grid_type = 0 ! gnomoic_ed
  character(len=128) :: topo_file = "orog"
  character(len=128) :: topo_field = "orog_filt"
  character(len=128) :: mask_field = "slmsk"
  character(len=128) :: grid_file = "atmos_mosaic.nc"
  namelist /filter_topo_nml/ topo_file, topo_field, mask_field, grid_file, zero_ocean, &
       zs_filter, cd4, n_del2_weak, peak_fac, max_slope, stretch_fac, nested, grid_type

  integer :: stdunit = 6 
  integer :: ntiles = 0

  real da_min

  real, allocatable :: oro(:,:,:), mask(:,:,:)
  real, allocatable :: dx(:,:,:), dy(:,:,:)
  real, allocatable :: dxa(:,:,:), dya(:,:,:)
  real, allocatable :: dxc(:,:,:), dyc(:,:,:)
  real, allocatable :: area(:,:,:)
  real, allocatable :: sin_sg(:,:,:,:)

  integer           :: is,ie,js,je,isd,ied,jsd,jed
  integer,parameter :: ng = 3
  integer           :: nx, ny, npx, npy

  !--- read namelist
  call read_namelist()

  !--- read the target grid.
  call read_grid_file()

  !--- read the topography data
  call read_topo_file

  !--- filter the data
  call FV3_zs_filter(is,ie,js,je,isd,ied,jsd,jed,npx,npy,npx,ntiles,grid_type, &
                     stretch_fac, nested, area, dxa, dya, dx, dy, dxc, dyc, sin_sg, oro )

  !--- write out the data
  call write_topo_file(is,ie,js,je,ntiles,oro(is:ie,js:je,:) )

contains

  !#####################################################################
  real function great_circle_dist( q1, q2, radius )
    real, intent(IN)           :: q1(2), q2(2)
    real, intent(IN), optional :: radius

    real (f_p):: p1(2), p2(2)
    real (f_p):: beta
    integer n

    do n=1,2
       p1(n) = q1(n)
       p2(n) = q2(n)
    enddo

    beta = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
         sin((p1(1)-p2(1))/2.)**2 ) ) * 2.

    if ( present(radius) ) then
       great_circle_dist = radius * beta
    else
       great_circle_dist = beta   ! Returns the angle
    endif

  end function great_circle_dist

  !####################################################################
  real function spherical_angle(p1, p2, p3)

    !           p3
    !         /
    !        /
    !       p1 ---> angle
    !         \
    !          \
    !           p2

    real p1(3), p2(3), p3(3)

    real (f_p):: e1(3), e2(3), e3(3)
    real (f_p):: px, py, pz
    real (f_p):: qx, qy, qz
    real (f_p):: angle, ddd
    integer n

    do n=1,3
       e1(n) = p1(n)
       e2(n) = p2(n)
       e3(n) = p3(n)
    enddo

    !-------------------------------------------------------------------
    ! Page 41, Silverman's book on Vector Algebra; spherical trigonmetry
    !-------------------------------------------------------------------
    ! Vector P:
    px = e1(2)*e2(3) - e1(3)*e2(2)
    py = e1(3)*e2(1) - e1(1)*e2(3)
    pz = e1(1)*e2(2) - e1(2)*e2(1)
    ! Vector Q:
    qx = e1(2)*e3(3) - e1(3)*e3(2)
    qy = e1(3)*e3(1) - e1(1)*e3(3)
    qz = e1(1)*e3(2) - e1(2)*e3(1)

    ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz)

    if ( ddd <= 0.0d0 ) then
       angle = 0.d0
    else
       ddd = (px*qx+py*qy+pz*qz) / sqrt(ddd)
       if ( abs(ddd)>1.d0) then
          angle = 2.d0*atan(1.0)    ! 0.5*pi
          !FIX (lmh) to correctly handle co-linear points (angle near pi or 0)
          if (ddd < 0.d0) then
             angle = 4.d0*atan(1.0d0) !should be pi
          else
             angle = 0.d0
          end if
       else
          angle = acos( ddd )
       endif
    endif

    spherical_angle = angle

  end function spherical_angle


  !####################################################################
  real function get_area(p1, p4, p2, p3, radius)
    !-----------------------------------------------
    real, intent(in), dimension(2):: p1, p2, p3, p4
    real, intent(in), optional:: radius
    !-----------------------------------------------
    real e1(3), e2(3), e3(3)
    real ang1, ang2, ang3, ang4

    ! S-W: 1
    call latlon2xyz(p1, e1)   ! p1
    call latlon2xyz(p2, e2)   ! p2
    call latlon2xyz(p4, e3)   ! p4
    ang1 = spherical_angle(e1, e2, e3)
    !----
    ! S-E: 2
    !----
    call latlon2xyz(p2, e1)
    call latlon2xyz(p3, e2)
    call latlon2xyz(p1, e3)
    ang2 = spherical_angle(e1, e2, e3)
    !----
    ! N-E: 3
    !----
    call latlon2xyz(p3, e1)
    call latlon2xyz(p4, e2)
    call latlon2xyz(p2, e3)
    ang3 = spherical_angle(e1, e2, e3)
    !----
    ! N-W: 4
    !----
    call latlon2xyz(p4, e1)
    call latlon2xyz(p3, e2)
    call latlon2xyz(p1, e3)
    ang4 = spherical_angle(e1, e2, e3)

    if ( present(radius) ) then
       get_area = (ang1 + ang2 + ang3 + ang4 - 2.*pi) * radius**2
    else
       get_area = ang1 + ang2 + ang3 + ang4 - 2.*pi
    endif

  end function get_area

  subroutine fill_AGRID_scalar_corners(q, ng, npx, npy, isd, jsd, fill)
    integer, intent(in)  :: ng, npx, npy, isd, jsd
    integer, intent(in)  :: fill
    real, DIMENSION(isd:,jsd:,:), intent(INOUT):: q

    integer :: i, j

    select case (FILL)
    case (XDir)
       do j=1,ng
          do i=1,ng
             q(1-i    ,1-j    ,:) = q(1-j    ,i        ,:)  !SW Corner 
             q(1-i    ,npy-1+j,:) = q(1-j    ,npy-1-i+1,:)  !NW Corner
             q(npx-1+i,1-j    ,:) = q(npx-1+j,i        ,:)  !SE Corner
             q(npx-1+i,npy-1+j,:) = q(npx-1+j,npy-1-i+1,:)  !NE Corner
          enddo
       enddo
    case (YDir)
       do j=1,ng
          do i=1,ng
             q(1-j    ,1-i    ,:) = q(i        ,1-j    ,:)  !SW Corner 
             q(1-j    ,npy-1+i,:) = q(i        ,npy-1+j,:)  !NW Corner
             q(npx-1+j,1-i    ,:) = q(npx-1-i+1,1-j    ,:)  !SE Corner
             q(npx-1+j,npy-1+i,:) = q(npx-1-i+1,npy-1+j,:)  !NE Corner
          enddo
       enddo
    case default
       do j=1,ng
          do i=1,ng
             q(1-j    ,1-i    ,:) = q(i        ,1-j    ,:)  !SW Corner 
             q(1-j    ,npy-1+i,:) = q(i        ,npy-1+j,:)  !NW Corner
             q(npx-1+j,1-i    ,:) = q(npx-1-i+1,1-j    ,:)  !SE Corner
             q(npx-1+j,npy-1+i,:) = q(npx-1-i+1,npy-1+j,:)  !NE Corner
          enddo
       enddo
    end select


  end subroutine fill_AGRID_scalar_corners


  !####################################################################
  subroutine fill_BGRID_scalar_corners(q, ng, npx, npy, isd, jsd, fill)
    integer, intent(in)  :: ng, npx, npy, isd, jsd
    integer, intent(in)  :: fill
    real, DIMENSION(isd:,jsd:,:), intent(INOUT):: q

    integer :: i, j

    select case (fill)
    case (XDir)
       do j=1,ng
          do i=1,ng
             q(1-i  ,1-j  ,:) = q(1-j  ,i+1    ,:)  !SW Corner 
             q(1-i  ,npy+j,:) = q(1-j  ,npy-i  ,:)  !NW Corner
             q(npx+i,1-j  ,:) = q(npx+j,i+1    ,:)  !SE Corner
             q(npx+i,npy+j,:) = q(npx+j,npy-i  ,:)  !NE Corner
          enddo
       enddo
    case (YDir)
       do j=1,ng
          do i=1,ng
             q(1-j  ,1-i  ,:) = q(i+1  ,1-j    ,:)  !SW Corner 
             q(1-j  ,npy+i,:) = q(i+1  ,npy+j  ,:)  !NW Corner
             q(npx+j,1-i  ,:) = q(npx-i,1-j    ,:)  !SE Corner
             q(npx+j,npy+i,:) = q(npx-i,npy+j  ,:)  !NE Corner
          enddo
       enddo
    case default
       do j=1,ng
          do i=1,ng
             q(1-i  ,1-j  ,:) = q(1-j  ,i+1    ,:)  !SW Corner 
             q(1-i  ,npy+j,:) = q(1-j  ,npy-i  ,:)  !NW Corner
             q(npx+i,1-j  ,:) = q(npx+j,i+1    ,:)  !SE Corner
             q(npx+i,npy+j,:) = q(npx+j,npy-i  ,:)  !NE Corner
          enddo
       enddo
    end select



  end subroutine fill_BGRID_scalar_corners

  !####################################################################
  subroutine fill_AGRID_xy_corners(x, y, ng, npx, npy, isd, jsd)
    integer, intent(in)                      :: ng, npx, npy, isd, jsd
    real, DIMENSION(isd:,jsd:,:), intent(INOUT):: x
    real, DIMENSION(isd:,jsd:,:), intent(INOUT):: y
    integer :: i,j

    do j=1,ng
       do i=1,ng
          x(1-i    ,1-j    ,:) = y(1-j    ,i        ,:)  !SW Corner
          x(1-i    ,npy-1+j,:) = y(1-j    ,npy-1-i+1,:)  !NW Corner
          x(npx-1+i,1-j    ,:) = y(npx-1+j,i        ,:)  !SE Corner
          x(npx-1+i,npy-1+j,:) = y(npx-1+j,npy-1-i+1,:)  !NE Corner

          y(1-j    ,1-i    ,:) = x(i        ,1-j    ,:)  !SW Corner
          y(1-j    ,npy-1+i,:) = x(i        ,npy-1+j,:)  !NW Corner
          y(npx-1+j,1-i    ,:) = x(npx-1-i+1,1-j    ,:)  !SE Corner
          y(npx-1+j,npy-1+i,:) = x(npx-1-i+1,npy-1+j,:)  !NE Corner
       enddo
    enddo

  end subroutine fill_AGRID_xy_corners



  !####################################################################
  subroutine fill_DGRID_xy_corners(x, y, ng, npx, npy, isd, jsd)
    integer, intent(in)                      :: ng, npx, npy, isd, jsd
    real, DIMENSION(isd:,jsd:,:), intent(INOUT):: x
    real, DIMENSION(isd:,jsd:,:), intent(INOUT):: y
    integer :: i,j

    do j=1,ng
       do i=1,ng
          x(1-i    ,1-j    , :) = y(1-j  ,i    , :)  !SW Corner 
          x(1-i    ,npy+j  , :) = y(1-j  ,npy-i, :)  !NW Corner
          x(npx-1+i,1-j    , :) = y(npx+j,i    , :)  !SE Corner
          x(npx-1+i,npy+j  , :) = y(npx+j,npy-i, :)  !NE Corner
          y(1-i    ,1-j    , :) = x(j    ,1-i  , :)  !SW Corner 
          y(1-i    ,npy-1+j, :) = x(j    ,npy+i, :)  !NW Corner
          y(npx+i  ,1-j    , :) = x(npx-j,1-i  , :)  !SE Corner
          y(npx+i  ,npy-1+j, :) = x(npx-j,npy+i, :)  !NE Corner
       enddo
    enddo

  end subroutine fill_DGRID_xy_corners

  !###############################################################
  subroutine mid_pt_sphere(p1, p2, pm)
    real, intent(IN)  :: p1(2), p2(2)
    real, intent(OUT) :: pm(2)
    !------------------------------------------
    real :: e1(3), e2(3), e3(3)

    call latlon2xyz(p1, e1)
    call latlon2xyz(p2, e2)
    call mid_pt3_cart(e1, e2, e3)
    call cart_to_latlon(1, e3, pm(1), pm(2))

  end subroutine mid_pt_sphere


  !#####################################################################
  subroutine mid_pt3_cart(p1, p2, e)
    real, intent(IN)  :: p1(3), p2(3)
    real, intent(OUT) :: e(3)
    !
    real (f_p):: q1(3), q2(3)
    real (f_p):: dd, e1, e2, e3
    integer k

    do k=1,3
       q1(k) = p1(k)
       q2(k) = p2(k)
    enddo

    e1 = q1(1) + q2(1)
    e2 = q1(2) + q2(2)
    e3 = q1(3) + q2(3)

    dd = sqrt( e1**2 + e2**2 + e3**2 )
    e1 = e1 / dd
    e2 = e2 / dd
    e3 = e3 / dd

    e(1) = e1
    e(2) = e2
    e(3) = e3

  end subroutine mid_pt3_cart

  subroutine latlon2xyz(p, e)
    !
    ! Routine to map (lon, lat) to (x,y,z)
    !
    real, intent(in) :: p(2)
    real, intent(out):: e(3)

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
    real, intent(inout):: q(3,np)
    real, intent(inout):: xs(np), ys(np)
    ! local
    real, parameter:: esl=1.d-10
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

  !#####################################################################
  real function cos_angle(p1, p2, p3)
    ! As spherical_angle, but returns the cos(angle)
    !       p3
    !       ^  
    !       |  
    !       | 
    !       p1 ---> p2
    !
    real, intent(in):: p1(3), p2(3), p3(3)

    real (f_p):: e1(3), e2(3), e3(3)
    real (f_p):: px, py, pz
    real (f_p):: qx, qy, qz
    real (f_p):: angle, ddd
    integer n

    do n=1,3
       e1(n) = p1(n)
       e2(n) = p2(n)
       e3(n) = p3(n)
    enddo

    !-------------------------------------------------------------------
    ! Page 41, Silverman's book on Vector Algebra; spherical trigonmetry
    !-------------------------------------------------------------------
    ! Vector P:= e1 X e2
    px = e1(2)*e2(3) - e1(3)*e2(2)
    py = e1(3)*e2(1) - e1(1)*e2(3)
    pz = e1(1)*e2(2) - e1(2)*e2(1)

    ! Vector Q: e1 X e3
    qx = e1(2)*e3(3) - e1(3)*e3(2)
    qy = e1(3)*e3(1) - e1(1)*e3(3)
    qz = e1(1)*e3(2) - e1(2)*e3(1)

    ! ddd = sqrt[ (P*P) (Q*Q) ]
    ddd = sqrt( (px**2+py**2+pz**2)*(qx**2+qy**2+qz**2) )
    if ( ddd > 0.d0 ) then
       angle = (px*qx+py*qy+pz*qz) / ddd
    else
       angle = 1.d0
    endif
    cos_angle = angle

  end function cos_angle


  !#####################################################################
  subroutine cell_center2(q1, q2, q3, q4, e2)
    real, intent(in ) :: q1(2), q2(2), q3(2), q4(2)
    real, intent(out) :: e2(2)
    ! Local
    real p1(3), p2(3), p3(3), p4(3)
    real ec(3)
    real dd
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


  !#####################################################################
  subroutine read_grid_file()

    integer :: fsize=65536
    integer :: status, ncid, id_dim, id_var, ncid2, t
    integer :: ni, nj, i, j, tw, te, ip
    real    :: g1(2), g2(2), g3(2), g4(2), g5(2)
    real    :: p1(3), p3(3)
    real    :: p_lL(2), p_uL(2), p_lR(2), p_uR(2)
    character(len=256) :: tile_file
    real, allocatable, dimension(:,:)   :: tmpvar
    real, allocatable, dimension(:,:,:) :: geolon_c, geolat_c
    real, allocatable, dimension(:,:,:) :: geolon_t, geolat_t, cos_sg, grid3
    integer :: start(4), nread(4)

    print*, "Read the grid from file "//trim(grid_file)

    status=NF__OPEN(trim(grid_file),NF_NOWRITE,fsize,ncid)
    call handle_err(status, 'Open file '//trim(grid_file) )

    status=nf_inq_dimid(ncid, 'ntiles', id_dim)
    call handle_err(status, 'inquire dimension ntiles from file '//trim(grid_file) )
    status=nf_inq_dimlen(ncid,id_dim,ntiles)
    call handle_err(status, 'inquire dimension ntiles length from file '//trim(grid_file) )

    !--- currently only support cubic sphere grid.
    if( ntiles .NE. 6 .and. ntiles .NE. 7) call handle_err(-1, "ntiles should be 6 or 7 for file "//trim(grid_file) )
    if( ntiles == 7 ) print*, " This grid is a nested grid "

    !--- loop through ntiles and make sure the grid size match between all the tiles.

    start(:) = 1
    nread(:) = 1

    do t = 1, ntiles
       start(2) = t; nread(1) = 255
       status =  nf_inq_varid(ncid, 'gridfiles', id_var)
       call handle_err(status, 'inquire varid of gridfiles from file '//trim(grid_file) )
       status = nf_get_vara_text(ncid, id_var, start, nread, tile_file )      
       call handle_err(status, 'get value of gridfiles from file '//trim(grid_file) )

       status=NF__OPEN(trim(tile_file),NF_NOWRITE,fsize,ncid2)
       call handle_err(status, 'Open file '//trim(tile_file) )

       status=nf_inq_dimid(ncid2, 'nx', id_dim)
       call handle_err(status, 'inquire dimension nx from file '//trim(grid_file) )
       status=nf_inq_dimlen(ncid2,id_dim,ni)
       call handle_err(status, 'inquire dimension nx length from file '//trim(grid_file) )
       status=nf_inq_dimid(ncid2, 'ny', id_dim)
       call handle_err(status, 'inquire dimension ny from file '//trim(grid_file) )
       status=nf_inq_dimlen(ncid2,id_dim,nj)
       call handle_err(status, 'inquire dimension ny length '//'from file '//trim(grid_file) )
       if( t == 1 ) then
          ! ni and nj must be even
          if(mod(ni,2) .NE. 0 .or. mod(nj,2) .NE. 0) &
               call handle_err(-1, "read_grid_file: ni and nj must be even")

          nx = ni/2
          ny = nj/2
          npx = nx + 1 
          npy = ny + 1
          is = 1  ; ie = nx
          js = 1  ; je = ny
          isd=is-ng; ied=ie+ng
          jsd=js-ng; jed=je+ng

          allocate(tmpvar(ni+1,nj+1))
          allocate(geolon_c(isd:ied+1,jsd:jed+1,6))
          allocate(geolat_c(isd:ied+1,jsd:jed+1,6))
       else if ( t == 7 ) then ! nested grid
          if(mod(ni,2) .NE. 0 .or. mod(nj,2) .NE. 0) &
               call handle_err(-1, "read_grid_file: ni and nj must be even")

          nx_nest = ni/2
          ny_nest = nj/2
          npx_nest = nx_nest + 1 
          npy_nest = ny_nest + 1
          is_nest = 1  ; ie_nest = nx
          js_nest = 1  ; je_nest = ny
          isd_nest=is_nest-ng; ied_nest=ie_nest+ng
          jsd_nest=js_nest-ng; jed_nest=je_nest+ng
          deallocate(tmpvar)
          allocate(tmpvar(ni+1,nj+1))
          allocate(geolon_c_nest(isd:ied+1,jsd:jed+1))
          allocate(geolat_c_nest(isd:ied+1,jsd:jed+1))
       else
          !-- make sure ni and nj match between tiles
          if(ni .ne. nx*2 .OR. nj .ne. ny*2) &
               call handle_err(-1, "mismatch of grid size between tiles")
       endif

       status=nf_inq_varid(ncid2, 'x', id_var)
       call handle_err(status, 'inquire varid of x from file '//trim(grid_file) )
       status=nf_get_var_double(ncid2, id_var, tmpvar)
       call handle_err(status, 'inquire data of x from file '//trim(grid_file) )
       if(t==7) then
         geolon_c_nest(1:npx,1:npy) = tmpvar(1:ni+1:2,1:nj+1:2)*PI/180.
       else
         geolon_c(1:npx,1:npy,t) = tmpvar(1:ni+1:2,1:nj+1:2)*PI/180.
       endif

       status=nf_inq_varid(ncid2, 'y', id_var)
       call handle_err(status, 'inquire varid of y from file '//trim(grid_file) )
       status=nf_get_var_double(ncid2, id_var, tmpvar)
       call handle_err(status, 'inquire data of y from file '//trim(grid_file) )
       if(t==7) then
         geolat_c_nest(1:npx,1:npy) = tmpvar(1:ni+1:2,1:nj+1:2)*PI/180.
       else
         geolat_c(1:npx,1:npy,t) = tmpvar(1:ni+1:2,1:nj+1:2)*PI/180.
       endif
       status = nf_close(ncid2)
       call handle_err(status, "close file "//trim(tile_file))      
    enddo

    deallocate(tmpvar)

    status = nf_close(ncid)
    call handle_err(status, "close file "//trim(grid_file))

    is = 1  ; ie = nx
    js = 1  ; je = ny
    isd=is-ng; ied=ie+ng
    jsd=js-ng; jed=je+ng

    call fill_cubic_grid_halo(geolon_c, geolon_c, ng, 1, 1, 1, 1)
    call fill_cubic_grid_halo(geolat_c, geolat_c, ng, 1, 1, 1, 1)    
    if(.not. nested) call fill_bgrid_scalar_corners(geolon_c, ng, npx, npy, isd, jsd, XDir)
    if(.not. nested) call fill_bgrid_scalar_corners(geolat_c, ng, npx, npy, isd, jsd, YDir)

    !--- compute grid cell center
    allocate(geolon_t(isd:ied,jsd:jed,ntiles), geolat_t(isd:ied,jsd:jed,ntiles))

    geolon_t(:,:,:) = -1.e25
    geolat_t(:,:,:) = -1.e25

    do t = 1, ntiles
       do j=js,je ; do i=is,ie
          g1(1) = geolon_c(i,j,t);     g1(2) = geolat_c(i,j,t)
          g2(1) = geolon_c(i+1,j,t);   g2(2) = geolat_c(i+1,j,t)
          g3(1) = geolon_c(i,j+1,t);   g3(2) = geolat_c(i,j+1,t)
          g4(1) = geolon_c(i+1,j+1,t); g4(2) = geolat_c(i+1,j+1,t)
          call cell_center2(g1, g2, g3, g4, g5 )
          geolon_t(i,j,t) = g5(1)
          geolat_t(i,j,t) = g5(2)
       enddo ; enddo
    enddo

    
    call fill_cubic_grid_halo(geolon_t, geolon_t, ng, 0, 0, 1, 1)
    call fill_cubic_grid_halo(geolat_t, geolat_t, ng, 0, 0, 1, 1)

    if (.not. nested) call fill_AGRID_scalar_corners(geolon_t, ng, npx, npy, isd, jsd, XDir)
    if (.not. nested) call fill_AGRID_scalar_corners(geolat_t, ng, npx, npy, isd, jsd, YDir)

    !--- compute dx, dy
    allocate(dx(isd:ied,jsd:jed+1,ntiles))
    allocate(dy(isd:ied+1,jsd:jed,ntiles))
    do t = 1, ntiles
       do j = js, je+1 ; do i = is, ie
          g1(1) = geolon_c(i  ,j,t)
          g1(2) = geolat_c(i  ,j,t)
          g2(1) = geolon_c(i+1,j,t)
          g2(2) = geolat_c(i+1,j,t)
          dx(i,j,t) = great_circle_dist( g2, g1, radius )
       enddo ; enddo
    enddo
    if( stretch_fac .NE. 1 ) then
       do t = 1, ntiles
          do j = js, je
             do i = is, ie+1
                g1(1) = geolon_c(i,j,  t)
                g1(2) = geolat_c(i,j,  t)
                g2(1) = geolon_c(i,j+1,t)
                g2(2) = geolat_c(i,j+1,t)
                dy(i,j,t) = great_circle_dist( g2, g1, radius )
             enddo
          enddo
       enddo
    else
       do t = 1, ntiles
          do j = js, je
             do i = is, ie+1
                dy(i,j,t) = dx(j,i,t)
             enddo
          enddo
       enddo
    endif

    !--- make sure it is consitent between tiles. The following maybe not necessary.
    do t = 1, ntiles
       if(mod(t,2) ==0) then ! tile 2 4 6
          tw = t - 1
          te = t + 2
          if(te > ntiles) te = te - ntiles
          dy(is, js:je,t) = dy(ie+1,js:je,tw)        ! west boundary
          dy(ie+1, js:je, t) = dx(ie:is:-1,js, te)  ! east boundary
       else
          tw = t - 2
          if( tw <= 0) tw = tw + ntiles
          te = t + 1  
          dy(is, js:je, t) = dx(ie:is:-1, je+1, tw)  ! west boundary
          dy(ie+1, js:je,t) = dy(1,js:je,te)        ! east boundary
       endif
    enddo

    call fill_cubic_grid_halo(dx, dy, ng, 0, 1, 1, 1)
    call fill_cubic_grid_halo(dy, dx, ng, 1, 0, 1, 1)

    if (.not. nested) call fill_dgrid_xy_corners(dx, dy, ng, npx, npy, isd, jsd)

    !--- compute dxa and dya -----
    allocate(dxa(isd:ied,jsd:jed,ntiles))
    allocate(dya(isd:ied,jsd:jed,ntiles))
    do t = 1, ntiles
       do j=js,je ; do i=is,ie
          g1(1) = geolon_c(i,j,t); g1(2) = geolat_c(i,j,t)
          g2(1) = geolon_c(i,j+1,t); g2(2) = geolat_c(i,j+1,t)
          call mid_pt_sphere(g1, g2, g3)
          g1(1) = geolon_c(i+1,j,t); g1(2) = geolat_c(i+1,j,t)
          g2(1) = geolon_c(i+1,j+1,t); g2(2) = geolat_c(i+1,j+1,t)
          call mid_pt_sphere(g1, g2, g4)
          dxa(i,j,t) = great_circle_dist( g4, g3, radius )
          g1(1) = geolon_c(i,j,t); g1(2) = geolat_c(i,j,t)
          g2(1) = geolon_c(i+1,j,t); g2(2) = geolat_c(i+1,j,t)
          call mid_pt_sphere(g1, g2, g3)
          g1(1) = geolon_c(i,j+1,t); g1(2) = geolat_c(i,j+1,t)
          g2(1) = geolon_c(i+1,j+1,t); g2(2) = geolat_c(i+1,j+1,t)
          call mid_pt_sphere(g1, g2, g4)
          dya(i,j,t) = great_circle_dist( g4, g3, radius )
       enddo; enddo
    enddo

    call fill_cubic_grid_halo(dxa, dya, ng, 0, 0, 1, 1)
    call fill_cubic_grid_halo(dya, dxa, ng, 0, 0, 1, 1)
    
    if (.not. nested) call fill_AGRID_xy_corners(dxa, dya, ng, npx, npy, isd, jsd)

    !--- compute dxc and dyc
    allocate(dxc(isd:ied+1,jsd:jed,ntiles))
    allocate(dyc(isd:ied,jsd:jed+1,ntiles))
    do t = 1, ntiles
       do j=jsd,jed
          do i=isd+1,ied
             g1(1) = geolon_c(i,j,t); g1(2) = geolat_c(i,j,t)
             g2(1) = geolon_c(i-1,j,t); g2(2) = geolat_c(i-1,j,t)
             dxc(i,j,t) = great_circle_dist(g1, g2, radius)
          enddo
          dxc(isd,j,t)   = dxc(isd+1,j,t)
          dxc(ied+1,j,t) = dxc(ied,j,t)
       enddo

       do j=jsd+1,jed
          do i=isd,ied
             g1(1) = geolon_c(i,j,t); g1(2) = geolat_c(i,j,t)
             g2(1) = geolon_c(i,j-1,t); g2(2) = geolat_c(i,j-1,t)
             dyc(i,j,t) = great_circle_dist(g1, g2, radius)
          enddo
       enddo
       do i=isd,ied
          dyc(i,jsd,t)   = dyc(i,jsd+1,t)
          dyc(i,jed+1,t) = dyc(i,jed,t)
       end do
    enddo

    !--- compute area
    allocate(area(isd:ied,jsd:jed,ntiles))
    do t = 1, ntiles
       do j=js,je
          do i=is,ie
             p_lL(1) = geolon_c(i  ,j  ,t) ; p_lL(2) = geolat_c(i  ,j  ,t)
             p_uL(1) = geolon_c(i  ,j+1,t) ; p_uL(2) = geolat_c(i  ,j+1,t)
             p_lR(1) = geolon_c(i+1,j  ,t) ; p_lR(2) = geolat_c(i+1,j  ,t)
             p_uR(1) = geolon_c(i+1,j+1,t) ; p_uR(2) = geolat_c(i+1,j+1,t)

             ! Spherical Excess Formula
             area(i,j,t) = get_area(p_lL, p_uL, p_lR, p_uR, radius)
          enddo
       enddo
    enddo

    call fill_cubic_grid_halo(area, area, ng, 0, 0, 1, 1)

    da_min = minval(area(is:ie,js:je,:))

    !--- compute sin_sg
    allocate(sin_sg(4,isd:ied,jsd:jed,ntiles))
    allocate(cos_sg(4,isd:ied,jsd:jed))
    allocate(grid3(3, npx, npy))
    cos_sg(:,:,:) =  big_number
    sin_sg(:,:,:,:) = tiny_number

    !     9---4---8
    !     |       |
    !     1   5   3
    !     |       |
    !     6---2---7
    do t = 1, ntiles
       do j=js,je+1
          do i = is,ie+1
             g1(1) = geolon_c(i,j,t)
             g1(2) = geolat_c(i,j,t)
             call latlon2xyz(g1, grid3(:,i,j))
          enddo
       enddo
       do j=js,je
          do i=is,ie
             g1(1) = geolon_t(i,j,t); g1(2) = geolat_t(i,j,t)
             call latlon2xyz(g1, p3)   ! righ-hand system consistent with grid3
             call mid_pt3_cart(grid3(1,i,j), grid3(1,i,j+1), p1)
             cos_sg(1,i,j) = cos_angle( p1, p3, grid3(1,i,j+1) )
             call mid_pt3_cart(grid3(1,i,j), grid3(1,i+1,j), p1)
             cos_sg(2,i,j) = cos_angle( p1, grid3(1,i+1,j), p3 )
             call mid_pt3_cart(grid3(1,i+1,j), grid3(1,i+1,j+1), p1)
             cos_sg(3,i,j) = cos_angle( p1, p3, grid3(1,i+1,j) )
             call mid_pt3_cart(grid3(1,i,j+1), grid3(1,i+1,j+1), p1)
             cos_sg(4,i,j) = cos_angle( p1, grid3(1,i,j+1), p3 )
          enddo
       enddo

       do ip=1,4
          do j=js,je
             do i=is,ie
                sin_sg(ip,i,j,t) = min(1.0, sqrt( max(0., 1.-cos_sg(ip,i,j)**2) ) )
             enddo
          enddo
       enddo
    enddo

    do ip=1,4
       call fill_cubic_grid_halo(sin_sg(ip,:,:,:), sin_sg(ip,:,:,:), ng, 0, 0, 1, 1)
    enddo

    deallocate(cos_sg, grid3, geolon_c, geolat_c, geolon_t, geolat_t)


  end subroutine read_grid_file


  !#####################################################################
  subroutine read_topo_file

    integer :: fsize=65536
    integer :: status, ncid, id_var, ndim, dimsiz
    character(len=256) :: tile_file
    character(len=32)   :: text
    integer :: len, t, dims(2)
    real :: tmp(is:ie,js:je)

    allocate(oro(isd:ied,jsd:jed,ntiles))
    allocate(mask(isd:ied,jsd:jed,ntiles))
    oro = -big_number
    mask = 0

    !--- make sure topo_file suffix is not ".nc"
    len = len_trim(topo_file)
    if( index(topo_file, '.nc', back=.true.) == len-2) then
      call handle_err(-1, "remove .nc from namelist topo_file="//trim(topo_file) )
    endif

    !--- loop through each tile file to get the orography
    do t = 1, ntiles
       write(text, '(i1.1)' ) t
       tile_file = trim(topo_file)//'.tile'//trim(text)//'.nc'
       status=NF__OPEN(trim(tile_file),NF_NOWRITE,fsize,ncid)
       call handle_err(status, 'Open file '//trim(tile_file) )

       status=nf_inq_varid(ncid, topo_field, id_var)
       call handle_err(status, 'inquire varid of '//trim(topo_field)//' from file '//trim(tile_file) )

       status = nf_inq_varndims(ncid, id_var, ndim)
       call handle_err(status, 'inquire ndims of '//trim(topo_field)//' from file '//trim(tile_file) )

       if(ndim .NE. 2) call handle_err(-1, 'ndims of '//trim(topo_field)//' from file '// &
            trim(tile_file)//' should be 2')

       ! get data dimension and should match grid file size
       status = nf_inq_vardimid(ncid, id_var,dims);
       call handle_err(status, 'inquire dimid of '//trim(topo_field)//' from file '//trim(tile_file) )

       status = nf_inq_dimlen(ncid, dims(1), dimsiz)
       call handle_err(status, 'inquire first dimension length of '//trim(topo_field)//' from file '//trim(tile_file) )
       if(dimsiz .NE. nx) call handle_err(-1, "mismatch of lon dimension size between "// &
            trim(grid_file)//' and '//trim(tile_file) )

       status = nf_inq_dimlen(ncid, dims(2), dimsiz)
       call handle_err(status, 'inquire second dimension length of '//trim(topo_field)//' from file '//trim(tile_file) )

       if(dimsiz .NE. ny) call handle_err(-1, "mismatch of lat dimension size between "// &
            trim(grid_file)//' and '//trim(tile_file) )

       status = nf_get_var_double(ncid, id_var, oro(is:ie,js:je,t))
       call handle_err(status, 'get the value of '//trim(topo_field)//' from file '//trim(tile_file) )

       status=nf_inq_varid(ncid, mask_field, id_var)
       call handle_err(status, 'inquire varid of '//trim(mask_field)//' from file '//trim(tile_file) )

       status = nf_get_var_double(ncid, id_var, tmp)
       call handle_err(status, 'get the value of '//trim(mask_field)//' from file '//trim(tile_file) )

       mask(is:ie,js:je,t) = tmp

       status = nf_close(ncid)
       call handle_err(status, "close file "//trim(tile_file))
    enddo

    !--- update halo
    call fill_cubic_grid_halo(oro, oro, ng, 0, 0, 1, 1)
    call fill_cubic_grid_halo(mask, mask, ng, 0, 0, 1, 1)



  end subroutine read_topo_file

  !##############################################################################
  !--- replace the topo_field
  subroutine write_topo_file(is,ie,js,je,ntiles,q)
     integer, intent(in) :: is,ie,js,je,ntiles                 
     real,    intent(in) :: q(is:ie,js:je,ntiles) 

     integer :: fsize=65536
     integer :: t, status, ncid, id_var
     character(len=256) :: tile_file
    character(len=3)   :: text
     !--- loop through each tile file to update topo_field

    do t = 1, ntiles
       write(text, '(i1.1)' ) t
       tile_file = trim(topo_file)//'.tile'//trim(text)//'.nc'
       status=NF__OPEN(trim(tile_file),NF_WRITE,fsize,ncid)
       call handle_err(status, 'write_topo_file: Open file '//trim(tile_file) )

       status=nf_inq_varid(ncid, topo_field, id_var)
       call handle_err(status, 'write_topo_file:inquire varid of '//trim(topo_field)//' from file '//trim(tile_file) )

       status = nf_put_var_double(ncid, id_var, q(:,:,t))
       call handle_err(status, 'write_topo_file: put the value of '//trim(topo_field)//' from file '//trim(tile_file) )       

       status = nf_close(ncid)
       call handle_err(status, "write_topo_file: close file "//trim(tile_file))
    enddo


  end subroutine write_topo_file

  !##############################################################################
  ! this routine fill the halo points for the cubic grid. ioff and joff is used to distinguish
  ! T, C, E, or N-cell
  subroutine fill_cubic_grid_halo(data, data2, halo, ioff, joff, sign1, sign2)
    integer, intent(in)                               :: halo
    real, dimension(1-halo:,1-halo:,:), intent(inout) :: data, data2
    integer,                            intent(in)    :: ioff, joff, sign1, sign2 
    integer :: lw, le, ls, ln
    integer :: i, tile

    ntiles = size(data,3)

    do tile = 1, ntiles
       if(mod(tile,2) == 0) then ! tile 2, 4, 6
          lw = tile - 1; le = tile + 2; ls = tile - 2; ln = tile + 1
          if(le > 6 ) le = le - 6
          if(ls < 1 ) ls = ls + 6
          if(ln > 6 ) ln = ln - 6
          data(1-halo:0, 1:ny+joff, tile) = data(nx-halo+1:nx, 1:ny+joff, lw) ! west 
          do i = 1, halo 
             data(nx+i+ioff, 1:ny+joff, tile)    = sign1*data2(nx+joff:1:-1, i+ioff, le) ! east 
          end do
          do i = 1, halo 
             data(1:nx+ioff, 1-i, tile)     = sign2*data2(nx-i+1, ny+ioff:1:-1, ls) ! south 
          end do
          data(1:nx+ioff, ny+1+joff:ny+halo+joff, tile) = data(1:nx+ioff, 1+joff:halo+joff, ln) ! north
       else ! tile 1, 3, 5
          lw = tile - 2; le = tile + 1; ls = tile - 1; ln = tile + 2
          if(lw < 1 ) lw = lw + 6
          if(ls < 1 ) ls = ls + 6
          if(ln > 6 ) ln = ln - 6
          do i = 1, halo 
             data(1-i, 1:ny+joff, tile)     = sign1*data2(nx+joff:1:-1, ny-i+1, lw) ! west 
          end do
          data(nx+1+ioff:nx+halo+ioff, 1:ny+joff, tile) = data(1+ioff:halo+ioff, 1:ny+joff, le) ! east 
          data(1:nx+ioff, 1-halo:0, tile)     = data(1:nx+ioff, ny-halo+1:ny, ls) ! south 
          do i = 1, halo 
             data(1:nx+ioff, ny+i+joff, tile)    = sign2*data2(i+joff, ny+ioff:1:-1, ln) ! north 
          end do
       end if
    enddo

  end subroutine fill_cubic_grid_halo

  !#####################################################################
  subroutine FV3_zs_filter (is, ie, js, je, isd, ied, jsd, jed, npx, npy, npx_global, ntiles,  &
       grid_type, stretch_fac, nested, area, dxa, dya, dx, dy, dxc, dyc, &
       sin_sg,  phis )
    integer, intent(in) :: is, ie, js, je, ntiles
    integer, intent(in) :: isd, ied, jsd, jed, npx, npy, npx_global, grid_type
    real, intent(in), dimension(isd:ied,jsd:jed, ntiles)::area, dxa, dya
    real, intent(in), dimension(isd:ied,  jsd:jed+1, ntiles):: dx, dyc
    real, intent(in), dimension(isd:ied+1,jsd:jed, ntiles):: dy, dxc

    real, intent(IN):: sin_sg(4,isd:ied,jsd:jed,ntiles)
    real, intent(IN):: stretch_fac
    logical, intent(IN) :: nested
    real, intent(inout):: phis(isd:ied,jsd,jed,ntiles)
    real:: cd2
    integer mdim, n_del2, n_del4

    mdim = nint( real(npx_global) * min(10., stretch_fac) )

    ! Del-2: high resolution only
    if ( npx_global<=97 ) then
       n_del2 = 0
    elseif ( npx_global<=193 ) then
       n_del2 = 1
    else
       n_del2 = 2
    endif
    cd2 = 0.16*da_min
    ! Applying strong 2-delta-filter:
    if ( n_del2 > 0 )   &
         call two_delta_filter(is,ie,js,je,isd,ied,jsd,jed, npx, npy, ntiles, phis, area, &
                               dx, dy, dxa, dya, dxc, dyc, sin_sg, cd2, zero_ocean,  &
                               .true.,0, grid_type, mask, nested, n_del2)

    ! MFCT Del-4:
    if ( mdim<=193 ) then
       n_del4 = 1
    elseif ( mdim<=1537 ) then
       n_del4 = 2
    else
       n_del4 = 3
    endif
    call del4_cubed_sphere(is,ie,js,je,isd,ied,jsd,jed,npx, npy, ntiles, &
          phis, area, dx, dy, dxc, dyc, sin_sg, n_del4, zero_ocean, mask, nested)
    ! Applying weak 2-delta-filter:
    cd2 = 0.12*da_min
    call two_delta_filter(is,ie,js,je,isd,ied,jsd,jed,npx, npy, ntiles, &
          phis, area, dx, dy, dxa, dya, dxc, dyc, sin_sg, cd2, zero_ocean,  &
          .true., 1, grid_type, mask, nested, n_del2_weak)

  end subroutine FV3_zs_filter

  !#####################################################################
  subroutine two_delta_filter(is, ie, js, je, isd, ied, jsd, jed, npx, npy, ntiles, &
       q, area, dx, dy, dxa, dya, dxc, dyc, sin_sg, cd, zero_ocean,  &
        check_slope, filter_type, grid_type, mask, nested, ntmax)
    integer, intent(in) :: is,  ie,  js,  je
    integer, intent(in) :: isd, ied, jsd, jed
    integer, intent(in) :: npx, npy, grid_type
    integer, intent(in) :: ntmax, ntiles
    integer, intent(in) :: filter_type    ! 0: strong,   1: weak
    real,    intent(in) :: cd
    ! INPUT arrays
    real, intent(in)::area(isd:ied,  jsd:jed, ntiles)
    real, intent(in)::  dx(isd:ied,  jsd:jed+1, ntiles)
    real, intent(in)::  dy(isd:ied+1,jsd:jed, ntiles)
    real, intent(in):: dxa(isd:ied,  jsd:jed, ntiles)
    real, intent(in):: dya(isd:ied,  jsd:jed, ntiles)
    real, intent(in):: dxc(isd:ied+1,jsd:jed, ntiles)
    real, intent(in):: dyc(isd:ied,  jsd:jed+1, ntiles)
    real, intent(in):: sin_sg(4,isd:ied,jsd:jed, ntiles)
    real, intent(in):: mask(isd:ied,  jsd:jed, ntiles)        ! 0==water, 1==land
    logical, intent(in):: zero_ocean, check_slope
    logical, intent(in):: nested
    ! OUTPUT arrays
    real, intent(inout):: q(isd:ied, jsd:jed,ntiles)
    ! Local:
    real, parameter:: p1 =  7./12.
    real, parameter:: p2 = -1./12.
    real, parameter:: c1 = -2./14.
    real, parameter:: c2 = 11./14.
    real, parameter:: c3 =  5./14.

    real:: ddx(is:ie+1,js:je), ddy(is:ie,js:je+1)
    logical:: extm(is-1:ie+1)
    logical:: ext2(is:ie,js-1:je+1)
    real::  a1(is-1:ie+2)
    real::  a2(is:ie,js-1:je+2)
    real::  a3(is:ie,js:je,ntiles)
    real:: smax, smin, m_slope, fac
    integer:: i,j, nt, t
    integer:: is1, ie2, js1, je2


    if ( .not. nested .and. grid_type<3 ) then
       is1 = max(3,is-1);  ie2 = min(npx-2,ie+2)
       js1 = max(3,js-1);  je2 = min(npy-2,je+2)
    else
       is1 = is-1;         ie2 = ie+2
       js1 = js-1;         je2 = je+2
    end if

    if ( check_slope ) then
        m_slope = max_slope
    else
        m_slope = 10.
    endif
         

    do nt=1, ntmax
       call fill_cubic_grid_halo(q, q, ng, 0, 0, 1, 1)

       ! Check slope
       if ( nt==1 .and. check_slope ) then
          do t = 1, ntiles
             do j=js,je
                do i=is,ie+1
                   ddx(i,j) = (q(i,j,t) - q(i-1,j,t))/dxc(i,j,t) 
                   ddx(i,j) = abs(ddx(i,j))
                enddo
             enddo
             do j=js,je+1
                do i=is,ie
                   ddy(i,j) = (q(i,j,t) - q(i,j-1,t))/dyc(i,j,t) 
                   ddy(i,j) = abs(ddy(i,j))
                enddo
             enddo
             do j=js,je
                do i=is,ie
                   a3(i,j,t) = max( ddx(i,j), ddx(i+1,j), ddy(i,j), ddy(i,j+1) )
                enddo
             enddo
          enddo
          smax = maxval(a3(is:ie,js:je,:))
          write(*,*) 'Before filter: Max_slope=', smax
       endif


       ! First step: average the corners:
       if ( .not. nested .and. nt==1 ) then
          do t = 1, ntiles
             q(1,1,t) = (q(1,1,t)*area(1,1,t)+q(0,1,t)*area(0,1,t)+q(1,0,t)*area(1,0,t))  &
                  / (       area(1,1,t)+       area(0,1,t)+       area(1,0,t) )
             q(0,1,t) =  q(1,1,t)
             q(1,0,t) =  q(1,1,t)

             q(ie, 1,t) = (q(ie,1,t)*area(ie,1,t)+q(npx,1,t)*area(npx,1,t)+q(ie,0,t)*area(ie,0,t)) &
                  / (        area(ie,1,t)+         area(npx,1,t)+        area(ie,0,t))
             q(npx,1,t) =  q(ie,1,t)
             q(ie, 0,t) =  q(ie,1,t)

             q(1, je,t) = (q(1,je,t)*area(1,je,t)+q(0,je,t)*area(0,je,t)+q(1,npy,t)*area(1,npy,t))   &
                  / (        area(1,je,t)+        area(0,je,t)+         area(1,npy,t))
             q(0, je,t) =  q(1,je,t)
             q(1,npy,t) =  q(1,je,t)

             q(ie, je,t) = (q(ie,je,t)*area(ie,je,t)+q(npx,je,t)*area(npx,je,t)+q(ie,npy,t)*area(ie,npy,t))  &
                  / (         area(ie,je,t)+          area(npx,je,t)+          area(ie,npy,t))
             q(npx,je,t) =  q(ie,je,t)
             q(ie,npy,t) =  q(ie,je,t)
          enddo
          call fill_cubic_grid_halo(q, q, ng, 0, 0, 1, 1)
       endif

       do t = 1, ntiles
          ! x-diffusive flux:
          do j=js,je

             do i=is1, ie2
                a1(i) = p1*(q(i-1,j,t)+q(i,j,t)) + p2*(q(i-2,j,t)+q(i+1,j,t))
             enddo

             if ( .not. nested .and. grid_type<3 ) then
                a1(0) = c1*q(-2,j,t) + c2*q(-1,j,t) + c3*q(0,j,t)
                a1(1) = 0.5*(((2.*dxa(0,j,t)+dxa(-1,j,t))*q(0,j,t)-dxa(0,j,t)*q(-1,j,t))/(dxa(-1,j,t)+dxa(0,j,t)) &
                     +      ((2.*dxa(1,j,t)+dxa( 2,j,t))*q(1,j,t)-dxa(1,j,t)*q( 2,j,t))/(dxa(1, j,t)+dxa(2,j,t)))
                a1(2) = c3*q(1,j,t) + c2*q(2,j,t) +c1*q(3,j,t)

                a1(npx-1) = c1*q(npx-3,j,t) + c2*q(npx-2,j,t) + c3*q(npx-1,j,t)
                a1(npx) = 0.5*(((2.*dxa(npx-1,j,t)+dxa(npx-2,j,t))*q(npx-1,j,t)-dxa(npx-1,j,t)*q(npx-2,j,t)) &
                          /(dxa(npx-2,j,t)+dxa(npx-1,j,t)) &
                     +      ((2.*dxa(npx,  j,t)+dxa(npx+1,j,t))*q(npx,  j,t)-dxa(npx,  j,t)*q(npx+1,j,t))/ &
                          (dxa(npx,  j,t)+dxa(npx+1,j,t)))
                a1(npx+1) = c3*q(npx,j,t) + c2*q(npx+1,j,t) + c1*q(npx+2,j,t)
             endif

             if ( filter_type == 0 ) then
                do i=is-1, ie+1
                   if( abs(3.*(a1(i)+a1(i+1)-2.*q(i,j,t))) > abs(a1(i)-a1(i+1)) ) then
                      extm(i) = .true.
                   else
                      extm(i) = .false.
                   endif
                enddo
             else
                do i=is-1, ie+1
                   if ( (a1(i)-q(i,j,t))*(a1(i+1)-q(i,j,t)) > 0. ) then
                      extm(i) = .true.
                   else
                      extm(i) = .false.
                   endif
                enddo
             endif

             do i=is,ie+1
                ddx(i,j) = (q(i-1,j,t)-q(i,j,t))/dxc(i,j,t)
                if ( extm(i-1).and.extm(i) ) then
                   ddx(i,j) = 0.5*(sin_sg(3,i-1,j,t)+sin_sg(1,i,j,t))*dy(i,j,t)*ddx(i,j)
                elseif ( abs(ddx(i,j)) > m_slope ) then
                   fac = min(1., max( 0.1, ( abs(ddx(i,j))-m_slope )/m_slope) )
                   ddx(i,j) = fac*0.5*(sin_sg(3,i-1,j,t)+sin_sg(1,i,j,t))*dy(i,j,t)*ddx(i,j)
                else
                   ddx(i,j) = 0.
                endif
             enddo
          enddo ! do j=js,je

          ! y-diffusive flux:
          do j=js1,je2
             do i=is,ie
                a2(i,j) = p1*(q(i,j-1,t)+q(i,j,t)) + p2*(q(i,j-2,t)+q(i,j+1,t))
             enddo
          enddo
          if ( .not. nested .and. grid_type<3 ) then
             do i=is,ie
                a2(i,0) = c1*q(i,-2,t) + c2*q(i,-1,t) + c3*q(i,0,t)
                a2(i,1) = 0.5*(((2.*dya(i,0,t)+dya(i,-1,t))*q(i,0,t)-dya(i,0,t)*q(i,-1,t))/(dya(i,-1,t)+dya(i,0,t))   &
                     +      ((2.*dya(i,1,t)+dya(i, 2,t))*q(i,1,t)-dya(i,1,t)*q(i, 2,t))/(dya(i, 1,t)+dya(i,2,t)))
                a2(i,2) = c3*q(i,1,t) + c2*q(i,2,t) + c1*q(i,3,t)
             enddo

             do i=is,ie
                a2(i,npy-1) = c1*q(i,npy-3,t) + c2*q(i,npy-2,t) + c3*q(i,npy-1,t)
                a2(i,npy) = 0.5*(((2.*dya(i,npy-1,t)+dya(i,npy-2,t))*q(i,npy-1,t)-dya(i,npy-1,t)*q(i,npy-2,t))/ &
                            (dya(i,npy-2,t)+dya(i,npy-1,t))  &
                     +      ((2.*dya(i,npy,t)+dya(i,npy+1,t))*q(i,npy,t)-dya(i,npy,t)*q(i,npy+1,t))/&
                           (dya(i,npy,t)+dya(i,npy+1,t)))
                a2(i,npy+1) = c3*q(i,npy,t) + c2*q(i,npy+1,t) + c1*q(i,npy+2,t)
             enddo
          endif

          if ( filter_type == 0 ) then
             do j=js-1,je+1
                do i=is,ie
                   if( abs(3.*(a2(i,j)+a2(i,j+1)-2.*q(i,j,t))) > abs(a2(i,j)-a2(i,j+1)) ) then
                      ext2(i,j) = .true.
                   else
                      ext2(i,j) = .false.
                   endif
                enddo
             enddo
          else
             do j=js-1,je+1
                do i=is,ie
                   if ( (a2(i,j)-q(i,j,t))*(a2(i,j+1)-q(i,j,t)) > 0. ) then
                      ext2(i,j) = .true.
                   else
                      ext2(i,j) = .false.
                   endif
                enddo
             enddo
          endif

          do j=js,je+1
             do i=is,ie
                ddy(i,j) = (q(i,j-1,t)-q(i,j,t))/dyc(i,j,t)
                if ( ext2(i,j-1) .and. ext2(i,j) ) then
                   ddy(i,j) = 0.5*(sin_sg(4,i,j-1,t)+sin_sg(2,i,j,t))*dx(i,j,t)*ddy(i,j)
                elseif ( abs(ddy(i,j))>m_slope ) then
                   fac = min(1., max(0.1,(abs(ddy(i,j))-m_slope)/m_slope))
                   ddy(i,j) = fac*0.5*(sin_sg(4,i,j-1,t)+sin_sg(2,i,j,t))*dx(i,j,t)*ddy(i,j)
                else
                   ddy(i,j) = 0.
                endif
             enddo
          enddo

          if ( zero_ocean ) then
             ! Limit diffusive flux over water cells:
             do j=js,je
                do i=is,ie+1
                   ddx(i,j) = max(0., min(mask(i-1,j,t), mask(i,j,t))) * ddx(i,j)
                enddo
             enddo
             do j=js,je+1
                do i=is,ie
                   ddy(i,j) = max(0., min(mask(i,j-1,t), mask(i,j,t))) * ddy(i,j)
                enddo
             enddo
          endif

          do j=js,je
             do i=is,ie
                q(i,j,t) = q(i,j,t) + cd/area(i,j,t)*(ddx(i,j)-ddx(i+1,j)+ddy(i,j)-ddy(i,j+1))
             enddo
          enddo
       enddo ! do t = 1, ntiles
    enddo ! nt=1, ntmax

! Check slope
    if ( check_slope ) then
         call fill_cubic_grid_halo(q, q, ng, 0, 0, 1, 1)
         do t = 1, ntiles
            do j=js,je
               do i=is,ie+1
                  ddx(i,j) = (q(i,j,t) - q(i-1,j,t))/dxc(i,j,t) 
                  ddx(i,j) = abs(ddx(i,j))
               enddo
            enddo
            do j=js,je+1
               do i=is,ie
                  ddy(i,j) = (q(i,j,t) - q(i,j-1,t))/dyc(i,j,t) 
                  ddy(i,j) = abs(ddy(i,j))
               enddo
            enddo
            do j=js,je
               do i=is,ie
                  a3(i,j,t) = max( ddx(i,j), ddx(i+1,j), ddy(i,j), ddy(i,j+1) )
               enddo
            enddo
         enddo
         smax = maxval(a3(is:ie,js:je,:))
         write(*,*) 'After filter: Max_slope=', smax
    endif

  end subroutine two_delta_filter

  !#####################################################################
  subroutine del2_cubed_sphere(is, ie, js, je, isd, ied, jsd, jed, npx, npy, ntiles,&
          q, area, dx, dy, dxc, dyc, sin_sg, nmax, cd, zero_ocean, mask, nested)
    integer, intent(in) :: is,  ie,  js,  je
    integer, intent(in) :: isd, ied, jsd, jed
    integer, intent(in):: npx, npy, ntiles
    integer, intent(in):: nmax
    real, intent(in):: cd
    logical, intent(in):: zero_ocean
    ! INPUT arrays
    real, intent(in)::area(isd:ied,  jsd:jed, ntiles)
    real, intent(in)::  dx(isd:ied,  jsd:jed+1, ntiles)
    real, intent(in)::  dy(isd:ied+1,jsd:jed, ntiles)
    real, intent(in):: dxc(isd:ied+1,jsd:jed, ntiles)
    real, intent(in):: dyc(isd:ied,  jsd:jed+1, ntiles)
    real, intent(IN):: sin_sg(4,isd:ied,jsd:jed, ntiles)
    real, intent(in):: mask(isd:ied,  jsd:jed, ntiles)        ! 0==water, 1==land
    logical, intent(IN) :: nested

    ! OUTPUT arrays
    real, intent(inout):: q(is-ng:ie+ng, js-ng:je+ng, ntiles)
    ! Local:
    real ddx(is:ie+1,js:je), ddy(is:ie,js:je+1)
    integer i,j,n,t

    call fill_cubic_grid_halo(q, q, ng, 0, 0, 1, 1)

    do t = 1, ntiles
       ! First step: average the corners:
       if ( .not. nested) then
          q(1,1,t) = (q(1,1,t)*area(1,1,t)+q(0,1,t)*area(0,1,t)+q(1,0,t)*area(1,0,t))  &
               / (       area(1,1,t)+       area(0,1,t)+       area(1,0,t) )
          q(0,1,t) =  q(1,1,t)
          q(1,0,t) =  q(1,1,t)
       endif
       if ( .not. nested) then
          q(ie, 1,t) = (q(ie,1,t)*area(ie,1,t)+q(npx,1,t)*area(npx,1,t)+q(ie,0,t)*area(ie,0,t)) &
               / (        area(ie,1,t)+         area(npx,1,t)+        area(ie,0,t))
          q(npx,1,t) =  q(ie,1,t)
          q(ie, 0,t) =  q(ie,1,t)
       endif
       if ( .not. nested ) then
          q(ie, je,t) = (q(ie,je,t)*area(ie,je,t)+q(npx,je,t)*area(npx,je,t)+q(ie,npy,t)*area(ie,npy,t))  &
               / (         area(ie,je,t)+          area(npx,je,t)+          area(ie,npy,t))
          q(npx,je,t) =  q(ie,je,t)
          q(ie,npy,t) =  q(ie,je,t)
       endif
       if ( .not. nested) then
          q(1, je,t) = (q(1,je,t)*area(1,je,t)+q(0,je,t)*area(0,je,t)+q(1,npy,t)*area(1,npy,t))   &
               / (        area(1,je,t)+        area(0,je,t)+         area(1,npy,t))
          q(0, je,t) =  q(1,je,t)
          q(1,npy,t) =  q(1,je,t)
       endif
    enddo

    do n=1,nmax
       if( n>1 ) call fill_cubic_grid_halo(q, q, ng, 0, 0, 1, 1)
       do t = 1, ntiles
          do j=js,je
             do i=is,ie+1
                ddx(i,j) = 0.5*(sin_sg(3,i-1,j,t)+sin_sg(1,i,j,t))*dy(i,j,t)*(q(i-1,j,t)-q(i,j,t))/dxc(i,j,t)
             enddo
          enddo
          do j=js,je+1
             do i=is,ie
                ddy(i,j) = dx(i,j,t)*(q(i,j-1,t)-q(i,j,t))/dyc(i,j,t) &
                     *0.5*(sin_sg(4,i,j-1,t)+sin_sg(2,i,j,t))
             enddo
          enddo

          if ( zero_ocean ) then
             ! Limit diffusive flux over ater cells:
             do j=js,je
                do i=is,ie+1
                   ddx(i,j) = max(0., min(mask(i-1,j,t), mask(i,j,t))) * ddx(i,j)
                enddo
             enddo
             do j=js,je+1
                do i=is,ie
                   ddy(i,j) = max(0., min(mask(i,j-1,t), mask(i,j,t))) * ddy(i,j)
                enddo
             enddo
          endif

          do j=js,je
             do i=is,ie
                q(i,j,t) = q(i,j,t) + cd/area(i,j,t)*(ddx(i,j)-ddx(i+1,j)+ddy(i,j)-ddy(i,j+1))
             enddo
          enddo
       enddo
    enddo

  end subroutine del2_cubed_sphere

  !#####################################################################
  subroutine del4_cubed_sphere(is, ie, js, je, isd, ied, jsd, jed, npx, npy, ntiles, &
       q, area, dx, dy, dxc, dyc, sin_sg, nmax, zero_ocean, mask, nested)
    integer, intent(in) :: is,  ie,  js,  je
    integer, intent(in) :: isd, ied, jsd, jed
    integer, intent(in) :: npx, npy, nmax, ntiles
    logical, intent(in) :: zero_ocean
    real, intent(in):: mask(isd:ied,  jsd:jed, ntiles)        ! 0==water, 1==land
    real, intent(in)::area(isd:ied,  jsd:jed, ntiles)
    real, intent(in)::  dx(isd:ied,  jsd:jed+1, ntiles)
    real, intent(in)::  dy(isd:ied+1,jsd:jed, ntiles)
    real, intent(in):: dxc(isd:ied+1,jsd:jed, ntiles)
    real, intent(in):: dyc(isd:ied,  jsd:jed+1, ntiles)
    real, intent(IN):: sin_sg(4,isd:ied,jsd:jed, ntiles)
    real, intent(inout):: q(isd:ied, jsd:jed, ntiles)
    logical, intent(IN) :: nested
    ! Local:
    ! diffusivity
    real :: diff(is-1:ie+1,js-1:je+1, ntiles)
    ! diffusive fluxes: 
    real :: fx1(is:ie+1,js:je), fy1(is:ie,js:je+1)
    real :: fx2(is:ie+1,js:je,ntiles), fy2(is:ie,js:je+1,ntiles)
    real :: fx4(is:ie+1,js:je,ntiles), fy4(is:ie,js:je+1,ntiles)
    real, dimension(isd:ied,jsd:jed,ntiles):: d2, win, wou 
    real, dimension(is:ie,js:je, ntiles) :: qlow, qmin, qmax
    real, parameter:: esl = 1.E-20
    integer i,j, n, t

    ! On a nested grid the haloes are not filled. Set to zero.
    d2 = 0.
    win = 0.
    wou = 0.

    do t = 1, ntiles
       do j=js-1,je+1 ; do i=is-1,ie+1
          diff(i,j,t) = cd4*area(i,j,t) ! area dependency is needed for stretched grid
       enddo; enddo

       do j=js,je ; do i=is,ie
          qmax(i,j,t) = q(i,j,t) * peak_fac
          qmin(i,j,t) = q(i,j,t) / peak_fac
       enddo; enddo
    enddo

    do n=1,nmax
       call fill_cubic_grid_halo(q, q, ng, 0, 0, 1, 1)

       ! First step: average the corners:
       if ( .not. nested .and. n==1 ) then
          do t = 1, ntiles
             q(1,1,t) = (q(1,1,t)*area(1,1,t)+q(0,1,t)*area(0,1,t)+q(1,0,t)*area(1,0,t))  &
                  / (       area(1,1,t)+       area(0,1,t)+       area(1,0,t) )
             q(0,1,t) = q(1,1,t)
             q(1,0,t) = q(1,1,t)
             q(0,0,t) = q(1,1,t)

             q(ie, 1,t) = (q(ie,1,t)*area(ie,1,t)+q(npx,1,t)*area(npx,1,t)+q(ie,0,t)*area(ie,0,t)) &
                  / (        area(ie,1,t)+         area(npx,1,t)+        area(ie,0,t))
             q(npx,1,t) = q(ie,1,t)
             q(ie, 0,t) = q(ie,1,t)
             q(npx,0,t) = q(ie,1,t)

             q(1, je,t) = (q(1,je,t)*area(1,je,t)+q(0,je,t)*area(0,je,t)+q(1,npy,t)*area(1,npy,t))   &
                  / (        area(1,je,t)+        area(0,je,t)+         area(1,npy,t))
             q(0, je,t) =  q(1,je,t)
             q(1,npy,t) =  q(1,je,t)
             q(0,npy,t) =  q(1,je,t)

             q(ie, je,t) = (q(ie,je,t)*area(ie,je,t)+q(npx,je,t)*area(npx,je,t)+q(ie,npy,t)*area(ie,npy,t))  &
                  / (         area(ie,je,t)+          area(npx,je,t)+          area(ie,npy,t))
             q(npx, je,t) = q(ie,je,t)
             q(ie, npy,t) = q(ie,je,t)
             q(npx,npy,t) = q(ie,je,t)
          enddo
          call fill_cubic_grid_halo(q, q, ng, 0, 0, 1, 1)
       endif

       do t = 1, ntiles

          !--------------
          ! Compute del-2
          !--------------
          !     call copy_corners(q, npx, npy, 1)
          do j=js,je
             do i=is,ie+1
                fx2(i,j,t) = 0.25*(diff(i-1,j,t)+diff(i,j,t))*dy(i,j,t)*(q(i-1,j,t)-q(i,j,t))/dxc(i,j,t)          &
                     *(sin_sg(1,i,j,t)+sin_sg(3,i-1,j,t))
             enddo
          enddo

          !     call copy_corners(q, npx, npy, 2)
          do j=js,je+1
             do i=is,ie
                fy2(i,j,t) = 0.25*(diff(i,j-1,t)+diff(i,j,t))*dx(i,j,t)*(q(i,j-1,t)-q(i,j,t))/dyc(i,j,t) &
                     *(sin_sg(2,i,j,t)+sin_sg(4,i,j-1,t))
             enddo
          enddo

          do j=js,je
             do i=is,ie
                d2(i,j,t) = (fx2(i,j,t)-fx2(i+1,j,t)+fy2(i,j,t)-fy2(i,j+1,t)) / area(i,j,t)
             enddo
          enddo

          ! qlow == low order monotonic solution
          if ( zero_ocean ) then
             ! Limit diffusive flux over water cells:
             do j=js,je
                do i=is,ie+1
                   fx1(i,j) = max(0., min(mask(i-1,j,t), mask(i,j,t))) * fx2(i,j,t)
                enddo
             enddo
             do j=js,je+1
                do i=is,ie
                   fy1(i,j) = max(0., min(mask(i,j-1,t), mask(i,j,t))) * fy2(i,j,t)
                enddo
             enddo
             do j=js,je
                do i=is,ie
                   qlow(i,j,t) = q(i,j,t) + (fx1(i,j)-fx1(i+1,j)+fy1(i,j)-fy1(i,j+1)) / area(i,j,t)
                   d2(i,j,t) = diff(i,j,t) * d2(i,j,t)
                enddo
             enddo
          else
             do j=js,je
                do i=is,ie
                   qlow(i,j,t) =    q(i,j,t) + d2(i,j,t)
                   d2(i,j,t) = diff(i,j,t) * d2(i,j,t)
                enddo
             enddo
          endif
       enddo
       call fill_cubic_grid_halo(d2, d2, ng, 0, 0, 1, 1)

       !---------------------
       ! Compute del4 fluxes:
       !---------------------
       !     call copy_corners(d2, npx, npy, 1)
       do t = 1, ntiles
          do j=js,je
             do i=is,ie+1
                fx4(i,j,t) = 0.5*(sin_sg(3,i-1,j,t)+sin_sg(1,i,j,t))*dy(i,j,t)*(d2(i,j,t)-d2(i-1,j,t))/dxc(i,j,t)-fx2(i,j,t)
             enddo
          enddo

          !     call copy_corners(d2, npx, npy, 2)
          do j=js,je+1
             do i=is,ie
                fy4(i,j,t) = dx(i,j,t)*(d2(i,j,t)-d2(i,j-1,t))/dyc(i,j,t) &
                     *0.5*(sin_sg(2,i,j,t)+sin_sg(4,i,j-1,t))-fy2(i,j,t)
             enddo
          enddo

          do j=js,je
             do i=is,ie
                qmin(i,j,t) = min(qmin(i,j,t), q(i-1,j-1,t), q(i,j-1,t), q(i+1,j-1,t),  &
                                               q(i-1,j  ,t), q(i,j  ,t), q(i+1,j  ,t),  &
                                               q(i-1,j+1,t), q(i,j+1,t), q(i+1,j+1,t) )
                qmax(i,j,t) = max(qmax(i,j,t), q(i-1,j-1,t), q(i,j-1,t), q(i+1,j-1,t),  &
                                               q(i-1,j  ,t), q(i,j  ,t), q(i+1,j  ,t),  &
                                               q(i-1,j+1,t), q(i,j+1,t), q(i+1,j+1,t) )
             enddo
          enddo

          !----------------
          ! Flux limitting:
          !----------------
          do j=js,je
             do i=is,ie
                win(i,j,t) = max(0.,fx4(i,  j,t)) - min(0.,fx4(i+1,j,t)) +   &
                     max(0.,fy4(i,  j,t)) - min(0.,fy4(i,j+1,t)) + esl
                wou(i,j,t) = max(0.,fx4(i+1,j,t)) - min(0.,fx4(i,  j,t)) +   &
                     max(0.,fy4(i,j+1,t)) - min(0.,fy4(i,  j,t)) + esl
                win(i,j,t) = max(0., qmax(i,j,t) - qlow(i,j,t)) / win(i,j,t)*area(i,j,t)
                wou(i,j,t) = max(0., qlow(i,j,t) - qmin(i,j,t)) / wou(i,j,t)*area(i,j,t)
             enddo
          enddo
       enddo
       call fill_cubic_grid_halo(win, win, ng, 0, 0, 1, 1)
       call fill_cubic_grid_halo(wou, wou, ng, 0, 0, 1, 1)
       do t = 1, ntiles
          do j=js,je
             do i=is,ie+1
                if ( fx4(i,j,t) > 0. ) then
                   fx4(i,j,t) = min(1., wou(i-1,j,t), win(i,j,t)) * fx4(i,j,t) 
                else
                   fx4(i,j,t) = min(1., win(i-1,j,t), wou(i,j,t)) * fx4(i,j,t) 
                endif
             enddo
          enddo
          do j=js,je+1
             do i=is,ie
                if ( fy4(i,j,t) > 0. ) then
                   fy4(i,j,t) = min(1., wou(i,j-1,t), win(i,j,t)) * fy4(i,j,t) 
                else
                   fy4(i,j,t) = min(1., win(i,j-1,t), wou(i,j,t)) * fy4(i,j,t) 
                endif
             enddo
          enddo


          if ( zero_ocean ) then
             ! Limit diffusive flux over ocean cells:
             do j=js,je
                do i=is,ie+1
                   fx4(i,j,t) = max(0., min(mask(i-1,j,t), mask(i,j,t))) * fx4(i,j,t)
                enddo
             enddo
             do j=js,je+1
                do i=is,ie
                   fy4(i,j,t) = max(0., min(mask(i,j-1,t), mask(i,j,t))) * fy4(i,j,t)
                enddo
             enddo
          endif

          ! Update:
          do j=js,je
             do i=is,ie
                q(i,j,t) = qlow(i,j,t) + (fx4(i,j,t)-fx4(i+1,j,t)+fy4(i,j,t)-fy4(i,j+1,t))/area(i,j,t)
             enddo
          enddo
       enddo
    enddo    ! end n-loop


  end subroutine del4_cubed_sphere

  !#####################################################################
  subroutine handle_err(status, string)
    integer,          intent(in) :: status
    character(len=*), intent(in) :: string
    character(len=256) :: errmsg

    if (status .ne. nf_noerr) then
       errmsg = nf_strerror(status)
       errmsg = trim(errmsg)//trim(string)
       print *, trim(errmsg)
       stop 'Stopped'
    endif

  end subroutine  handle_err


  !#######################################################################
  ! reads the namelist file, write namelist to log file,

  subroutine read_namelist

    !  read namelist
    integer :: unit=7, io_status
    logical :: opened

    do
       inquire( unit=unit, opened=opened )
       if( .NOT.opened )exit
       unit = unit + 1
       if( unit.EQ.100 )call handle_err(-1, 'Unable to locate unit number.' )
    end do
    open( unit=unit, file='input.nml', iostat=io_status )
    read( unit,filter_topo_nml, iostat=io_status )
    close(unit)

    if (io_status > 0) call handle_err(-1, 'Error reading input.nml')

    write (stdunit, nml=filter_topo_nml)  

  end subroutine read_namelist


end program filter_topo
