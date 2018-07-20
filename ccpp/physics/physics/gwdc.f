!> \file gwdc.f This file is the original code for parameterization of
!! stationary convection forced gravity wave drag based on 
!! \cite chun_and_baik_1998 .

      module gwdc_pre
      contains

! \brief Brief description of the subroutine
!
!> \section arg_table_gwdc_pre_init Argument Table
!!
      subroutine gwdc_pre_init()
      end subroutine gwdc_pre_init

! \brief Brief description of the subroutine
!!
!! \section arg_table_gwdc_pre_run Argument Table
!! | local_name     | standard_name                                                            | long_name                                                     | units | rank | type      | kind      | intent | optional |
!! |----------------|--------------------------------------------------------------------------|---------------------------------------------------------------|-------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                   | horizontal loop extent                                        | count |    0 | integer   |           | in     | F        |
!! | cgwf           | multiplication_factors_for_convective_gravity_wave_drag                  | multiplication factors for convective gravity wave drag       | none  |    1 | real      | kind_phys | in     | F        |
!! | dx             | cell_size                                                                | grid size in zonal direction                                  | m     |    1 | real      | kind_phys | in     | F        |
!! | work1          | grid_size_related_coefficient_used_in_scale-sensitive_schemes            | grid size related coefficient used in scale-sensitive schemes | none  |    1 | real      | kind_phys | in     | F        |
!! | work2          | grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement | complement to work1                                           | none  |    1 | real      | kind_phys | in     | F        |
!! | dlength        | characteristic_grid_length_scale                                         | representative horizontal length scale of grid box            | m     |    1 | real      | kind_phys | out    | F        |
!! | cldf           | cloud_area_fraction                                                      | fraction of grid box area in which updrafts occur             | frac  |    1 | real      | kind_phys | out    | F        |
!! | levs           | vertical_dimension                                                       | number of vertical layers                                     | count |    0 | integer   |           | in     | F        |
!! | kbot           | vertical_index_at_cloud_base                                             | vertical index at cloud base                                  | index |    1 | integer   |           | in     | F        |
!! | ktop           | vertical_index_at_cloud_top                                              | vertical index at cloud top                                   | index |    1 | integer   |           | in     | F        |
!! | dtp            | time_step_for_physics                                                    | physics time step                                             | s     |    0 | real      | kind_phys | in     | F        |
!! | gt0            | air_temperature_updated_by_physics                                       | updated air temperature                                       | K     |    2 | real      | kind_phys | in     | F        |
!! | gt0_init       | air_temperature_save                                                     | air temperature before entering convection scheme             | K     |    2 | real      | kind_phys | in     | F        |
!! | del            | air_pressure_difference_between_midlayers                                | difference between mid-layer pressures                        | Pa    |    2 | real      | kind_phys | in     | F        |
!! | cumabs         | maximum_column_heating_rate                                              | maximum heating rate in column                                | K s-1 |    1 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                                            | error message for error handling in CCPP                      | none  |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                               | error flag for error handling in CCPP                         | flag  |    0 | integer   |           | out    | F        |
!!
      subroutine gwdc_pre_run (                                         &
     &  im, cgwf, dx, work1, work2, dlength, cldf,                      &
     &  levs, kbot, ktop, dtp, gt0, gt0_init, del, cumabs,              &
     &  errmsg, errflg )

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs
      integer, intent(in) :: kbot(:), ktop(:)
      real(kind=kind_phys), intent(in) :: dtp
      real(kind=kind_phys), intent(in) :: cgwf(2)
      real(kind=kind_phys), intent(in) :: dx(:), work1(:), work2(:)
      real(kind=kind_phys), intent(in) ::                               &
     &  gt0(:,:), gt0_init(:,:), del(:,:)

      real(kind=kind_phys), intent(out) ::                              &
     &  dlength(:), cldf(:), cumabs(:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k
      real(kind=kind_phys) :: tem1, tem2
      real(kind=kind_phys) :: work3(im)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i = 1, im
        tem1       = dx(i)
        tem2       = tem1
        dlength(i) = sqrt( tem1*tem1+tem2*tem2 )
        cldf(i)    = cgwf(1)*work1(i) + cgwf(2)*work2(i)
      enddo

!  --- ...  calculate maximum convective heating rate
!           cuhr = temperature change due to deep convection

      cumabs(:) = 0.0
      work3(:)  = 0.0
      do k = 1, levs
        do i = 1, im
          if (k >= kbot(i) .and. k <= ktop(i)) then
            cumabs(i)                                                   &
     &        = cumabs(i) + (gt0(i,k) - gt0_init(i,k)) * del(i,k)
            work3(i)  = work3(i)  + del(i,k)
          endif
        enddo
      enddo
      do i=1,im
        if (work3(i) > 0.0) cumabs(i) = cumabs(i) / (dtp*work3(i))
      enddo

      end subroutine gwdc_pre_run

! \brief Brief description of the subroutine
!
!> \section arg_table_gwdc_pre_finalize Argument Table
!!
      subroutine gwdc_pre_finalize ()
      end subroutine gwdc_pre_finalize

      end module gwdc_pre



      module gwdc

      contains

! \brief Brief description of the subroutine
!
!> \section arg_table_gwdc_init Argument Table
!!
      subroutine gwdc_init()
      end subroutine gwdc_init

! \defgroup GFS_cgwd GFS Convective Gravity Wave Drag
!> \defgroup GFS_gwdc_run GFS gwdc Main
!! \brief This subroutine is the parameterization of convective gravity wave
!! drag based on the theory given by 
!! \cite chun_and_baik_1998 modified for implementation into the
!! GFS/CFS by Ake Johansson(Aug 2005).
!!
!> \section arg_table_gwdc_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                          | units      | rank | type      | kind      | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                 | horizontal loop extent                                             | count      |    0 | integer   |           | in     | F        |
!! | ix             | horizontal_dimension                                   | horizontal dimension                                               | count      |    0 | integer   |           | in     | F        |
!! | km             | vertical_dimension                                     | number of vertical layers                                          | count      |    0 | integer   |           | in     | F        |
!! | lat            | latitude_index_in_debug_printouts                      | latitude index in debug printouts                                  | index      |    0 | integer   |           | in     | F        |
!! | u1             | x_wind                                                 | zonal wind                                                         | m s-1      |    2 | real      | kind_phys | in     | F        |
!! | v1             | y_wind                                                 | meridional wind                                                    | m s-1      |    2 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature                                        | mid-layer temperature                                              | K          |    2 | real      | kind_phys | in     | F        |
!! | q1             | water_vapor_specific_humidity                          | mid-layer specific humidity of water vapor                         | kg kg-1    |    2 | real      | kind_phys | in     | F        |
!! | deltim         | time_step_for_physics                                  | physics time step                                                  | s          |    0 | real      | kind_phys | in     | F        |
!! | pmid1          | air_pressure                                           | mid-layer pressure                                                 | Pa         |    2 | real      | kind_phys | in     | F        |
!! | pint1          | air_pressure_at_interface                              | interface pressure                                                 | Pa         |    2 | real      | kind_phys | in     | F        |
!! | dpmid1         | air_pressure_difference_between_midlayers              | difference between mid-layer pressures                             | Pa         |    2 | real      | kind_phys | in     | F        |
!! | qmax           | maximum_column_heating_rate                            | maximum heating rate in column                                     | K s-1      |    1 | real      | kind_phys | in     | F        |
!! | ktop           | vertical_index_at_cloud_top                            | vertical index at cloud top                                        | index      |    1 | integer   |           | in     | F        |
!! | kbot           | vertical_index_at_cloud_base                           | vertical index at cloud base                                       | index      |    1 | integer   |           | in     | F        |
!! | kcnv           | flag_deep_convection                                   | flag indicating whether convection occurs in column (0 or 1)       | flag       |    1 | integer   |           | in     | F        |
!! | cldf           | cloud_area_fraction                                    | fraction of grid box area in which updrafts occur                  | frac       |    1 | real      | kind_phys | in     | F        |
!! | grav           | gravitational_acceleration                             | gravitational acceleration                                         | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | cp             | specific_heat_of_dry_air_at_constant_pressure          | specific heat of dry air at constant pressure                      | J kg-1 K-1 |    0 | real      | kind_phys | in     | F        |
!! | rd             | gas_constant_dry_air                                   | ideal gas constant for dry air                                     | J kg-1 K-1 |    0 | real      | kind_phys | in     | F        |
!! | fv             | ratio_of_vapor_to_dry_air_gas_constants_minus_one      | rv/rd - 1 (rv = ideal gas constant for water vapor)                | none       |    0 | real      | kind_phys | in     | F        |
!! | pi             | pi                                                     | ratio of a circle's circumference to its diameter                  | radians    |    0 | real      | kind_phys | in     | F        |
!! | dlength        | characteristic_grid_length_scale                       | representative horizontal length scale of grid box                 | m          |    1 | real      | kind_phys | in     | F        |
!! | lprnt          | flag_print                                             | flag for debugging printouts                                       | flag       |    0 | logical   |           | in     | F        |
!! | ipr            | horizontal_index_of_printed_column                     | horizontal index of column used in debugging printouts             | index      |    0 | integer   |           | in     | F        |
!! | fhour          | forecast_time                                          | forecast hour                                                      | h          |    0 | real      | kind_phys | in     | F        |
!! | utgwc          | tendency_of_x_wind_due_to_convective_gravity_wave_drag | zonal wind tendency due to convective gravity wave drag            | m s-2      |    2 | real      | kind_phys | out    | F        |
!! | vtgwc          | tendency_of_y_wind_due_to_convective_gravity_wave_drag | meridional wind tendency due to convective gravity wave drag       | m s-2      |    2 | real      | kind_phys | out    | F        |
!! | tauctx         | instantaneous_x_stress_due_to_gravity_wave_drag        | zonal stress at cloud top due to convective gravity wave drag      | Pa         |    1 | real      | kind_phys | out    | F        |
!! | taucty         | instantaneous_y_stress_due_to_gravity_wave_drag        | meridional stress at cloud top due to convective gravity wave drag | Pa         |    1 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                           | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                              | flag       |    0 | integer   |           | out    | F        |
!!
!>\section gen_gwdc GFS Convective GWD Scheme General Algorithm
!! Parameterizing subgrid-scale convection-induced gravity wave
!! momentum flux for use in large-scale models inherently requires
!! some information from subgrid-scale cumulus parameterization.
!! The methodology for parameterizing the zonal momentum flux induced
!! by thermal forcing can be summarized as follows. From the cloud-base
!! to cloud-top height, the effect of the momentum flux
!! induced by subgrid-scale diabatic forcing is not considered because
!! subgrid-scale cumulus convection in large-scale models is only
!! activated in a conditionally unstable atmosphere. Below the cloud
!! base, the momentum flux is also not considered because of the wave
!! momentum cancellation.
!!
!! The formulation of the momentum flux at the cloud top that can be broken
!! in the upper atmosphere is similar to that of the surface drag of mountain
!! waves by using the nonlinearity factor of thermally induced waves that is
!! analogous to the inverse Frounde number in mountain waves. A vertical
!! profiles of redistributed momentum flux above the cloud top can be approximated
!! either by specifying a functional form or by using the wave saturation hypothesis
!! in terms of the Richardson number criterion. The formulation of the minimum
!! Richardson number including wave impact is similar to that in the mountain
!! wave case.
!!
!> \section al_gwdc GFS Convective GWD Scheme Detailed Algorithm
!> @{
      subroutine gwdc_run (im,ix,km,lat,u1,v1,t1,q1,deltim,             &
     &           pmid1,pint1,dpmid1,qmax,ktop,kbot,kcnv,cldf,           &
     &           grav,cp,rd,fv,pi,dlength,lprnt,ipr,fhour,              &
     &           utgwc,vtgwc,tauctx,taucty,errmsg,errflg)


!***********************************************************************
! aug   2005 Ake Johansson - ORIGINAL CODE FOR PARAMETERIZATION OF CONVECTIVELY FORCED
!                            GRAVITY WAVE DRAG FROM YONSEI UNIVERSITY, KOREA
!                            BASED ON THE THEORY GIVEN BY CHUN AND BAIK (JAS, 1998)
!                            MODIFIED FOR IMPLEMENTATION INTO THE GFS/CFSD BY
!        2013   S. Moorthi - Updated and optimized code for T1534 GFS implementation
! ??? ?? 2015   J. Alpert  - reducing the magnitude of tauctmax to fix blow up in L64 GFS
!               S. Kar & M. Young
! aug 15 2016 - S. Moorthi - Fix for exessive dissipation which led to blow up in
!                            128 level runs with NEMS/GSM
!***********************************************************************

      USE MACHINE , ONLY : kind_phys
      implicit none

! DH* TODO - add intent information for all variables
!---------------------------- Arguments --------------------------------
!
!  Input variables
!
!  u        : midpoint zonal wind
!  v        : midpoint meridional wind
!  t        : midpoint temperatures
!  pmid     : midpoint pressures
!  pint     : interface pressures
!  dpmid    : midpoint delta p ( pi(k)-pi(k-1) )
!  lat      : latitude index
!  qmax     : deep convective heating
!  kcldtop  : Vertical level index for cloud top    ( mid level )
!  kcldbot  : Vertical level index for cloud bottom ( mid level )
!  kcnv     : (0,1) dependent on whether convection occur or not
!
!  Output variables
!
!  utgwc    : zonal wind tendency
!  vtgwc    : meridional wind tendency
!
!-----------------------------------------------------------------------

      integer, intent(in) :: im, ix, km, lat, ipr
      integer, intent(in) :: ktop(im),kbot(im),kcnv(im)
      real(kind=kind_phys), intent(in) :: grav,cp,rd,fv,fhour,deltim,pi
      real(kind=kind_phys), dimension(im), intent(in) :: qmax
      real(kind=kind_phys), dimension(im), intent(out) :: tauctx,taucty
      real(kind=kind_phys), dimension(im), intent(in) :: cldf,dlength
      real(kind=kind_phys), dimension(ix,km), intent(in) :: u1,v1,t1,   &
     &                                                  q1,pmid1,dpmid1
      real(kind=kind_phys), dimension(ix,km), intent(out) :: utgwc,vtgwc
      real(kind=kind_phys), dimension(ix,km+1), intent(in) :: pint1
!
      logical, intent(in) :: lprnt
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!------------------------- Local workspace -----------------------------
!
!  i, k     : Loop index
!  kk       : Loop index
!  cldf     : Deep convective cloud fraction at the cloud top.
!  ugwdc    : Zonal wind after GWDC paramterization
!  vgwdc    : Meridional wind after GWDC parameterization
!  plnmid   : Log(pmid) ( mid level )
!  plnint   : Log(pint) ( interface level )
!  dpint    : Delta pmid ( interface level )
!  tauct    : Wave stress at the cloud top calculated using basic-wind
!             parallel to the wind vector at the cloud top ( mid level )
!  tauctx   : Wave stress at the cloud top projected in the east
!  taucty   : Wave stress at the cloud top projected in the north
!  qmax     : Maximum deep convective heating rate ( K s-1 ) in a
!             horizontal grid point calculated from cumulus para-
!             meterization. ( mid level )
!  wtgwc    : Wind tendency in direction to the wind vector at the cloud top level
!             due to convectively generated gravity waves ( mid level )
!  utgwcl   : Zonal wind tendency due to convectively generated
!             gravity waves ( mid level )
!  vtgwcl   : Meridional wind tendency due to convectively generated
!             gravity waves ( mid level )
!  taugwci  : Profile of wave stress calculated using basic-wind
!             parallel to the wind vector at the cloud top
!  taugwcxi : Profile of zonal component of gravity wave stress
!  taugwcyi : Profile of meridional component of gravity wave stress
!
!  taugwci, taugwcxi, and taugwcyi are defined at the interface level
!
!  bruni    : Brunt-Vaisala frequency ( interface level )
!  brunm    : Brunt-Vaisala frequency ( mid level )
!  rhoi     : Air density ( interface level )
!  rhom     : Air density ( mid level )
!  ti       : Temperature ( interface level )
!  basicum  : Basic-wind profile. Basic-wind is parallel to the wind
!             vector at the cloud top level. (mid level)
!  basicui  : Basic-wind profile. Basic-wind is parallel to the wind
!             vector at the cloud top level. ( interface level )
!  riloc    : Local Richardson number ( interface level )
!  rimin    : Minimum Richardson number including both the basic-state
!             and gravity wave effects ( interface level )
!  gwdcloc  : Horizontal location where the GWDC scheme is activated.
!  break    : Horizontal location where wave breaking is occurred.
!  critic   : Horizontal location where critical level filtering is
!             occurred.
!  dogwdc   : Logical flag whether the GWDC parameterization is
!             calculated at a grid point or not.
!
!  dogwdc is used in order to lessen CPU time for GWDC calculation.
!
!-----------------------------------------------------------------------

      integer i,ii,k,k1,kk,kb,ilev,npt,kcb,kcldm,npr
      integer, dimension(im) :: ipt

      real(kind=kind_phys) tem, tem1,  tem2, qtem, wtgwc, tauct,        &
     &                     windcltop,  shear, nonlinct, nonlin, nonlins,&
     &                     n2,   dtdp,  crit1, crit2,     p1, p2,       &
!    &                     n2,   dtdp,  crit1, crit2, pi, p1, p2,
     &                     gsqr,  onebg
!    &                     taus, n2,   dtdp,  crit1, crit2, pi, p1, p2

      integer,              allocatable :: kcldtop(:),kcldbot(:)
      logical,              allocatable :: do_gwc(:)
      real(kind=kind_phys), allocatable :: tauctxl(:), tauctyl(:),
     &                                     gwdcloc(:), break(:),
!    &                                     critic(:),
!    &                                     critic(:),  angle(:),
     &                                     cosphi(:),  sinphi(:),
     &                                     xstress(:), ystress(:),
     &                                     ucltop(:),  vcltop(:),
     &                                     wrk(:),     dtfac(:),
     &                                     dlen(:),       gqmcldlen(:)
!     real(kind=kind_phys), allocatable :: plnint(:,:),   dpint(:,:),
!    &                                     taugwci(:,:),  taugwcxi(:,:),
!    &                                     taugwcyi(:,:), bruni(:,:),
!    &                                     taugwcyi(:,:), bruni(:,:),
      real(kind=kind_phys), allocatable :: plnint(:,:),   velco(:,:),
     &                                     taugwci(:,:),  bruni(:,:),
     &                                     rhoi(:,:),     basicui(:,:),
     &                                     ti(:,:),       riloc(:,:),
     &                                     rimin(:,:),    pint(:,:)
!     real(kind=kind_phys), allocatable :: ugwdc(:,:),    vgwdc(:,:),
      real(kind=kind_phys), allocatable ::
!    &                                     plnmid(:,:),   wtgwc(:,:),
     &                                     plnmid(:,:),   taugw(:,:),
     &                                     utgwcl(:,:),   vtgwcl(:,:),
     &                                     basicum(:,:),  u(:,:),v(:,:),
     &                                     t(:,:),        spfh(:,:),
     &                                     pmid(:,:),     dpmid(:,:),
!    &                                     pmid(:,:),     cumchr(:,:),
     &                                     brunm(:,:),    rhom(:,:)

!-----------------------------------------------------------------------
!
!  ucltop    : Zonal wind at the cloud top ( mid level )
!  vcltop    : Meridional wind at the cloud top ( mid level )
!  windcltop : Wind speed at the cloud top ( mid level )
!  shear     : Vertical shear of basic wind
!  cosphi    : Cosine of angle of wind vector at the cloud top
!  sinphi    : Sine   of angle of wind vector at the cloud top
!  c1        : Tunable parameter
!  c2        : Tunable parameter
!  dlength   : Grid spacing in the direction of basic wind at the cloud top
!  nonlinct  : Nonlinear parameter at the cloud top
!  nonlin    : Nonlinear parameter above the cloud top
!  nonlins   : Saturation nonlinear parameter
!  taus      : Saturation gravity wave drag == taugwci(i,k)
!  n2        : Square of Brunt-Vaisala frequency
!  dtdp      : dT/dp
!  xstress   : Vertically integrated zonal momentum change due to GWDC
!  ystress   : Vertically integrated meridional momentum change due to GWDC
!  crit1     : Variable 1 for checking critical level
!  crit2     : Variable 2 for checking critical level
!
!-----------------------------------------------------------------------

      real(kind=kind_phys), parameter ::
     &                      c1=1.41,          c2=-0.38,     ricrit=0.25
     &,                     n2min=1.e-32,     zero=0.0,     one=1.0
     &,                     taumin=1.0e-20,   tauctmax=-20.
!    &,                     taumin=1.0e-20,   tauctmax=-5.
     &,                     qmin=1.0e-10,     shmin=1.0e-20
     &,                     rimax=1.0e+20,    rimaxm=0.99e+20
     &,                     rimaxp=1.01e+20,  rilarge=0.9e+20
     &,                     riminx=-1.0e+20,  riminm=-1.01e+20
     &,                     riminp=-0.99e+20, rismall=-0.9e+20

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!
!> - Determine if deep convection occurs and activate convection-induced GWD scheme.
      npt = 0
      do i = 1,im
        ipt(i) = 0
        if (kcnv(i) /= 0 .and. qmax(i) > zero) then
          npt      = npt + 1
          ipt(npt) = i
        endif
      enddo
      do k=1,km
        do i=1,im
          utgwc(i,k) = 0.0
          vtgwc(i,k) = 0.0
!         brunm(i,k) = 0.0
!         rhom(i,k)  = 0.0
        enddo
      enddo
      do i=1,im
        tauctx(i) = 0.0
        taucty(i) = 0.0
      enddo
      if (npt == 0) return      ! No gwdc calculation done!

!***********************************************************************
!
!  Begin GWDC
!
!***********************************************************************

!-----------------------------------------------------------------------
!        Write out incoming variables
!-----------------------------------------------------------------------

!     fhourpr = zero
!     if (lprnt) then
!       if (fhour >= fhourpr) then
!         print *,' '
!         write(*,*) 'Inside GWDC raw input start print at fhour = ',
!    &               fhour
!         write(*,*) 'IX  IM  KM  ',ix,im,km
!         write(*,*) 'KBOT KTOP QMAX DLENGTH kcnv  ',
!    +     kbot(ipr),ktop(ipr),qmax(ipr),dlength(ipr),kcnv(ipr)
!         write(*,*) 'grav  cp  rd  ',grav,cp,rd

!-------- Pressure levels ----------
!         write(*,9100)
!         ilev=km+1
!         write(*,9110) ilev,(10.*pint1(ipr,ilev))
!         do ilev=km,1,-1
!           write(*,9120) ilev,(10.*pmid1(ipr,ilev)),
!    &                         (10.*dpmid1(ipr,ilev))
!           write(*,9110) ilev,(10.*pint1(ipr,ilev))
!         enddo

!-------- U1 V1 T1 ----------
!         write(*,9130)
!         do ilev=km,1,-1
!           write(*,9140) ilev,U1(ipr,ilev),V1(ipr,ilev),T1(ipr,ilev)
!         enddo

!         print *,' '
!         print *,' Inside GWDC raw input end print'
!       endif
!     endif

!9100 format(//,14x,'PRESSURE LEVELS',//,
!    +' ILEV',6x,'PINT1',7x,'PMID1',6x,'DPMID1',/)
!9110 format(i4,2x,f10.3)
!9120 format(i4,12x,2(2x,f10.3))
!9130 format(//,' ILEV',7x,'U1',10x,'V1',10x,'T1',/)
!9140 format(i4,3(2x,f10.3))

!     Allocate local arrays

      allocate (kcldtop(npt), kcldbot(npt), do_gwc(npt))
      allocate (tauctxl(npt), tauctyl(npt), dtfac(npt),
     &          gwdcloc(npt), break(npt),                cosphi(npt),
!    &          gwdcloc(npt), break(npt), critic(npt),   cosphi(npt),
     &          sinphi(npt),  xstress(npt),  ystress(npt), wrk(npt),
     &          ucltop(npt),  vcltop(npt),dlen(npt),     gqmcldlen(npt))

!     allocate (plnint(npt,2:km+1), dpint(npt,km+1),
!    &          taugwci(npt,km+1),  taugwcxi(npt,km+1),
!    &          taugwcyi(npt,km+1), bruni(npt,km+1),
      allocate (plnint(npt,2:km+1),
     &          taugwci(npt,km+1),  bruni(npt,km+1),
     &          rhoi(npt,km+1),     basicui(npt,km+1),
     &          ti(npt,km+1),       riloc(npt,km+1),
     &          rimin(npt,km+1),    pint(npt,km+1))

!     allocate (ugwdc(npt,km),   vgwdc(npt,km),
      allocate
!    &         (plnmid(npt,km),  wtgwc(npt,km),
     &         (plnmid(npt,km),  velco(npt,km),
     &          utgwcl(npt,km),  vtgwcl(npt,km),
     &          basicum(npt,km), u(npt,km),    v(npt,km),
     &          t(npt,km),       spfh(npt,km), pmid(npt,km),
     &          dpmid(npt,km),   taugw(npt,km),
!    &          dpmid(npt,km),   cumchr(npt,km),
     &          brunm(npt,km),   rhom(npt,km))

!-----------------------------------------------------------------------
!> - Create local arrays with reversed vertical indices
!!   and Initialize local variables.
!-----------------------------------------------------------------------
      gsqr  = grav * grav
      onebg = one / grav

      if (lprnt) then
        npr = 1
        do i=1,npt
          if (ipr == ipt(i))then
            npr = i
            exit
          endif
        enddo
      endif

      do k=1,km
        k1 = km - k + 1
        do i=1,npt
          ii = ipt(i)
          u(i,k)        = u1(ii,k1)
          v(i,k)        = v1(ii,k1)
          t(i,k)        = t1(ii,k1)
          spfh(i,k)     = max(q1(ii,k1),qmin)
          pmid(i,k)     = pmid1(ii,k1)
          dpmid(i,k)    = dpmid1(ii,k1) * onebg
!         cumchr(i,k)   = cumchr1(ii,k1)

          rhom(i,k)     = pmid(i,k) / (rd*t(i,k)*(1.0+fv*spfh(i,k)))
          plnmid(i,k)   = log(pmid(i,k))
          utgwcl(i,k)   = zero
          vtgwcl(i,k)   = zero
!         ugwdc(i,k)    = zero
!         vgwdc(i,k)    = zero
          brunm(i,k)    = zero
          basicum(i,k)  = zero
        enddo
      enddo

       do k=1,km+1
        k1 = km - k + 2
        do i=1,npt
          ii = ipt(i)
          pint(i,k)     = pint1(ii,k1)
          taugwci(i,k)  = zero
          bruni(i,k)    = zero
          rhoi(i,k)     = zero
          ti(i,k)       = zero
          basicui(i,k)  = zero
          riloc(i,k)    = zero
          rimin(i,k)    = zero
        enddo
      enddo
      do k=2,km+1
        do i=1,npt
          plnint(i,k)   = log(pint(i,k))
        enddo
      enddo

      do i = 1, npt
        ii = ipt(i)
        kcldtop(i)   = km - ktop(ii) + 1
        kcldbot(i)   = km - kbot(ii) + 1
        dlen(i)      = dlength(ii)
!                                    (g*qmax(ii)*cldf(ii)*dlength(ii))
        gqmcldlen(i) = grav*qmax(ii)*cldf(ii)*dlen(i)
      enddo
!     if (lprnt) write(7000,*)' ktop=',ktop(ipr),' kbot=',kbot(ipr)
!    &,' kcldtop=',kcldtop(npr),' kcldbot=',kcldbot(npr),
!    &' dlength=',dlength(ipr),' qmax=',qmax(ipr),' cldf=',cldf(ipr)

!     if (lprnt) then
!       if (fhour.ge.fhourpr) then
!         write(*,9200)
!         do i=1,im
!           write(*,9201) kcnv(i),kcldbot(i),kcldtop(i)
!         enddo
!       endif
!     endif

!9200 format(//,'  Inside GWDC local variables start print',//,
!    +2x,'kcnv',2x,'KCLDBOT',2x,'KCLDTOP',//)
!9201 format(i4,2x,i5,4x,i5)

!***********************************************************************

!     pi     = 2.*asin(1.)

!-----------------------------------------------------------------------
!
!                              PRESSURE VARIABLES
!
!  Interface 1 ======== pint(1)           *********
!  Mid-Level 1 --------          pmid(1)            dpmid(1)
!            2 ======== pint(2)           dpint(2)
!            2 --------          pmid(2)            dpmid(2)
!            3 ======== pint(3)           dpint(3)
!            3 --------          pmid(3)            dpmid(3)
!            4 ======== pint(4)           dpint(4)
!            4 --------          pmid(4)            dpmid(4)
!              ........
!           17 ======== pint(17)          dpint(17)
!           17 --------          pmid(17)           dpmid(17)
!           18 ======== pint(18)          dpint(18)
!           18 --------          pmid(18)           dpmid(18)
!           19 ======== pint(19)          *********
!
!-----------------------------------------------------------------------

      do i = 1, npt
        tauctxl(i)    = zero
        tauctyl(i)    = zero

!-----------------------------------------------------------------------
!                              THERMAL VARIABLES
!
!  Interface 1 ========       TI(1)           RHOI(1)            BRUNI(1)
!            1 -------- T(1)         RHOM(1)           BRUNM(1)
!            2 ========       TI(2)           RHOI(2)            BRUNI(2)
!            2 -------- T(2)         RHOM(2)           BRUNM(2)
!            3 ========       TI(3)           RHOI(3)            BRUNI(3)
!            3 -------- T(3)         RHOM(3)           BRUNM(3)
!            4 ========       TI(4)           RHOI(4)            BRUNI(4)
!            4 -------- T(4)         RHOM(4)           BRUNM(4)
!              ........
!           17 ========
!           17 -------- T(17)        RHOM(17)          BRUNM(17)
!           18 ========       TI(18)          RHOI(18)           BRUNI(18)
!           18 -------- T(18)        RHOM(18)          BRUNM(18)
!           19 ========       TI(19)          RHOI(19)           BRUNI(19)
!

!
!>  - The top interface temperature, density, and Brunt-Vaisala
!!    frequencies (\f$N\f$) are calculated assuming an isothermal
!!    atmosphere above the top mid level.

        ti(i,1)    = t(i,1)
        rhoi(i,1)  = pint(i,1)/(rd*ti(i,1))
        bruni(i,1) = sqrt ( gsqr / (cp*ti(i,1)) )
!
!>  - The bottom interface temperature, density, and Brunt-Vaisala
!!    frequencies (\f$N\f$) are calculated assuming an isothermal
!!    atmosphere below the bottom mid level.

        ti(i,km+1)    = t(i,km)
        rhoi(i,km+1)  = pint(i,km+1)/(rd*ti(i,km+1)*(1.0+fv*spfh(i,km)))
        bruni(i,km+1) = sqrt ( gsqr / (cp*ti(i,km+1)) )
      enddo

!-----------------------------------------------------------------------
!
!>  - The interface level temperature, density, and Brunt-Vaisala
!!    frequencies (\f$N\f$) are calculated based on linear interpolation
!!    of temperature in ln(P).
!
!-----------------------------------------------------------------------

      do k = 2, km
        do i = 1, npt
          tem1 = (plnmid(i,k)-plnint(i,k)) / (plnmid(i,k)-plnmid(i,k-1))
          tem2 = one - tem1
          ti(i,k)    = t(i,k-1)    * tem1 + t(i,k)    * tem2
          qtem       = spfh(i,k-1) * tem1 + spfh(i,k) * tem2
          rhoi(i,k)  = pint(i,k) / ( rd * ti(i,k)*(1.0+fv*qtem) )
          dtdp       = (t(i,k)-t(i,k-1)) / (pmid(i,k)-pmid(i,k-1))
          n2         = gsqr / ti(i,k) * ( 1./cp - rhoi(i,k)*dtdp )
          bruni(i,k) = sqrt (max (n2min, n2))
        enddo
      enddo

      deallocate (spfh)
!-----------------------------------------------------------------------
!
!>  - The mid-level Brunt-Vaisala frequencies (\f$N\f$) are calculated
!!    based on interpolated interface temperatures.
!-----------------------------------------------------------------------

      do k = 1, km
        do i = 1, npt
          dtdp       = (ti(i,k+1)-ti(i,k)) / (pint(i,k+1)-pint(i,k))
          n2         = gsqr / t(i,k) * ( 1./cp - rhom(i,k)*dtdp )
          brunm(i,k) = sqrt (max (n2min, n2))
        enddo
      enddo

!-----------------------------------------------------------------------
!        PRINTOUT
!-----------------------------------------------------------------------

!     if (lprnt) then
!       if (fhour.ge.fhourpr) then

!-------- Pressure levels ----------
!         write(*,9101)
!         do ilev=1,km
!           write(*,9111) ilev,(0.01*pint(ipr,ilev)),
!    &                         (0.01*dpint(ipr,ilev)),plnint(ipr,ilev)
!           write(*,9121) ilev,(0.01*pmid(ipr,ilev)),
!    &                         (0.01*dpmid(ipr,ilev)),plnmid(ipr,ilev)
!         enddo
!         ilev=km+1
!         write(*,9111) ilev,(0.01*pint(ipr,ilev)),
!    &                       (0.01*dpint(ipr,ilev)),plnint(ipr,ilev)

!                2
!-------- U V T N  ----------
!         write(*,9102)
!         do ilev=1,km
!           write(*,9112) ilev,ti(ipr,ilev),(100.*bruni(ipr,ilev))
!           write(*,9122) ilev,u(ipr,ilev),v(ipr,ilev),
!    +                    t(ipr,ilev),(100.*brunm(ipr,ilev))
!         enddo
!         ilev=km+1
!         write(*,9112) ilev,ti(ipr,ilev),(100.*bruni(ipr,ilev))

!       endif
!     endif

!9101 format(//,14x,'PRESSURE LEVELS',//,
!    +' ILEV',4x,'PINT',4x,'PMID',4x,'DPINT',3x,'DPMID',5x,'LNP',/)
!9111 format(i4,1x,f8.2,9x,f8.2,9x,f8.2)
!9121 format(i4,9x,f8.2,9x,f8.2,1x,f8.2)
!9102 format(//' ILEV',5x,'U',7x,'V',5x,'TI',7x,'T',
!    +5x,'BRUNI',3x,'BRUNM',//)
!9112 format(i4,16x,f8.2,8x,f8.3)
!9122 format(i4,2f8.2,8x,f8.2,8x,f8.3)


!***********************************************************************
!
!        Big loop over grid points                    ONLY done if kcnv=1
!
!***********************************************************************

      kcldm = 1
      do i = 1, npt
        kk        = kcldtop(i)
        kb        = kcldbot(i)
        kcldm     = max(kcldm,kk)

!-----------------------------------------------------------------------
!
!> - Calculate the cloud top wind components, speed and direction.
! Here, ucltop, vcltop, and windcltop are wind components and
!  wind speed at mid-level cloud top index
!
!-----------------------------------------------------------------------

        ucltop(i) = u(i,kk)
        vcltop(i) = v(i,kk)
!       windcltop = sqrt( ucltop(i)*ucltop(i) + vcltop(i)*vcltop(i) )
        windcltop = 1.0 / sqrt( ucltop(i)*ucltop(i)
     &                        + vcltop(i)*vcltop(i) )
        cosphi(i) = ucltop(i)*windcltop
        sinphi(i) = vcltop(i)*windcltop
!       angle(i)  = acos(cosphi)*180./pi
      enddo

!-----------------------------------------------------------------------
!
!> - Calculate the basic state wind profiles projected in the direction of the
!!  cloud top wind at mid level and interface level. 
! \n  U : Basic-wind speed profile. Basic-wind is parallel to the wind
!             vector at the cloud top level. (mid level)
! \n  UI: Basic-wind speed profile. Basic-wind is parallel to the wind
!             vector at the cloud top level. ( interface level )
!  Input u(i,k) and v(i,k) is defined at mid level
!
!-----------------------------------------------------------------------

      do k=1,km
        do i=1,npt
          basicum(i,k) = u(i,k)*cosphi(i) + v(i,k)*sinphi(i)
        enddo
      enddo

!-----------------------------------------------------------------------
!
!  Basic state wind at interface level is also calculated
!  based on linear interpolation in ln(Pressure)
!
!  In the top and bottom boundaries, basic-state wind at interface level
!  is assumed to be vertically uniform.
!
!-----------------------------------------------------------------------

      do i=1,npt
        basicui(i,1)    = basicum(i,1)
        basicui(i,km+1) = basicum(i,km)
      enddo
      do k=2,km
        do i=1,npt
          tem1 = (plnmid(i,k)-plnint(i,k)) / (plnmid(i,k)-plnmid(i,k-1))
          tem2 = one - tem1
          basicui(i,k) = basicum(i,k)*tem2 + basicum(i,k-1)*tem1
        enddo
      enddo

!-----------------------------------------------------------------------
!
!> - Calculate the local Richardson number:
!! \f[
!!    Ri=N^2/\eta^2
!! \f]
!! where \f$\eta\f$ is the vertical shear (\f$dU/dz\f$).

!  basicum   : U at mid level
!  basicui   : UI at interface level
!
!  Interface 1 ========       UI(1)            rhoi(1)  bruni(1)  riloc(1)
!  Mid-level 1 -------- U(1)
!            2 ========       UI(2)  dpint(2)  rhoi(2)  bruni(2)  riloc(2)
!            2 -------- U(2)
!            3 ========       UI(3)  dpint(3)  rhoi(3)  bruni(3)  riloc(3)
!            3 -------- U(3)
!            4 ========       UI(4)  dpint(4)  rhoi(4)  bruni(4)  riloc(4)
!            4 -------- U(4)
!              ........
!           17 ========       UI(17) dpint(17) rhoi(17) bruni(17) riloc(17)
!           17 -------- U(17)
!           18 ========       UI(18) dpint(18) rhoi(18) bruni(18) riloc(18)
!           18 -------- U(18)
!           19 ========       UI(19)           rhoi(19) bruni(19) riloc(19)
!
!-----------------------------------------------------------------------

      do k=2,km
        do i=1,npt
          shear = grav*rhoi(i,k) * (basicum(i,k) - basicum(i,k-1))
     &                           / (pmid(i,k) - pmid(i,k-1))
          if ( abs(shear) < shmin ) then
            riloc(i,k) = rimax
          else
            tem = bruni(i,k) / shear
            riloc(i,k)  = tem * tem
            if (riloc(i,k) >= rimax ) riloc(i,k) = rilarge
          end if
        enddo
      enddo

      do i=1,npt
        riloc(i,1)    = riloc(i,2)
        riloc(i,km+1) = riloc(i,km)
      enddo

!     if (lprnt.and.(i.eq.ipr)) then
!       if (fhour.ge.fhourpr) then
!         write(*,9104) ucltop,vcltop,windcltop,angle,kk
!         do ilev=1,km
!           write(*,9114) ilev,basicui(ipr,ilev),dpint(ipr,ilev),
!    +      rhoi(ipr,ilev),(100.*bruni(ipr,ilev)),riloc(ilev)
!           write(*,9124) ilev,(basicum(ipr,ilev))
!         enddo
!         ilev=km+1
!         write(*,9114) ilev,basicui(ipr,ilev),dpint(ipr,ilev),
!    +      rhoi(ipr,ilev),(100.*bruni(ipr,ilev)),riloc(ilev)
!       endif
!     endif

!9104 format(//,'WIND VECTOR AT CLOUDTOP = (',f6.2,' , ',f6.2,' ) = ',
!    +f6.2,' IN DIRECTION ',f6.2,4x,'KK = ',i2,//,
!    +' ILEV',2x,'BASICUM',2x,'BASICUI',4x,'DPINT',6x,'RHOI',5x,
!    +'BRUNI',6x,'RI',/)
!9114 format(i4,10x,f8.2,4(2x,f8.2))
!9124 format(i4,1x,f8.2)

!-----------------------------------------------------------------------
!
!> - Calculate the gravity wave stress at the interface level cloud top.
!
!  kcldtopi  : The interface level cloud top index
!  kcldtop   : The midlevel cloud top index
!  kcldbot   : The midlevel cloud bottom index
!
!  A : Find deep convective heating rate maximum
!
!      If kcldtop(i) is less than kcldbot(i) in a horizontal grid point,
!      it can be thought that there is deep convective cloud. However,
!      deep convective heating between kcldbot and kcldtop is sometimes
!      zero in spite of kcldtop less than kcldbot. In this case,
!      maximum deep convective heating is assumed to be 1.e-30.
!
!  B : kk is the vertical index for interface level cloud top
!
!  C : Total convective fractional cover (cldf) is used as the
!      convective cloud cover for GWDC calculation instead of
!      convective cloud cover in each layer (concld).
!                       a1 = cldf*dlength
!      You can see the difference between cldf(i) and concld(i)
!      in (4.a.2) in Description of the NCAR Community Climate
!      Model (CCM3).
!      In NCAR CCM3, cloud fractional cover in each layer in a deep
!      cumulus convection is determined assuming total convective
!      cloud cover is randomly overlapped in each layer in the
!      cumulus convection.
!
!  D : Wave stress at cloud top is calculated when the atmosphere
!      is dynamically stable at the cloud top
!
!  E : Cloud top wave stress and nonlinear parameter are calculated
!      using density, temperature, and wind that are defined at mid
!      level just below the interface level in which cloud top wave
!      stress is defined.
!      Nonlinct is defined at the interface level.
!
!  F : If the atmosphere is dynamically unstable at the cloud top,
!      GWDC calculation in current horizontal grid is skipped.
!
!  G : If mean wind at the cloud top is less than zero, GWDC
!      calculation in current horizontal grid is skipped.
!
!-----------------------------------------------------------------------
!D
!>  - Wave stress at cloud top is calculated when the atmosphere
!!    is dynamically stable at the cloud top.
      do i=1,npt
        kk = kcldtop(i)
        if ( abs(basicui(i,kk)) > zero .and. riloc(i,kk) > ricrit) then
!E
!>   - The cloud top wave stress and nonlinear parameter are calculated
!!      using density, temperature, and wind that are defined at mid
!!      level just below the interface level in which cloud top wave
!!      stress is defined.
!! The parameter \f$\mu\f$ is the nonlinearity factor of thermally
!! induced internal gravity waves defined by eq.(17) in \cite chun_and_baik_1998 :
!! \f[
!!  \mu=\frac{gQ_{0}a_{1}}{c_{p}T_{0}NU^{2}}
!! \f]
!! where \f$Q_{0}\f$ is the maximum deep convective heating rate in a
!! horizontal grid point calculated from cumulus parameterization.
!! \f$a_{1}\f$ is the half-width of
!! the forcing function.\f$g\f$ is gravity. \f$c_{p}\f$ is specific
!! heat at constant pressure. \f$T_{0}\f$ is the layer mean
!! temperature (T1). As eqs.(18) and (19) \cite chun_and_baik_1998,
!! the zonal momentum flux is given by
!! \f[
!! \tau_{x}=-[\rho U^{3}/(N\triangle x)]G(\mu)
!! \f]
!! where
!! \f[
!! G(\mu)=c_{1}c_2^2 \mu^{2}
!! \f]
!! wher \f$\rho\f$ is the local density.
!! The tunable parameter \f$c_1\f$ is related to the horizontal
!! structure of thermal forcing. The tunable parameter \f$c_2\f$ is
!! related to the basic-state wind and stability and the bottom and
!! top heights of thermal forcing. If the atmosphere is dynamically
!! unstable at the cloud top, the convective GWD calculation is
!! skipped at that grid point.

          tem       = basicum(i,kk)
          tem1      = tem * tem
          nonlinct  = gqmcldlen(i) / (bruni(i,kk)*t(i,kk)*tem1)    ! Mu
          tem2      = c2*nonlinct
!                                  RhoU^3c1(c2mu)^2/Ndx
          tauct     = - rhom(i,kk) * tem * tem1 * c1 * tem2 * tem2
     &              /  (bruni(i,kk)*dlen(i))
!>   - The stress is capped at \f$\tau_{x} =  - 5n/m^2\f$
!!      in order to prevent numerical instability.

          tauct         = max(tauctmax, tauct)
          tauctxl(i)    = tauct * cosphi(i)           ! X stress at cloud top
          tauctyl(i)    = tauct * sinphi(i)           ! Y stress at cloud top
          taugwci(i,kk) = tauct                                    !  *1
          do_gwc(i)     = .true.
        else
!F
!>  - If the atmosphere is dynamically unstable at the cloud top,
!!  GWDC calculation in current horizontal grid is skipped.
!>  - If mean wind at the cloud top is less than zero, GWDC
!!      calculation in current horizontal grid is skipped.


          tauctxl(i) = zero
          tauctyl(i) = zero
          do_gwc(i) = .false.
        end if
!H
      enddo

!       if (lprnt.and.(i.eq.ipr)) then
!         if (fhour.ge.fhourpr) then
!            write(*,9210) tauctx(ipr),taucty(ipr),tauct(ipr),angle,kk
!         endif
!       endif

!9210 format(/,5x,'STRESS VECTOR = ( ',f8.3,' , ',f8.3,' ) = ',f8.3,
!    +' IN DIRECTION ',f6.2,4x,'KK = ',i2,/)

!-----------------------------------------------------------------------
!
!  At this point, mean wind at the cloud top is larger than zero and
!  local RI at the cloud top is larger than ricrit (=0.25)
!
!  Calculate minimum of Richardson number including both basic-state
!  condition and wave effects.
!
!          g*Q_0*alpha*dx                  RI_loc*(1 - mu*|c2|)
!  mu  =  ----------------  RI_min =  -----------------------------
!           c_p*N*T*U^2                (1 + mu*RI_loc^(0.5)*|c2|)^2
!
!  Minimum RI is calculated for the following two cases
!
!  (1)   RIloc < 1.e+20
!  (2)   Riloc = 1.e+20  ----> Vertically uniform basic-state wind
!
!  RIloc cannot be smaller than zero because N^2 becomes 1.E-32 in the
!  case of N^2 < 0.. Thus the sign of RINUM is determined by
!  1 - nonlin*|c2|.
!
!-----------------------------------------------------------------------
!> - Calculate the minimum Richardson number including both the
!! basic-state condition and wave effects.
!!\f[
!! Ri_{min}\approx\frac{Ri(1-\mu|c_{2}|)}{(1+\mu Ri^{1/2}|c_{2}|)^{2}}
!!\f]

      do k=kcldm,1,-1

        do i=1,npt
          if (do_gwc(i)) then
            kk = kcldtop(i)
            if (k > kk) cycle
            if ( k /= 1 ) then
              tem1  = (u(i,k)+u(i,k-1))*0.5
              tem2  = (v(i,k)+v(i,k-1))*0.5
              crit1 = ucltop(i)*tem1
              crit2 = vcltop(i)*tem2
              velco(i,k) = tem1 * cosphi(i) + tem2 * sinphi(i)
            else
              crit1 = ucltop(i)*u(i,1)
              crit2 = vcltop(i)*v(i,1)
              velco(i,1) = u(i,1) * cosphi(i) + v(i,1) * sinphi(i)
            end if
!     if (lprnt .and. i == npr) write(7000,*)' k=',k,' crit1=',
!    &crit1,' crit2=',crit2,' basicui=',basicui(i,k)

            if ( abs(basicui(i,k)) > zero .and. crit1 > zero
     &                                    .and. crit2 > zero ) then
              tem = basicui(i,k) * basicui(i,k)
              nonlin   = gqmcldlen(i) / (bruni(i,k)*ti(i,k)*tem)
              tem  = nonlin*abs(c2)
              if ( riloc(i,k)  <  rimaxm ) then
                tem1 = 1 + tem*sqrt(riloc(i,k))
                rimin(i,k) = riloc(i,k) * (1-tem) / (tem1*tem1)
              else if((riloc(i,k) > rimaxm) .and.
     &                (riloc(i,k) < rimaxp)) then
                rimin(i,k) = ( 1 - tem) / (tem*tem)
              end if
              if ( rimin(i,k) <= riminx ) then
                rimin(i,k) = rismall
              end if
            else
              rimin(i,k) = riminx
            end if
!     if (lprnt .and. i == npr) write(7000,*)' rimin=',rimin(i,k)

!-----------------------------------------------------------------------
!
! If the minimum \f$R_{i}\f$ at interface cloud top is less than or equal to 1/4,
!  the convective GWD calculation is skipped at that grid point.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!> - Calculate the gravity wave stress profile using the wave
!!  saturation hypothesis of Lindzen (1981) \cite lindzen_1981.
!
!  Assuming kcldtop(i)=10 and kcldbot=16,
!
!                             TAUGWCI  RIloc  RImin   UTGWC
!
!  Interface 1 ========       - 0.001         -1.e20
!            1 --------                               0.000
!            2 ========       - 0.001         -1.e20
!            2 --------                               0.000
!            3 ========       - 0.001         -1.e20
!            3 --------                               -.xxx
!            4 ========       - 0.001  2.600  2.000
!            4 --------                               0.000
!            5 ========       - 0.001  2.500  2.000
!            5 --------                               0.000
!            6 ========       - 0.001  1.500  0.110
!            6 --------                               +.xxx
!            7 ========       - 0.005  2.000  3.000
!            7 --------                               0.000
!            8 ========       - 0.005  1.000  0.222
!            8 --------                               +.xxx
!            9 ========       - 0.010  1.000  2.000
!            9 --------                               0.000
! kcldtopi  10 ========  $$$  - 0.010
! kcldtop   10 --------  $$$                          yyyyy
!           11 ========  $$$  0
!           11 --------  $$$
!           12 ========  $$$  0
!           12 --------  $$$
!           13 ========  $$$  0
!           13 --------  $$$
!           14 ========  $$$  0
!           14 --------  $$$
!           15 ========  $$$  0
!           15 --------  $$$
!           16 ========  $$$  0
! kcldbot   16 --------  $$$
!           17 ========       0
!           17 --------
!           18 ========       0
!           18 --------
!           19 ========       0
!
!-----------------------------------------------------------------------
!
!   Even though the cloud top level obtained in deep convective para-
!   meterization is defined in mid-level, the cloud top level for
!   the GWDC calculation is assumed to be the interface level just
!   above the mid-level cloud top vertical level index.
!
!-----------------------------------------------------------------------

!>  - When \f$Ri_{min}\f$ is set to 1/4 based on Lindzen's (1981)
!! \cite lindzen_1981 saturation hypothesis, the nonlinearity factor
!! for wave saturation can be derived by
!! \f[
!! \mu_{s}=\frac{1}{|c_{2}|}[2\sqrt{2+\frac{1}{\sqrt{Ri}}}-(2+\frac{1}{\sqrt{Ri}})]
!! \f]
!! Then the saturation zonal momentum flux is given by
!! \f[
!! \tau_{s}=-[\rho U^{3}/(N\triangle x)]c_{1}c_2^2\mu_s^2
!! \f]

            if (k < kk .and. k > 1) then
              if ( abs(taugwci(i,k+1)) > taumin ) then                  ! TAUGWCI
                if ( riloc(i,k) > ricrit ) then                         ! RIloc
                  if ( rimin(i,k) > ricrit ) then                       ! RImin
                    taugwci(i,k) = taugwci(i,k+1)
                  elseif (rimin(i,k) > riminp) then
                    tem = 2.0 + 1.0 / sqrt(riloc(i,k))
                    nonlins = (1.0/abs(c2)) * (2.*sqrt(tem) - tem)
                    tem1 = basicui(i,k)
                    tem2 = c2*nonlins*tem1
                    taugwci(i,k) = - rhoi(i,k) * c1 * tem1 * tem2 * tem2
     &                           /  (bruni(i,k)*dlen(i))
                  elseif (rimin(i,k) > riminm) then
                    taugwci(i,k) = zero
!                   taugwci(i,k) = taugwci(i,k+1)
                  end if                                              ! RImin
                else

!>  - If the minimum \f$R_{i}\f$ at interface cloud top is less than
!! or equal to 1/4, the convective GWD calculation is skipped at that
!! grid point.

                  taugwci(i,k) = zero
                end if                                                ! RIloc
              else
                 taugwci(i,k) = zero
              end if                                                  ! TAUGWCI

              if ( (basicum(i,k+1)*basicum(i,k) ) < 0. ) then
                 taugwci(i,k+1) = zero
                 taugwci(i,k)   = zero
              endif

              if (abs(taugwci(i,k)) > abs(taugwci(i,k+1))) then
                taugwci(i,k) = taugwci(i,k+1)
              end if

            elseif (k == 1) then

!>  - As an upper boundary condition, upward propagation of gravity
!!    wave energy is permitted.

              taugwci(i,1) = taugwci(i,2)
            endif

!     if(lprnt .and. i == npr) then
!     write(7000,*)'k=',k,'  taugwci=',taugwci(i,k),
!    &'riloc',riloc(i,k),'riminp=',riminp,' ricrit=',ricrit
!    &,'bruni(i,k)=',bruni(i,k),' deln=',bruni(i,k)
!    &,'basicui(i,k)=',basicui(i,k),' rimin=',rimin(i,k)
!    &,' dlen=',dlen(i),' rhoi=',rhoi(i,k)
!     endif

          endif
        enddo                     ! end of i=1,npt loop
      enddo                       ! end of k=kcldm,1,-1 loop

      do i=1,npt
        dtfac(i)   = 1.0
      enddo
      do k=1,km
        do i=1,npt
          if (do_gwc(i)) then
            kk = kcldtop(i)
            if (k < kk) then
              taugw(i,k) = (taugwci(i,k+1) - taugwci(i,k)) / dpmid(i,k)
              if (taugw(i,k) /= 0.0) then
                tem  = deltim * taugw(i,k)
                dtfac(i) = min(dtfac(i), abs(velco(i,k)/tem))
              endif
            else
              taugw(i,k) = 0.0
            endif
          else
            taugw(i,k) = 0.0
          endif
        enddo
      enddo

!!!!!! Vertical differentiation
!!!!!!
!> - Calculate wind tendency in direction to the wind vector,zonal
!! wind tendency and meridional wind tendency above the cloud top
!! level due to convectively generated gravity waves.

      do k=1,km
        do i=1,npt
          if (do_gwc(i)) then
            kk = kcldtop(i)
            if (k < kk) then
!             wtgwc       = (taugwci(i,k+1) - taugwci(i,k)) / dpmid(i,k)
              wtgwc       = taugw(i,k) * dtfac(i)
              utgwcl(i,k) = wtgwc * cosphi(i)
              vtgwcl(i,k) = wtgwc * sinphi(i)
            else
              utgwcl(i,k) = zero
              vtgwcl(i,k) = zero
            endif
!     if(lprnt .and. i == npr) then
!     write(7000,*)'k=',k,' wtgwc=',wtgwc,' taugwci=',taugwci(i,k),
!    &taugwci(i,k+1),' dpmid=',dpmid(i,k),' cosphi=',cosphi(i),
!    & ' sinphi=',sinphi(i),' utgwcl=',utgwcl(i,k),
!    &'vtgwcl=',vtgwcl(i,k),' dtfac=',dtfac(i)
!     endif
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!
!  Calculate momentum flux = stress deposited above cloup top
!  Apply equal amount with opposite sign within cloud
!
!-----------------------------------------------------------------------

      do i=1,npt
        xstress(i) = zero
        ystress(i) = zero
      enddo
      do k=1,kcldm
        do i=1,npt
          if (do_gwc(i)) then
            xstress(i) = xstress(i) + utgwcl(i,k)*dpmid(i,k)
            ystress(i) = ystress(i) + vtgwcl(i,k)*dpmid(i,k)
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!        ALT 1      ONLY UPPERMOST LAYER
!-----------------------------------------------------------------------

!     kk = kcldtop(i)
!     tem1 = g / dpmid(i,kk)
!     utgwc(i,kk) = - tem1 * xstress
!     vtgwc(i,kk) = - tem1 * ystress

!-----------------------------------------------------------------------
!        ALT 2      SIN(KT-KB)
!-----------------------------------------------------------------------

      do i=1,npt
        if (do_gwc(i)) then
          wrk(i) = 0.5 * pi / (pint(i,kcldbot(i)+1)-pint(i,kcldtop(i)))
        endif
      enddo
      do k=1,km
        do i=1,npt
          if (do_gwc(i)) then
            kk = kcldtop(i)
            if (k >= kk .and. k <= kcldbot(i)) then
              p1 = sin(wrk(i) * (pint(i,k)  -pint(i,kk)))
              p2 = sin(wrk(i) * (pint(i,k+1)-pint(i,kk)))
              tem = - (p2-p1) / dpmid(i,k)
              utgwcl(i,k) = tem*xstress(i)
              vtgwcl(i,k) = tem*ystress(i)
            endif
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!        ALT 3      FROM KT to KB  PROPORTIONAL TO CONV HEATING
!-----------------------------------------------------------------------

!     do k=kcldtop(i),kcldbot(i)
!     p1=cumchr(i,k)
!     p2=cumchr(i,k+1)
!     utgwcl(i,k) = - g*xstress*(p1-p2)/dpmid(i,k)
!     enddo

!-----------------------------------------------------------------------
!
!  The GWDC should accelerate the zonal and meridional wind in the
!  opposite direction of the previous zonal and meridional wind,
!  respectively
!
!-----------------------------------------------------------------------

!     do k=1,kcldtop(i)-1

!      if (utgwcl(i,k)*u(i,k) .gt. 0.0) then

!-------------------- x-component-------------------

!       write(6,'(a)')
!    +  '(GWDC) WARNING: The GWDC should accelerate the zonal wind '
!       write(6,'(a,a,i3,a,i3)')
!    +  'in the opposite direction of the previous zonal wind',
!    +  ' at I = ',i,' and J = ',lat
!       write(6,'(4(1x,e17.10))') u(i,kk),v(i,kk),u(i,k),v(i,k)
!       write(6,'(a,1x,e17.10))') 'Vcld . V =',
!    +  u(i,kk)*u(i,k)+v(i,kk)*v(i,k)

!       if(u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k).gt.0.0)then
!       do k1=1,km
!         write(6,'(i2,36x,2(1x,e17.10))')
!    +             k1,taugwcxi(i,k1),taugwci(i,k1)
!         write(6,'(i2,2(1x,e17.10))') k1,utgwcl(i,k1),u(i,k1)
!       end do
!       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwcxi(i,km+1)
!       end if

!-------------------- Along wind at cloud top -----

!       do k1=1,km
!         write(6,'(i2,36x,2(1x,e17.10))')
!    +             k1,taugwci(i,k1)
!         write(6,'(i2,2(1x,e17.10))') k1,wtgwc(i,k1),basicum(i,k1)
!       end do
!       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwci(i,km+1)

!      end if

!      if (vtgwc(i,k)*v(i,k) .gt. 0.0) then
!       write(6,'(a)')
!    +  '(GWDC) WARNING: The GWDC should accelerate the meridional wind'
!       write(6,'(a,a,i3,a,i3)')
!    +  'in the opposite direction of the previous meridional wind',
!    +  ' at I = ',i,' and J = ',lat
!       write(6,'(4(1x,e17.10))') u(i,kcldtop(i)),v(i,kcldtop(i)),
!    +                            u(i,k),v(i,k)
!       write(6,'(a,1x,e17.10))') 'Vcld . V =',
!    +                    u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k)
!       if(u(i,kcldtop(i))*u(i,k)+v(i,kcldtop(i))*v(i,k).gt.0.0)then
!       do k1=1,km
!         write(6,'(i2,36x,2(1x,e17.10))')
!    +                        k1,taugwcyi(i,k1),taugwci(i,k1)
!         write(6,'(i2,2(1x,e17.10))') k1,vtgwc(i,k1),v(i,k1)
!       end do
!       write(6,'(i2,36x,1x,e17.10)') (km+1),taugwcyi(i,km+1)
!       end if
!      end if

!     enddo

!1000 continue


!***********************************************************************

!     if (lprnt) then
!       if (fhour.ge.fhourpr) then
!-------- UTGWC VTGWC ----------
!         write(*,9220)
!         do ilev=1,km
!           write(*,9221) ilev,(86400.*utgwcl(ipr,ilev)),
!    +                         (86400.*vtgwcl(ipr,ilev))
!         enddo
!       endif
!     endif

!9220 format(//,14x,'TENDENCY DUE TO GWDC',//,
!    +' ILEV',6x,'UTGWC',7x,'VTGWC',/)
!9221 format(i4,2(2x,f10.3))

!-----------------------------------------------------------------------
!
!  For GWDC performance analysis
!
!-----------------------------------------------------------------------

!     do k = 1, kk-1
!       do i = 1, nct

!         kk = kcldtop(i)

!         if ( (abs(taugwci(i,kk)) > taumin) ) then

!           gwdcloc(i) = one

!        if ( abs(taugwci(i,k)-taugwci(i,kk)) > taumin ) then
!         break(i) = 1.0
!         go to 2000
!        endif
!       enddo
!2000   continue

!       do k = 1, kk-1

!        if ( ( abs(taugwci(i,k)).lt.taumin ) .and.
!    &        ( abs(taugwci(i,k+1)).gt.taumin ) .and.
!    &        ( basicum(i,k+1)*basicum(i,k) .lt. 0. ) ) then
!         critic(i) = 1.0
!         print *,i,k,' inside GWDC  taugwci(k) = ',taugwci(i,k)
!         print *,i,k+1,' inside GWDC  taugwci(k+1) = ',taugwci(i,k+1)
!         print *,i,k,' inside GWDC  basicum(k) = ',basicum(i,k)
!         print *,i,k+1,' inside GWDC  basicum(k+1) = ',basicum(i,k+1)
!         print *,i,' inside GWDC  critic = ',critic(i)
!         goto 2010
!        endif
!       enddo
!2010   continue

!      endif

!     enddo

!-----------------------------------------------------------------------
!> - Convert back local convective GWD tendency arrays to GFS model
!!    vertical indices.
!  Outgoing (FU1,FV1)=(utgwc,vtgwc)
!-----------------------------------------------------------------------

      do k=1,km
        k1 = km - k + 1
        do i=1,npt
          ii = ipt(i)
          utgwc(ii,k1) = utgwcl(i,k)

          vtgwc(ii,k1) = vtgwcl(i,k)

!         brunm(ii,kk) = brunm(i,k)
!         brunm(i,k)  = tem

!         rhom(ii,kk) = rhom(i,k)

        enddo
!       if (lprnt) write(7000,*)' k=',k,' k1=',k1,' utgwc='
!    &, utgwc(ipr,k1),' vtgwc=',vtgwc(ipr,k1)
      enddo
      do i=1,npt
        ii = ipt(i)
        tauctx(ii) = tauctxl(i)
        taucty(ii) = tauctyl(i)
       enddo

!     if (lprnt) then
!       if (fhour.ge.fhourpr) then
!-------- UTGWC VTGWC ----------
!         write(*,9225)
!         do ilev=km,1,-1
!           write(*,9226) ilev,(86400.*fu1(ipr,ilev)),
!    +                         (86400.*fv1(ipr,ilev))
!         enddo
!       endif
!     endif

!9225 format(//,14x,'TENDENCY DUE TO GWDC - TO GBPHYS',//,
!    +' ILEV',6x,'UTGWC',7x,'VTGWC',/)
!9226 format(i4,2(2x,f10.3))

      deallocate (kcldtop,kcldbot,do_gwc)
      deallocate (tauctxl,  tauctyl, dtfac,
!    &            gwdcloc, break, critic,   cosphi,
     &            gwdcloc, break,           cosphi,
     &            sinphi,         xstress,  ystress,
     &            dlen,    ucltop, vcltop,  gqmcldlen, wrk)

      deallocate (plnint,          taugwci, velco,
     &            bruni, rhoi,     basicui,
     &            ti,       riloc, rimin,    pint)

      deallocate (plnmid, utgwcl, vtgwcl, basicum, u, v, t,
     &            pmid,   dpmid,  brunm,  rhom, taugw)

      return
      end subroutine gwdc_run
!> @}

! \brief Brief description of the subroutine
!
!> \section arg_table_gwdc_finalize Argument Table
!!
      subroutine gwdc_finalize()
      end subroutine gwdc_finalize

      end module gwdc



      module gwdc_post

      contains

! \brief Brief description of the subroutine
!
!> \section arg_table_gwdc_post_init Argument Table
!!
      subroutine gwdc_post_init()
      end subroutine gwdc_post_init

! \brief Brief description of the subroutine
!
!> \section arg_table_gwdc_post_run Argument Table
!! | local_name     | standard_name                                                   | long_name                                                                | units      | rank | type      | kind      | intent | optional |
!! |----------------|-----------------------------------------------------------------|--------------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                          | horizontal loop extent                                                   | count      |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                              | number of vertical layers                                                | count      |    0 | integer   |           | in     | F        |
!! | lssav          | flag_diagnostics                                                | flag for calculating diagnostic fields                                   | flag       |    0 | logical   |           | in     | F        |
!! | ldiag3d        | flag_diagnostics_3D                                             | flag for calculating 3-D diagnostic fields                               | flag       |    0 | logical   |           | in     | F        |
!! | dtf            | time_step_for_dynamics                                          | dynamics time step                                                       | s          |    0 | real      | kind_phys | in     | F        |
!! | dtp            | time_step_for_physics                                           | physics time step                                                        | s          |    0 | real      | kind_phys | in     | F        |
!! | con_cp         | specific_heat_of_dry_air_at_constant_pressure                   | specific heat of dry air at constant pressure                            | J kg-1 K-1 |    0 | real      | kind_phys | in     | F        |
!! | tauctx         | instantaneous_x_stress_due_to_gravity_wave_drag                 | zonal stress at cloud top due to convective gravity wave drag            | Pa         |    1 | real      | kind_phys | in     | F        |
!! | taucty         | instantaneous_y_stress_due_to_gravity_wave_drag                 | meridional stress at cloud top due to convective gravity wave drag       | Pa         |    1 | real      | kind_phys | in     | F        |
!! | gwdcu          | tendency_of_x_wind_due_to_convective_gravity_wave_drag          | zonal wind tendency due to convective gravity wave drag                  | m s-2      |    2 | real      | kind_phys | in     | F        |
!! | gwdcv          | tendency_of_y_wind_due_to_convective_gravity_wave_drag          | meridional wind tendency due to convective gravity wave drag             | m s-2      |    2 | real      | kind_phys | in     | F        |
!! | dugwd          | time_integral_of_x_stress_due_to_gravity_wave_drag              | integral over time of zonal stress due to gravity wave drag              | Pa s       |    1 | real      | kind_phys | inout  | F        |
!! | dvgwd          | time_integral_of_y_stress_due_to_gravity_wave_drag              | integral over time of meridional stress due to gravity wave drag         | Pa s       |    1 | real      | kind_phys | inout  | F        |
!! | du3dt          | cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag | cumulative change in zonal wind due to convective gravity wave drag      | m s-1      |    2 | real      | kind_phys | inout  | F        |
!! | dv3dt          | cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag | cumulative change in meridional wind due to convective gravity wave drag | m s-1      |    2 | real      | kind_phys | inout  | F        |
!! | gu0            | x_wind_updated_by_physics                                       | updated zonal wind                                                       | m s-1      |    2 | real      | kind_phys | inout  | F        |
!! | gv0            | y_wind_updated_by_physics                                       | updated meridional wind                                                  | m s-1      |    2 | real      | kind_phys | inout  | F        |
!! | gt0            | air_temperature_updated_by_physics                              | updated air temperature                                                  | K          |    2 | real      | kind_phys | inout  | F        |
!! | errmsg         | error_message                                                   | error message for error handling in CCPP                                 | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                      | error flag for error handling in CCPP                                    | flag       |    0 | integer   |           | out    | F        |
!!
      subroutine gwdc_post_run(                                         &
     &  im, levs, lssav, ldiag3d, dtf, dtp, con_cp,                     &
     &  tauctx, taucty, gwdcu, gwdcv,                                   &
     &  dugwd, dvgwd, du3dt, dv3dt, gu0, gv0, gt0,                      &
     &  errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs
      logical, intent(in) :: lssav, ldiag3d
      real(kind=kind_phys), intent(in) :: dtf, dtp, con_cp
      real(kind=kind_phys), intent(in) ::                               &
     &  tauctx(:), taucty(:), gwdcu(:,:), gwdcv(:,:)

      real(kind=kind_phys), intent(inout) ::                            &
     &  dugwd(:), dvgwd(:), du3dt(:,:), dv3dt(:,:),                     &
     &  gu0(:,:), gv0(:,:), gt0(:,:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k
      real(kind=kind_phys) :: eng0, eng1

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- ...  write out cloud top stress and wind tendencies

      if (lssav) then
        dugwd(:) = dugwd(:) + tauctx(:)*dtf
        dvgwd(:) = dvgwd(:) + taucty(:)*dtf

        if (ldiag3d) then
          du3dt(:,:) = du3dt(:,:) + gwdcu(:,:)  * dtf
          dv3dt(:,:) = dv3dt(:,:) + gwdcv(:,:)  * dtf
        endif
      endif   ! end if_lssav

!  --- ...  update the wind components with  gwdc tendencies

      do k = 1, levs
        do i = 1, im
          eng0     = 0.5*(gu0(i,k)*gu0(i,k) + gv0(i,k)*gv0(i,k))
          gu0(i,k) = gu0(i,k) + gwdcu(i,k) * dtp
          gv0(i,k) = gv0(i,k) + gwdcv(i,k) * dtp
          eng1     = 0.5*(gu0(i,k)*gu0(i,k) + gv0(i,k)*gv0(i,k))
          gt0(i,k) = gt0(i,k) + (eng0-eng1)/(dtp*con_cp)
        enddo
!         if (lprnt) write(7000,*)' gu0=',gu0(ipr,k),' gwdcu=',
!    &gwdcu(ipr,k), ' gv0=', gv0(ipr,k),' gwdcv=',gwdcv(ipr,k)
!    &,' k=',k
      enddo

      end subroutine gwdc_post_run

! \brief Brief description of the subroutine
!
!> \section arg_table_gwdc_post_finalize Argument Table
!!
      subroutine gwdc_post_finalize()
      end subroutine gwdc_post_finalize

      end module gwdc_post

