!>  \file mfshalcnv.f
!!  This file contains the Scale-Aware mass flux Shallow Convection scheme.

      module sasas_shal_post
      contains

!! \section arg_table_sasas_shal_post_run Argument Table
!! | local_name     | standard_name                                            | long_name                                                            | units   | rank | type             |    kind   | intent | optional |
!! |----------------|----------------------------------------------------------|----------------------------------------------------------------------|---------|------|------------------|-----------|--------|----------|
!! | frain          | dynamics_to_physics_timestep_ratio                       | ratio of dynamics timestep to physics timestep                       | none    |    0 | real             | kind_phys | in     | F        |
!! | rain1          | lwe_thickness_of_shallow_convective_precipitation_amount | shallow convective rainfall amount on physics timestep               | m       |    1 | real             | kind_phys | in     | F        |
!! | cnvc           | convective_cloud_cover                                   | convective cloud cover                                               | frac    |    2 | real             | kind_phys | in     | F        |
!! | cnvw           | convective_cloud_water_specific_humidity                 | convective cloud water specific humidity                             | kg kg-1 |    2 | real             | kind_phys | in     | F        |
!! | Model          | FV3-GFS_Control_type                                     | Fortran DDT containing FV3-GFS model control parameters              | DDT     |    0 | GFS_control_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                        | Fortran DDT containing FV3-GFS grid and interpolation related data   | DDT     |    0 | GFS_grid_type    |           | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                        | Fortran DDT containing FV3-GFS fields targeted for diagnostic output | DDT     |    0 | GFS_diag_type    |           | inout  | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                         | Fortran DDT containing FV3-GFS miscellaneous data                    | DDT     |    0 | GFS_tbd_type     |           | inout  | F        |
!! | errmsg         | error_message                                            | error message for error handling in CCPP                             | none    |    0 | character        | len=*     | out    | F        |
!! | errflg         | error_flag                                               | error flag for error handling in CCPP                                | flag    |    0 | integer          |           | out    | F        |
!!
      subroutine sasas_shal_post_run (frain, rain1, cnvc, cnvw, Model,   &
     &                                Grid, Diag, Tbd, errmsg, errflg)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type,              &
     &     GFS_grid_type, GFS_diag_type, GFS_tbd_type
        implicit none
!
        type(GFS_grid_type),            intent(in) :: Grid
        type(GFS_control_type),         intent(in) :: Model
        type(GFS_diag_type),         intent(inout) :: Diag
        type(GFS_tbd_type),          intent(inout) :: Tbd

        real(kind=kind_phys), intent(in) :: frain
        real(kind=kind_phys), dimension(size(Grid%xlon,1)),             &
     &     intent(in) :: rain1
        real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs),  &
     &     intent(in) :: cnvw, cnvc

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        real(kind=kind_phys), dimension(size(Grid%xlon,1)) :: raincs
        integer :: num2, num3

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        raincs(:)     = frain * rain1(:)
        Diag%rainc(:) = Diag%rainc(:) + raincs(:)
        if (Model%lssav) then
          Diag%cnvprcp(:) = Diag%cnvprcp(:) + raincs(:)
        endif
        if ((Model%shcnvcw) .and. (Model%num_p3d == 4) .and.            &
     &                             (Model%npdf3d == 3)) then
          num2 = Model%num_p3d + 2
          num3 = num2 + 1
          Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
          Tbd%phy_f3d(:,:,num3) = cnvc(:,:)
        elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then
          num2 = Model%num_p3d + 1
          Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
        endif

      end subroutine sasas_shal_post_run

!! \section arg_table_sasas_shal_post_init Argument Table
!!
      subroutine sasas_shal_post_init ()
      end subroutine sasas_shal_post_init

!! \section arg_table_sasas_shal_post_finalize Argument Table
!!
      subroutine sasas_shal_post_finalize ()
      end subroutine sasas_shal_post_finalize

      end module sasas_shal_post


      module sasas_shal
      contains


!> \section arg_table_sasas_shal_init Argument Table
!!
      subroutine sasas_shal_init
      end subroutine sasas_shal_init

! \defgroup SAMF_shal GFS Scale-Aware Mass-Flux Shallow Convection
!>\defgroup SAMF_shal_main GFS mfshalcnv Main
!! @{
!> \brief The subroutine contains the entirety of the SAMF shallow convection scheme.
!! The scale-aware mass-flux shallow  convection
!! scheme is an updated version of the previous mass-flux shallow
!! convection scheme with scale and aerosol awareness and
!! parameterizes the effect of shallow convection on the environment.
!!
!! \section arg_table_sasas_shal_run Argument Table
!! | local_name     | standard_name                                             | long_name                                              | units   | rank | type      |    kind   | intent | optional |
!! |----------------|-----------------------------------------------------------|--------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                    | horizontal loop extent                                 | count   |    0 | integer   |           | in     | F        |
!! | ix             | horizontal_dimension                                      | horizontal dimension                                   | count   |    0 | integer   |           | in     | F        |
!! | km             | vertical_dimension                                        | vertical layer dimension                               | count   |    0 | integer   |           | in     | F        |
!! | delt           | time_step_for_physics                                     | physics time step                                      | s       |    0 | real      | kind_phys | in     | F        |
!! | delp           | air_pressure_difference_between_midlayers                 | pres(k) - pres(k+1)                                    | Pa      |    2 | real      | kind_phys | in     | F        |
!! | prslp          | air_pressure                                              | mean layer pressure                                    | Pa      |    2 | real      | kind_phys | in     | F        |
!! | psp            | surface_air_pressure                                      | surface pressure                                       | Pa      |    1 | real      | kind_phys | in     | F        |
!! | phil           | geopotential                                              | layer geopotential                                     | m2 s-2  |    2 | real      | kind_phys | in     | F        |
!! | ql1            | cloud_ice_specific_humidity                               | cloud ice specific humidity                            | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | ql2            | cloud_liquid_water_specific_humidity                      | cloud water specific humidity                          | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | q1             | water_vapor_specific_humidity_updated_by_physics          | updated vapor specific humidity                        | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | t1             | air_temperature_updated_by_physics                        | updated temperature                                    | K       |    2 | real      | kind_phys | inout  | F        |
!! | u1             | x_wind_updated_by_physics                                 | updated x-direction wind                               | m s-1   |    2 | real      | kind_phys | inout  | F        |
!! | v1             | y_wind_updated_by_physics                                 | updated y-direction wind                               | m s-1   |    2 | real      | kind_phys | inout  | F        |
!! | rn             | lwe_thickness_of_shallow_convective_precipitation_amount  | shallow convective rainfall amount on physics timestep | m       |    1 | real      | kind_phys | out    | F        |
!! | kbot           | vertical_index_at_cloud_base                              | index at cloud base                                    | index   |    1 | integer   |           | inout  | F        |
!! | ktop           | vertical_index_at_cloud_top                               | index at cloud top                                     | index   |    1 | integer   |           | inout  | F        |
!! | kcnv           | flag_deep_convection                                      | deep convection: 0=no, 1=yes                           | flag    |    1 | integer   |           | inout  | F        |
!! | islimsk        | sea_land_ice_mask                                         | landmask: sea/land/ice=0/1/2                           | flag    |    1 | integer   |           | in     | F        |
!! | garea          | cell_area                                                 | grid cell area                                         | m2      |    1 | real      | kind_phys | in     | F        |
!! | dot            | omega                                                     | layer mean vertical velocity                           | Pa s-1  |    2 | real      | kind_phys | in     | F        |
!! | ncloud         | number_of_hydrometeors                                    | number of hydrometeors                                 | count   |    0 | integer   |           | in     | F        |
!! | hpbl           | atmosphere_boundary_layer_thickness                       | PBL top height                                         | m       |    1 | real      | kind_phys | in     | F        |
!! | ud_mf          | instantaneous_atmosphere_updraft_convective_mass_flux     | (updraft mass flux) * delt                             | kg m-2  |    2 | real      | kind_phys | out    | F        |
!! | dt_mf          | instantaneous_atmosphere_detrainment_convective_mass_flux | (detrainment mass flux) * delt                         | kg m-2  |    2 | real      | kind_phys | out    | F        |
!! | cnvw           | convective_cloud_water_specific_humidity                  | convective cloud water specific humidity               | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | cnvc           | convective_cloud_cover                                    | convective cloud cover                                 | frac    |    2 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                             | error message for error handling in CCPP               | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                | error flag for error handling in CCPP                  | flag    |    0 | integer   |           | out    | F        |
!!
!>  \section general_mfshal GFS SAMF Shallow Convection Scheme General Algorithm
!!  -# Compute preliminary quantities needed for the static and feedback
!!  control portions of the algorithm.
!!  -# Perform calculations related to the updraft of the entraining/detraining
!!  cloud model ("static control").
!!  -# The cloud base mass flux is obtained using the cumulus updraft
!!  velocity averaged ove the whole cloud depth.
!!  -# Calculate the tendencies of the state variables (per unit cloud
!!  base mass flux) and the cloud base mass flux.
!!  -# For the "feedback control", calculate updated values of the state
!!  variables by multiplying the cloud base mass flux and the tendencies
!! calculated per unit cloud base mass flux from the static control.
!!
!!  \section detailed_mfshal GFS SAMF Shallow Convection Scheme Detailed Algorithm
!!  @{
      subroutine sasas_shal_run (im,ix,km,delt,delp,prslp,psp,phil,ql1, &
     &     ql2,q1,t1,u1,v1,rn,kbot,ktop,kcnv,islimsk,garea,             &
     &     dot,ncloud,hpbl,ud_mf,dt_mf,cnvw,cnvc,errmsg,errflg)

      use machine , only : kind_phys
      use funcphys , only : fpvs
      use physcons, grav => con_g, cp => con_cp, hvap => con_hvap
     &,             rv => con_rv, fv => con_fvirt, t0c => con_t0c
     &,             rd => con_rd, cvap => con_cvap, cliq => con_cliq
     &,             eps => con_eps, epsm1 => con_epsm1
      implicit none
!
! In the current NCEP spectral model im <= ix for reduced grid numbers
! near the pole and a parallel computing. For FV3, im=ix.
!
! Interface variables
      integer, intent(in) :: im, ix,  km, ncloud
      integer, intent(inout) :: kbot(im), ktop(im), kcnv(im)
!
      real(kind=kind_phys), intent(in) :: delt
      real(kind=kind_phys), intent(in) :: psp(im), delp(ix,km),         &
     &                     prslp(ix,km), phil(ix,km)
      real(kind=kind_phys), intent(inout) :: ql1(ix,km), ql2(ix,km),    &
     &                     q1(ix,km), t1(ix,km), u1(ix,km), v1(ix,km)
      real(kind=kind_phys), intent(out) :: rn(im)
!
      integer, intent(in) :: islimsk(im)
      real(kind=kind_phys), intent(in) :: garea(im), dot(ix,km),        &
     &                     hpbl(im)
      real(kind=kind_phys), intent(out) :: cnvw(ix,km), cnvc(ix,km),    &
     &                     ud_mf(im,km), dt_mf(im,km)                     ! hchuang code change mass flux output
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
! Local variables
      real(kind=kind_phys) ps(im), del(ix,km), prsl(ix,km)
!
      integer              i,j,indx, k, kk, km1, n
      integer              kpbl(im)
!
      real(kind=kind_phys) dellat,  delta,
     &                     c0l,     c0s,     d0,
     &                     c1,      asolfac,
     &                     desdt,   dp,
     &                     dq,      dqsdp,   dqsdt,   dt,
     &                     dt2,     dtmax,   dtmin,   dxcrt,
     &                     dv1h,    dv2h,    dv3h,
     &                     dv1q,    dv2q,    dv3q,
     &                     dz,      dz1,     e1,      clam,
     &                     el2orc,  elocp,   aafac,   cm,
     &                     es,      etah,    h1,
     &                     evef,    evfact,  evfactl, fact1,
     &                     fact2,   factor,  dthk,
     &                     g,       gamma,   pprime,  betaw,
     &                     qlk,     qrch,    qs,
     &                     rfact,   shear,   tfac,
     &                     val,     val1,    val2,
     &                     w1,      w1l,     w1s,     w2,
     &                     w2l,     w2s,     w3,      w3l,
     &                     w3s,     w4,      w4l,     w4s,
     &                     rho,     tem,     tem1,    tem2,
     &                     ptem,    ptem1,
     &                     pgcon
!
      integer              kb(im), kbcon(im), kbcon1(im),
     &                     ktcon(im), ktcon1(im), ktconn(im),
     &                     kbm(im), kmax(im)
!
      real(kind=kind_phys) aa1(im),     cina(im),
     &                     umean(im),  tauadv(im),  gdx(im),
     &                     delhbar(im), delq(im),   delq2(im),
     &                     delqbar(im), delqev(im), deltbar(im),
     &                     deltv(im),   dtconv(im), edt(im),
     &                     pdot(im),    po(im,km),
     &                     qcond(im),   qevap(im),  hmax(im),
     &                     rntot(im),   vshear(im),
     &                     xlamud(im),  xmb(im),    xmbmax(im),
     &                     delubar(im), delvbar(im)
!
      real(kind=kind_phys) c0(im)
c
      real(kind=kind_phys) crtlamd
!
      real(kind=kind_phys) cinpcr,  cinpcrmx,  cinpcrmn,
     &                     cinacr,  cinacrmx,  cinacrmn
!
!  parameters for updraft velocity calculation
      real(kind=kind_phys) bet1,    cd1,     f1,      gam1,
     &                     bb1,     bb2,     wucb
cc
c  physical parameters
      parameter(g=grav,asolfac=0.89)
      parameter(elocp=hvap/cp,
     &          el2orc=hvap*hvap/(rv*cp))
      parameter(c0s=0.002,c1=5.e-4,d0=.01)
      parameter(c0l=c0s*asolfac)
!
! asolfac: aerosol-aware parameter based on Lim & Hong (2012)
!      asolfac= cx / c0s(=.002)
!      cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
!      Nccn: CCN number concentration in cm^(-3)
!      Until a realistic Nccn is provided, typical Nccns are assumed
!      as Nccn=100 for sea and Nccn=7000 for land
!
      parameter(cm=1.0,delta=fv)
      parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(dthk=25.)
      parameter(cinpcrmx=180.,cinpcrmn=120.)
!     parameter(cinacrmx=-120.,cinacrmn=-120.)
      parameter(cinacrmx=-120.,cinacrmn=-80.)
      parameter(crtlamd=3.e-4)
      parameter(dtmax=10800.,dtmin=600.)
      parameter(bet1=1.875,cd1=.506,f1=2.0,gam1=.5)
      parameter(betaw=.03,dxcrt=15.e3)
      parameter(h1=0.33333333)
c  local variables and arrays
      real(kind=kind_phys) pfld(im,km),    to(im,km),     qo(im,km),
     &                     uo(im,km),      vo(im,km),     qeso(im,km)
!  for updraft velocity calculation
      real(kind=kind_phys) wu2(im,km),     buo(im,km),    drag(im,km)
      real(kind=kind_phys) wc(im),         scaldfunc(im), sigmagfm(im)
!
c  cloud water
!     real(kind=kind_phys) qlko_ktcon(im), dellal(im,km), tvo(im,km),
      real(kind=kind_phys) qlko_ktcon(im), dellal(im,km),
     &                     dbyo(im,km),    zo(im,km),     xlamue(im,km),
     &                     heo(im,km),     heso(im,km),
     &                     dellah(im,km),  dellaq(im,km),
     &                     dellau(im,km),  dellav(im,km), hcko(im,km),
     &                     ucko(im,km),    vcko(im,km),   qcko(im,km),
     &                     qrcko(im,km),   eta(im,km),
     &                     zi(im,km),      pwo(im,km),    c0t(im,km),
     &                     sumx(im),       tx1(im),       cnvwt(im,km)
!
      logical totflg, cnvflg(im), flg(im)
!
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))
!
c-----------------------------------------------------------------------
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!************************************************************************
!     convert input Pa terms to Cb terms  -- Moorthi
!>  ## Compute preliminary quantities needed for the static and feedback control portions of the algorithm.
!>  - Convert input pressure terms to centibar units.
      ps   = psp   * 0.001
      prsl = prslp * 0.001
      del  = delp  * 0.001
!************************************************************************
!
      km1 = km - 1
c
c  initialize arrays
c
!>  - Initialize column-integrated and other single-value-per-column variable arrays.
      do i=1,im
        cnvflg(i) = .true.
        if(kcnv(i) == 1) cnvflg(i) = .false.
        if(cnvflg(i)) then
          kbot(i)=km+1
          ktop(i)=0
        endif
        rn(i)=0.
        kbcon(i)=km
        ktcon(i)=1
        ktconn(i)=1
        kb(i)=km
        pdot(i) = 0.
        qlko_ktcon(i) = 0.
        edt(i)  = 0.
        aa1(i)  = 0.
        cina(i) = 0.
        vshear(i) = 0.
        gdx(i) = sqrt(garea(i))
      enddo
!!
!>  - Return to the calling routine if deep convection is present or the surface buoyancy flux is negative.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if (totflg) then
        cnvw  = 0.
        cnvc  = 0.
        ud_mf = 0.
        dt_mf = 0.
        return
      endif
!!
!>  - determine aerosol-aware rain conversion parameter over land
      do i=1,im
        if(islimsk(i) == 1) then
           c0(i) = c0l
        else
           c0(i) = c0s
        endif
      enddo
!
!>  - determine rain conversion parameter above the freezing level
!! which exponentially decreases with decreasing temperature from
!! \cite han_et_al_2017 equation 8.
      do k = 1, km
        do i = 1, im
          if(t1(i,k) > 273.16) then
            c0t(i,k) = c0(i)
          else
            tem = d0 * (t1(i,k) - 273.16)
            tem1 = exp(tem)
            c0t(i,k) = c0(i) * tem1
          endif
        enddo
      enddo
!
!>  - Initialize convective cloud water and cloud cover to zero.
      do k = 1, km
        do i = 1, im
          cnvw(i,k) = 0.
          cnvc(i,k) = 0.
        enddo
      enddo
! hchuang code change
!>  - Initialize updraft mass fluxes to zero.
      do k = 1, km
        do i = 1, im
          ud_mf(i,k) = 0.
          dt_mf(i,k) = 0.
        enddo
      enddo
c
      dt2   = delt
!
c  model tunable parameters are all here
      clam    = .3
      aafac   = .1
c     evef    = 0.07
      evfact  = 0.3
      evfactl = 0.3
!
!     pgcon   = 0.7     ! Gregory et al. (1997, QJRMS)
      pgcon   = 0.55    ! Zhang & Wu (2003,JAS)
      w1l     = -8.e-3
      w2l     = -4.e-2
      w3l     = -5.e-3
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
c
c  define top layer for search of the downdraft originating layer
c  and the maximum thetae for updraft
c
!>  - Determine maximum indices for the parcel starting point (kbm) and cloud top (kmax).
      do i=1,im
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!
      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) > 0.70) kbm(i)   = k + 1
          if (prsl(i,k)*tx1(i) > 0.60) kmax(i)  = k + 1
        enddo
      enddo
      do i=1,im
        kbm(i)   = min(kbm(i),kmax(i))
      enddo
c
c  hydrostatic height assume zero terr and compute
c  updraft entrainment rate as an inverse function of height
c
!>  - Calculate hydrostatic height at layer centers assuming a flat
!! surface (no terrain) from the geopotential.
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / g
        enddo
      enddo
!>  - Calculate interface height and the entrainment rate as an inverse
!! function of height.
      do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
        enddo
      enddo
      do i=1,im
        xlamue(i,km) = xlamue(i,km1)
      enddo
c
c  pbl height
c
!>  - Find the index for the PBL top using the PBL height; enforce that
!! it is lower than the maximum parcel starting level.
      do i=1,im
        flg(i) = cnvflg(i)
        kpbl(i)= 1
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i) .and. zo(i,k) <= hpbl(i)) then
            kpbl(i) = k
          else
            flg(i) = .false.
          endif
        enddo
      enddo
      do i=1,im
        kpbl(i)= min(kpbl(i),kbm(i))
      enddo
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c   convert surface pressure to mb from cb
c
!>  - Convert prsl from centibar to millibar, set normalized mass flux
!! to 1, cloud properties to 0, and save model state variables (after
!! advection/turbulence).
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k <= kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
            hcko(i,k) = 0.
            qcko(i,k) = 0.
            qrcko(i,k)= 0.
            ucko(i,k) = 0.
            vcko(i,k) = 0.
            dbyo(i,k) = 0.
            pwo(i,k)  = 0.
            dellal(i,k) = 0.
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k)
            vo(i,k)   = v1(i,k)
!           uo(i,k)   = u1(i,k) * rcs(i)
!           vo(i,k)   = v1(i,k) * rcs(i)
            wu2(i,k)  = 0.
            buo(i,k)  = 0.
            drag(i,k) = 0.
            cnvwt(i,k) = 0.
          endif
        enddo
      enddo
c
c  column variables
c  p is pressure of the layer (mb)
c  t is temperature at t-dt (k)..tn
c  q is specific humidity at t-dt (kg/kg)..qn
c  to is temperature at t+dt (k)... this is after advection and turbulence
c  qo is specific humidity at t+dt (kg/kg)..q1
c
!>  - Calculate saturation specific humidity and enforce minimum moisture values.
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
!           tvo(i,k)  = to(i,k) + delta * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
c
c  compute moist static energy
c
!>  - Calculate moist static energy (heo) and saturation moist static energy (heso).
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)) then
!           tem       = g * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
c           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
c
c  determine level with largest moist static energy within pbl
c  this is the level where updraft starts
c
!> ## Perform calculations related to the updraft of the entraining/detraining cloud model ("static control").
!> - Search in the PBL for the level of maximum moist static energy to start the ascending parcel.
      do i=1,im
         if (cnvflg(i)) then
            hmax(i) = heo(i,1)
            kb(i) = 1
         endif
      enddo
      do k = 2, km
        do i=1,im
          if (cnvflg(i) .and. k <= kpbl(i)) then
            if(heo(i,k) > hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
c
!> - Calculate the temperature, water vapor mixing ratio, and pressure at interface levels.
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
            es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (g * dz + hvap * dqsdp * dp) / (cp * (1. + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
!
!> - Recalculate saturation specific humidity, moist static energy, 
!! saturation moist static energy, and horizontal momentum on interface
!! levels. Enforce minimum specific humidity.
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            heo(i,k)  = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * g * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
          endif
        enddo
      enddo
c
c  look for the level of free convection as cloud base
c
!> - Search below the index "kbm" for the level of free convection (LFC)
!! where the condition \f$h_b > h^*\f$ is first met, where \f$h_b, h^*\f$ 
!! are the state moist static energy at the parcel's starting level and
!! saturation moist static energy, respectively. Set "kbcon" to the
!! index of the LFC.
      do i=1,im
        flg(i)   = cnvflg(i)
        if(flg(i)) kbcon(i) = kmax(i)
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i) .and. k < kbm(i)) then
            if(k > kb(i) .and. heo(i,kb(i)) > heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
c
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon(i) == kmax(i)) cnvflg(i) = .false.
        endif
      enddo
!!
!> - If no LFC, return to the calling routine without modifying state variables.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!> - Determine the vertical pressure velocity at the LFC. After 
!! \cite han_and_pan_2011 , determine the maximum pressure
!! thickness between a parcel's starting level and the LFC. If a parcel
!! doesn't reach the LFC within the critical thickness, then the convective
!! inhibition is deemed too great for convection to be triggered, and the
!! subroutine returns to the calling routine without modifying the state
!! variables.
      do i=1,im
        if(cnvflg(i)) then
!         pdot(i)  = 10.* dot(i,kbcon(i))
          pdot(i)  = 0.01 * dot(i,kbcon(i)) ! Now dot is in Pa/s
        endif
      enddo
c
c   turn off convection if pressure depth between parcel source level
c      and cloud base is larger than a critical value, cinpcr
c
      do i=1,im
        if(cnvflg(i)) then
          if(islimsk(i) == 1) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i) <= w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i) >= -w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = 0.
          endif
          val1    =            -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          ptem = 1. - tem
          ptem1= .5*(cinpcrmx-cinpcrmn)
          cinpcr = cinpcrmx - ptem * ptem1
          tem1 = pfld(i,kb(i)) - pfld(i,kbcon(i))
          if(tem1 > cinpcr) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  specify the detrainment rate for the updrafts
c
!> - The updraft detrainment rate is set constant and equal to the entrainment rate at cloud base.
      do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
!         xlamud(i) = crtlamd
        endif
      enddo
c
c  determine updraft mass flux for the subcloud layers
c
!> - Calculate the normalized mass flux for subcloud and in-cloud layers
!! according to \cite pan_and_wu_1995 equation 1:
!!  \f[
!!  \frac{1}{\eta}\frac{\partial \eta}{\partial z} = \lambda_e - \lambda_d
!!  \f]
!!  where \f$\eta\f$ is the normalized mass flux, \f$\lambda_e\f$ is the
!! entrainment rate and \f$\lambda_d\f$ is the detrainment rate. The
!! normalized mass flux increases upward below the cloud base and decreases
!! upward above.
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k < kbcon(i) .and. k >= kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
              eta(i,k) = eta(i,k+1) / (1. + ptem * dz)
            endif
          endif
        enddo
      enddo
c
c  compute mass flux above cloud base
c
      do i = 1, im
        flg(i) = cnvflg(i)
      enddo
      do k = 2, km1
        do i = 1, im
         if(flg(i))then
           if(k > kbcon(i) .and. k < kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
              if(eta(i,k) <= 0.) then
                kmax(i) = k
                ktconn(i) = k
                kbm(i) = min(kbm(i),kmax(i))
                flg(i) = .false.
              endif
           endif
         endif
        enddo
      enddo
c
c  compute updraft cloud property
c
!> - Set cloud properties equal to the state variables at updraft starting level (kb).
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
        endif
      enddo
c
!  cm is an enhancement factor in entrainment rates for momentum
!
!> - Calculate the cloud properties as a parcel ascends, modified by
!! entrainment and detrainment. Discretization follows Appendix B of
!! \cite grell_1993 . Following 
!! \cite han_and_pan_2006, the convective momentum transport is reduced
!! by the convection-induced pressure gradient force by the constant
!! "pgcon", currently set to 0.55 after 
!! \cite zhang_and_wu_2003 .
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*
     &                     (heo(i,k)+heo(i,k-1)))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
!
              tem  = 0.5 * cm * tem
              factor = 1. + tem
              ptem = tem + pgcon
              ptem1= tem - pgcon
              ucko(i,k) = ((1.-tem)*ucko(i,k-1)+ptem*uo(i,k)
     &                     +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem)*vcko(i,k-1)+ptem*vo(i,k)
     &                     +ptem1*vo(i,k-1))/factor
            endif
          endif
        enddo
      enddo
c
c   taking account into convection inhibition due to existence of
c    dry layers below cloud base
c
!> - With entrainment, recalculate the LFC as the first level where
!! buoyancy is positive. The difference in pressure levels between LFCs
!! calculated with/without entrainment must be less than a threshold
!! (currently 25 hPa). Otherwise, convection is inhibited and the scheme
!! returns to the calling routine without modifying the state variables.
!! This is the subcloud dryness trigger modification discussed in 
!! \cite han_and_pan_2011.
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i) .and. k < kbm(i)) then
          if(k >= kbcon(i) .and. dbyo(i,k) > 0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i) == kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(tem > dthk) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  calculate convective inhibition
c
!> - Calculate additional trigger condition of the convective inhibition
!! (CIN) according to \cite han_et_al_2017 equation 13.
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < kbcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma
     &                 * to(i,k) / hvap
              cina(i) = cina(i) +
!    &                 dz1 * eta(i,k) * (g / (cp * to(i,k)))
     &                 dz1 * (g / (cp * to(i,k)))
     &                 * dbyo(i,k) / (1. + gamma)
     &                 * rfact
              val = 0.
              cina(i) = cina(i) +
!    &                 dz1 * eta(i,k) * g * delta *
     &                 dz1 * g * delta *
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
!> - Turn off convection if the CIN is less than a critical value (cinacr)
!! which is inversely proportional to the large-scale vertical velocity.
      do i = 1, im
        if(cnvflg(i)) then
!
          if(islimsk(i) == 1) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i) <= w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i) >= -w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = 0.
          endif

          val1    =            -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          tem = 1. - tem
          tem1= .5*(cinacrmx-cinacrmn)
          cinacr = cinacrmx - tem * tem1
!
!         cinacr = cinacrmx
          if(cina(i) < cinacr) cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  determine first guess cloud top as the level of zero buoyancy
c    limited to the level of P/Ps=0.7
c
!> - Calculate the cloud top as the first level where parcel buoyancy
!! becomes negative; the maximum possible value is at \f$p=0.7p_{sfc}\f$.
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) ktcon(i) = kbm(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i) .and. k < kbm(i)) then
          if(k > kbcon1(i) .and. dbyo(i,k) < 0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
c
c  specify upper limit of mass flux at cloud base
c
!> - Calculate the maximum value of the cloud base mass flux using the
!! CFL-criterion-based formula of \cite han_and_pan_2011,
!! equation 7.
      do i = 1, im
        if(cnvflg(i)) then
!         xmbmax(i) = .1
!
          k = kbcon(i)
          dp = 1000. * del(i,k)
          xmbmax(i) = dp / (g * dt2)
!
!         tem = dp / (g * dt2)
!         xmbmax(i) = min(tem, xmbmax(i))
        endif
      enddo
c
c  compute cloud moisture property and precipitation
c
!> - Set cloud moisture property equal to the enviromental moisture at
!! updraft starting level (kb).
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          qcko(i,kb(i)) = qo(i,kb(i))
          qrcko(i,kb(i)) = qo(i,kb(i))
        endif
      enddo
!> - Calculate the moisture content of the entraining/detraining parcel
!! (qcko) and the value it would have if just saturated (qrch), according
!! to equation A.14 in \cite grell_1993 . Their difference
!! is the amount of convective cloud water (qlk = rain + condensate).
!! Determine the portion of convective cloud water that remains suspended
!! and the portion that is converted into convective precipitation (pwo).
!! Calculate and save the negative cloud work function (aa1) due to water
!! loading. Above the level of minimum moist static energy, some of the
!! cloud water is detrained into the grid-scale cloud water from every
!! cloud layer with a rate of 0.0005 \f$m^{-1}\f$ (dellal).
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*
     &                     (qo(i,k)+qo(i,k-1)))/factor
              qrcko(i,k) = qcko(i,k)
cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
c
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
c
c  below lfc check if there is excess moisture to release latent heat
c
              if(k >= kbcon(i) .and. dq > 0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud > 0) then
                  dp = 1000. * del(i,k)
                  ptem = c0t(i,k) + c1
                  qlk = dq / (eta(i,k) + etah * ptem * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0t(i,k) * dz)
                endif
                buo(i,k) = buo(i,k) - g * qlk
                qcko(i,k)= qlk + qrch
                pwo(i,k) = etah * c0t(i,k) * dz * qlk
                cnvwt(i,k) = etah * qlk * g / dp
              endif
!
!  compute buoyancy and drag for updraft velocity
!
              if(k >= kbcon(i)) then
                rfact =  1. + delta * cp * gamma
     &                   * to(i,k) / hvap
                buo(i,k) = buo(i,k) + (g / (cp * to(i,k)))
     &                   * dbyo(i,k) / (1. + gamma)
     &                   * rfact
                val = 0.
                buo(i,k) = buo(i,k) + g * delta *
     &                     max(val,(qeso(i,k) - qo(i,k)))
                drag(i,k) = max(xlamue(i,k),xlamud(i))
              endif
!
            endif
          endif
        enddo
      enddo
c
c  calculate cloud work function
c
!     do k = 2, km1
!       do i = 1, im
!         if (cnvflg(i)) then
!           if(k >= kbcon(i) .and. k < ktcon(i)) then
!             dz1 = zo(i,k+1) - zo(i,k)
!             gamma = el2orc * qeso(i,k) / (to(i,k)**2)
!             rfact =  1. + delta * cp * gamma
!    &                 * to(i,k) / hvap
!             aa1(i) = aa1(i) +
!!   &                 dz1 * eta(i,k) * (g / (cp * to(i,k)))
!    &                 dz1 * (g / (cp * to(i,k)))
!    &                 * dbyo(i,k) / (1. + gamma)
!    &                 * rfact
!             val = 0.
!             aa1(i) = aa1(i) +
!!   &                 dz1 * eta(i,k) * g * delta *
!    &                 dz1 * g * delta *
!    &                 max(val,(qeso(i,k) - qo(i,k)))
!           endif
!         endif
!       enddo
!     enddo
!     do i = 1, im
!       if(cnvflg(i) .and. aa1(i) <= 0.) cnvflg(i) = .false.
!     enddo
!
!  calculate cloud work function
!
!> - Calculate the cloud work function according to 
!! \cite pan_and_wu_1995 equation 4:
!!  \f[
!!  A_u=\int_{z_0}^{z_t}\frac{g}{c_pT(z)}\frac{\eta}{1 + \gamma}[h(z)-h^*(z)]dz
!!  \f]
!! (discretized according to \cite grell_1993 equation B.10
!! using B.2 and B.3 of \cite arakawa_and_schubert_1974
!! and assuming \f$\eta=1\f$) where \f$A_u\f$ is the updraft cloud work function,
!! \f$z_0\f$ and \f$z_t\f$ are cloud base and cloud top, respectively,
!! \f$\gamma = \frac{L}{c_p}\left(\frac{\partial \overline{q_s}}{\partial T}\right)_p\f$
!! and other quantities are previously defined.
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k >= kbcon(i) .and. k < ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              aa1(i) = aa1(i) + buo(i,k) * dz1
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i) .and. aa1(i) <= 0.) cnvflg(i) = .false.
      enddo
!!
!> - If the updraft cloud work function is negative, convection does not
!! occur, and the scheme returns to the calling routine.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  estimate the onvective overshooting as the level
c    where the [aafac * cloud work function] becomes zero,
c    which is the final cloud top
c    limited to the level of P/Ps=0.7
c
!> - Continue calculating the cloud work function past the point of
!! neutral buoyancy to represent overshooting according to 
!! \cite han_and_pan_2011 . Convective overshooting stops when
!! \f$ cA_u < 0\f$ where \f$c\f$ is currently 10%, or when 10% of the
!! updraft cloud work function has been consumed by the stable buoyancy
!! force. Overshooting is also limited to the level where \f$p=0.7p_{sfc}\f$.
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = aafac * aa1(i)
        endif
      enddo
c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon1(i) = kbm(i)
      enddo
      do k = 2, km1
        do i = 1, im
          if (flg(i)) then
            if(k >= ktcon(i) .and. k < kbm(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + delta * cp * gamma
     &                 * to(i,k) / hvap
              aa1(i) = aa1(i) +
!    &                 dz1 * eta(i,k) * (g / (cp * to(i,k)))
     &                 dz1 * (g / (cp * to(i,k)))
     &                 * dbyo(i,k) / (1. + gamma)
     &                 * rfact
!             val = 0.
!             aa1(i) = aa1(i) +
!!   &                 dz1 * eta(i,k) * g * delta *
!    &                 dz1 * g * delta *
!    &                 max(val,(qeso(i,k) - qo(i,k)))
              if(aa1(i) < 0.) then
                ktcon1(i) = k
                flg(i) = .false.
              endif
            endif
          endif
        enddo
      enddo
c
c  compute cloud moisture property, detraining cloud water
c    and precipitation in overshooting layers
c
!> - For the overshooting convection, calculate the moisture content of
!! the entraining/detraining parcel as before. Partition convective cloud
!! water and precipitation and detrain convective cloud water in the
!! overshooting layers.
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k >= ktcon(i) .and. k < ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
cj
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              qcko(i,k) = ((1.-tem1)*qcko(i,k-1)+tem*0.5*
     &                     (qo(i,k)+qo(i,k-1)))/factor
              qrcko(i,k) = qcko(i,k)
cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
c
c  check if there is excess moisture to release latent heat
c
              if(dq > 0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                if(ncloud > 0) then
                  dp = 1000. * del(i,k)
                  ptem = c0t(i,k) + c1
                  qlk = dq / (eta(i,k) + etah * ptem * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * g / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0t(i,k) * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0t(i,k) * dz * qlk
                cnvwt(i,k) = etah * qlk * g / dp
              endif
            endif
          endif
        enddo
      enddo
!
!  compute updraft velocity square(wu2)
!> - Calculate updraft velocity square(wu2) according to \cite han_et_al_2017 equation 7.
!
!     bb1 = 2. * (1.+bet1*cd1)
!     bb2 = 2. / (f1*(1.+gam1))
!
!     bb1 = 3.9
!     bb2 = 0.67
!
!     bb1 = 2.0
!     bb2 = 4.0
!
      bb1 = 4.0
      bb2 = 0.8
!
      do i = 1, im
        if (cnvflg(i)) then
          k = kbcon1(i)
          tem = po(i,k) / (rd * to(i,k))
          wucb = -0.01 * dot(i,k) / (tem * g)
          if(wucb > 0.) then
            wu2(i,k) = wucb * wucb
          else
            wu2(i,k) = 0.
          endif
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kbcon1(i) .and. k < ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              tem  = 0.25 * bb1 * (drag(i,k)+drag(i,k-1)) * dz
              tem1 = 0.5 * bb2 * (buo(i,k)+buo(i,k-1)) * dz
              ptem = (1. - tem) * wu2(i,k-1)
              ptem1 = 1. + tem
              wu2(i,k) = (ptem + tem1) / ptem1
              wu2(i,k) = max(wu2(i,k), 0.)
            endif
          endif
        enddo
      enddo
!
!  compute updraft velocity averaged over the whole cumulus
!
!> - Calculate the mean updraft velocity within the cloud (wc).
      do i = 1, im
        wc(i) = 0.
        sumx(i) = 0.
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kbcon1(i) .and. k < ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem = 0.5 * (sqrt(wu2(i,k)) + sqrt(wu2(i,k-1)))
              wc(i) = wc(i) + tem * dz
              sumx(i) = sumx(i) + dz
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          if(sumx(i) == 0.) then
             cnvflg(i)=.false.
          else
             wc(i) = wc(i) / sumx(i)
          endif
          val = 1.e-4
          if (wc(i) < val) cnvflg(i)=.false.
        endif
      enddo
c
c exchange ktcon with ktcon1
c
      do i = 1, im
        if(cnvflg(i)) then
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
        endif
      enddo
c
c  this section is ready for cloud water
c
      if(ncloud > 0) then
c
c  compute liquid and vapor separation at cloud top
c
!> - => Separate the total updraft cloud water at cloud top into vapor and condensate.
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k)
     &         + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k) - qrch
c
c  check if there is excess moisture to release latent heat
c
          if(dq > 0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
          endif
        endif
      enddo
      endif
c
c--- compute precipitation efficiency in terms of windshear
c
!! - Calculate the wind shear and precipitation efficiency according to
!! equation 58 in \cite fritsch_and_chappell_1980 :
!! \f[
!! E = 1.591 - 0.639\frac{\Delta V}{\Delta z} + 0.0953\left(\frac{\Delta V}{\Delta z}\right)^2 - 0.00496\left(\frac{\Delta V}{\Delta z}\right)^3
!! \f]
!! where \f$\Delta V\f$ is the integrated horizontal shear over the cloud
!! depth, \f$\Delta z\f$, (the ratio is converted to units of
!! \f$10^{-3} s^{-1}\f$). The variable "edt" is \f$1-E\f$ and is
!! constrained to the range \f$[0,0.9]\f$.
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 0.
        endif
      enddo
      do k = 2, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2
     &                  + (vo(i,k)-vo(i,k-1)) ** 2)
              vshear(i) = vshear(i) + shear
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
          e1=1.591-.639*vshear(i)
     &       +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
          edt(i)=1.-e1
          val =         .9
          edt(i) = min(edt(i),val)
          val =         .0
          edt(i) = max(edt(i),val)
        endif
      enddo
c
c--- what would the change be, that a cloud with unit mass
c--- will do to the environment?
c
!> ## Calculate the tendencies of the state variables (per unit cloud base mass flux) and the cloud base mass flux.
!> - Calculate the change in moist static energy, moisture mixing ratio,
!! and horizontal winds per unit cloud base mass flux for all layers
!! below cloud top from equations B.14 and B.15 from 
!! \cite grell_1993, and for the cloud top from B.16 and B.17.
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k <= kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
c
c--- changed due to subsidence and entrainment
c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < ktcon(i)) then
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
c
              dv1h = heo(i,k)
              dv2h = .5 * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
              dv1q = qo(i,k)
              dv2q = .5 * (qo(i,k) + qo(i,k-1))
              dv3q = qo(i,k-1)
c
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)
cj
              dellah(i,k) = dellah(i,k) +
     &     ( eta(i,k)*dv1h - eta(i,k-1)*dv3h
     &    -  tem*eta(i,k-1)*dv2h*dz
     &    +  tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz
     &         ) *g/dp
cj
              dellaq(i,k) = dellaq(i,k) +
     &     ( eta(i,k)*dv1q - eta(i,k-1)*dv3q
     &    -  tem*eta(i,k-1)*dv2q*dz
     &    +  tem1*eta(i,k-1)*.5*(qrcko(i,k)+qcko(i,k-1))*dz
     &         ) *g/dp
cj
              tem1=eta(i,k)*(uo(i,k)-ucko(i,k))
              tem2=eta(i,k-1)*(uo(i,k-1)-ucko(i,k-1))
              dellau(i,k) = dellau(i,k) + (tem1-tem2) * g/dp
cj
              tem1=eta(i,k)*(vo(i,k)-vcko(i,k))
              tem2=eta(i,k-1)*(vo(i,k-1)-vcko(i,k-1))
              dellav(i,k) = dellav(i,k) + (tem1-tem2) * g/dp
cj
            endif
          endif
        enddo
      enddo
c
c------- cloud top
c
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dv1h = heo(i,indx-1)
          dellah(i,indx) = eta(i,indx-1) *
     &                     (hcko(i,indx-1) - dv1h) * g / dp
          dv1q = qo(i,indx-1)
          dellaq(i,indx) = eta(i,indx-1) *
     &                     (qcko(i,indx-1) - dv1q) * g / dp
          dellau(i,indx) = eta(i,indx-1) *
     &             (ucko(i,indx-1) - uo(i,indx-1)) * g / dp
          dellav(i,indx) = eta(i,indx-1) *
     &             (vcko(i,indx-1) - vo(i,indx-1)) * g / dp
c
c  cloud water
c
          dellal(i,indx) = eta(i,indx-1) *
     &                     qlko_ktcon(i) * g / dp
        endif
      enddo
!
!
!  compute convective turn-over time
!
!> - Following \cite bechtold_et_al_2008,
!! calculate the convective turnover time using the mean updraft velocity
!! (wc) and the cloud depth. It is also proportional to the grid size (gdx).
      do i= 1, im
        if(cnvflg(i)) then
          tem = zi(i,ktcon1(i)) - zi(i,kbcon1(i))
          dtconv(i) = tem / wc(i)
          tfac = 1. + gdx(i) / 75000.
          dtconv(i) = tfac * dtconv(i)
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = max(dtconv(i),dt2)
          dtconv(i) = min(dtconv(i),dtmax)
        endif
      enddo
!
!> - Calculate advective time scale (tauadv) using a mean cloud layer wind speed.
      do i= 1, im
        if(cnvflg(i)) then
          sumx(i) = 0.
          umean(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if(cnvflg(i)) then
            if(k >= kbcon1(i) .and. k < ktcon1(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem = sqrt(u1(i,k)*u1(i,k)+v1(i,k)*v1(i,k))
              umean(i) = umean(i) + tem * dz
              sumx(i) = sumx(i) + dz
            endif
          endif
        enddo
      enddo
      do i= 1, im
        if(cnvflg(i)) then
           umean(i) = umean(i) / sumx(i)
           umean(i) = max(umean(i), 1.)
           tauadv(i) = gdx(i) / umean(i)
        endif
      enddo
c
c  compute cloud base mass flux as a function of the mean
c     updraft velcoity
c
!> - From \cite han_et_al_2017 equation 6, calculate cloud base mass flux as a function of the mean updraft velcoity.
!!  As discussed in \cite han_et_al_2017 , when dtconv
!! is larger than tauadv, the convective mixing is not fully conducted
!! before the cumulus cloud is advected out of the grid cell. In this
!! case, therefore, the cloud base mass flux is further reduced in
!! proportion to the ratio of tauadv to dtconv.
      do i= 1, im
        if(cnvflg(i)) then
          k = kbcon(i)
          rho = po(i,k)*100. / (rd*to(i,k))
          tfac = tauadv(i) / dtconv(i)
          tfac = min(tfac, 1.)
          xmb(i) = tfac*betaw*rho*wc(i)
        endif
      enddo
!
!> - For scale-aware parameterization, the updraft fraction (sigmagfm)
!! is first computed as a function of the lateral entrainment rate at
!! cloud base (see \cite han_et_al_2017 equation 4
!! and 5), following the study by 
!! \cite grell_and_freitas_2014.
      do i = 1, im
        if(cnvflg(i)) then
          tem = min(max(xlamue(i,kbcon(i)), 2.e-4), 6.e-4)
          tem = 0.2 / tem
          tem1 = 3.14 * tem * tem
          sigmagfm(i) = tem1 / garea(i)
          sigmagfm(i) = max(sigmagfm(i), 0.001)
          sigmagfm(i) = min(sigmagfm(i), 0.999)
        endif
      enddo
!
!> - Then, calculate the reduction factor (scaldfunc) of the vertical
!! convective eddy transport of mass flux as a function of updraft
!! fraction from the studies by \cite arakawa_and_wu_2013
!! (also see \cite han_et_al_2017 equation 1 and 2).
!! The final cloud base mass flux with scale-aware parameterization is
!! obtained from the mass flux when sigmagfm << 1, multiplied by the
!! reduction factor (\cite han_et_al_2017 equation 2).
      do i = 1, im
        if(cnvflg(i)) then
          if (gdx(i) < dxcrt) then
            scaldfunc(i) = (1.-sigmagfm(i)) * (1.-sigmagfm(i))
            scaldfunc(i) = max(min(scaldfunc(i), 1.0), 0.)
          else
            scaldfunc(i) = 1.0
          endif
          xmb(i) = xmb(i) * scaldfunc(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
!> ## For the "feedback control", calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
!! - Recalculate saturation specific humidity.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k <= kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
!> - Calculate the temperature tendency from the moist static energy and
!! specific humidity tendencies.
!> - Update the temperature, specific humidity, and horiztonal wind
!! state variables by multiplying the cloud base mass flux-normalized
!! tendencies by the cloud base mass flux.
!> - Accumulate column-integrated tendencies.
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        delubar(i) = 0.
        delvbar(i) = 0.
        qcond(i) = 0.
      enddo
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
!             tem = 1./rcs(i)
!             u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
!             v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2
              dp = 1000. * del(i,k)
              delhbar(i) = delhbar(i) + dellah(i,k)*xmb(i)*dp/g
              delqbar(i) = delqbar(i) + dellaq(i,k)*xmb(i)*dp/g
              deltbar(i) = deltbar(i) + dellat*xmb(i)*dp/g
              delubar(i) = delubar(i) + dellau(i,k)*xmb(i)*dp/g
              delvbar(i) = delvbar(i) + dellav(i,k)*xmb(i)*dp/g
            endif
          endif
        enddo
      enddo
!
!> - Recalculate saturation specific humidity using the updated temperature.
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
c
!> - Add up column-integrated convective precipitation by multiplying
!! the normalized value by the cloud base mass flux.
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k < ktcon(i) .and. k > kb(i)) then
              rntot(i) = rntot(i) + pwo(i,k) * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
c
c evaporating rain
c
!> - Determine the evaporation of the convective precipitation and
!! update the integrated convective precipitation.
!> - Update state temperature and moisture to account for evaporation
!! of convective precipitation.
!> - Update column-integrated tendencies to account for evaporation of
!! convective precipitation.
      do k = km, 1, -1
        do i = 1, im
          if (k <= kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i)) then
              if(k < ktcon(i) .and. k > kb(i)) then
                rn(i) = rn(i) + pwo(i,k) * xmb(i) * .001 * dt2
              endif
            endif
            if(flg(i) .and. k < ktcon(i)) then
              evef = edt(i) * evfact
              if(islimsk(i) == 1) evef=edt(i) * evfactl
!             if(islimsk(i) == 1) evef=.07
c             if(islimsk(i) == 1) evef = 0.
              qcond(i) = evef * (q1(i,k) - qeso(i,k))
     &                 / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              if(rn(i) > 0. .and. qcond(i) < 0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*g/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * dp / g
              endif
              if(rn(i) > 0. .and. qcond(i) < 0. .and.
     &           delq2(i) > rntot(i)) then
                qevap(i) = 1000.* g * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i) > 0. .and. qevap(i) > 0.) then
                tem  = .001 * dp / g
                tem1 = qevap(i) * tem
                if(tem1 > rn(i)) then
                  qevap(i) = rn(i) / tem
                  rn(i) = 0.
                else
                  rn(i) = rn(i) - tem1
                endif
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + .001*dp*qevap(i)/g
              endif
              delqbar(i) = delqbar(i) + delq(i)*dp/g
              deltbar(i) = deltbar(i) + deltv(i)*dp/g
            endif
          endif
        enddo
      enddo
cj
!     do i = 1, im
!     if(me == 31 .and. cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' shallow delhbar, delqbar, deltbar = ',
!    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!       print *, ' shallow delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip =', hvap*rn(i)*1000./dt2
!       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!     endif
!     enddo
cj
      do i = 1, im
        if(cnvflg(i)) then
          if(rn(i) < 0. .or. .not.flg(i)) rn(i) = 0.
          ktop(i) = ktcon(i)
          kbot(i) = kbcon(i)
          kcnv(i) = 2
        endif
      enddo
c
c  convective cloud water
c
!> - Calculate shallow convective cloud water.
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if (k >= kbcon(i) .and. k < ktcon(i)) then
              cnvw(i,k) = cnvwt(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo

c
c  convective cloud cover
c
!> - Calculate convective cloud cover, which is used when pdf-based cloud
!! fraction is used (i.e., pdfcld=.true.).
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if (k >= kbcon(i) .and. k < ktcon(i)) then
              cnvc(i,k) = 0.04 * log(1. + 675. * eta(i,k) * xmb(i))
              cnvc(i,k) = min(cnvc(i,k), 0.2)
              cnvc(i,k) = max(cnvc(i,k), 0.0)
            endif
          endif
        enddo
      enddo

c
c  cloud water
c
!> - Separate detrained cloud water into liquid and ice species as a
!! function of temperature only.
      if (ncloud > 0) then
!
      do k = 1, km1
        do i = 1, im
          if (cnvflg(i)) then
!           if (k > kb(i) .and. k <= ktcon(i)) then
            if (k >= kbcon(i) .and. k <= ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
              tem1 = max(0.0, min(1.0, (tcr-t1(i,k))*tcrf))
              if (ql2(i,k) > -999.0) then
                ql1(i,k) = ql1(i,k) + tem * tem1            ! ice
                ql2(i,k) = ql2(i,k) + tem *(1.0-tem1)       ! water
              else
                ql1(i,k) = ql1(i,k) + tem
              endif
            endif
          endif
        enddo
      enddo
!
      endif
!
! hchuang code change
!
!> - Calculate and retain the updraft mass flux for dust transport by
!! cumulus convection.
!
!> - Calculate the updraft convective mass flux.
      do k = 1, km
        do i = 1, im
          if(cnvflg(i)) then
            if(k >= kb(i) .and. k < ktop(i)) then
              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo
!> - save the updraft convective mass flux at cloud top.
      do i = 1, im
        if(cnvflg(i)) then
           k = ktop(i)-1
           dt_mf(i,k) = ud_mf(i,k)
        endif
      enddo
!!
      return
      end subroutine sasas_shal_run
!! @}
!! @}

!!
!! \section arg_table_sasas_shal_init Argument Table
!!
      subroutine sasas_shal_finalize
      end subroutine sasas_shal_finalize

      end module sasas_shal

