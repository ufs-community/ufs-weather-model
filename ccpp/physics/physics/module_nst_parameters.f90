!>\file module_nst_parameters.f90
!>\defgroup nst_parameters GFS NSST Parameter Module
!! \ingroup GFS_NSST
!!\brief This module includes constants and parameters used in GFS
!! near sea surface temperature scheme.
!! @{
module module_nst_parameters
  use machine, only :  kind_phys &
       ,kind_rad ! for astronomy (date) calculations
  !
  ! air constants and coefficients from the atmospehric model
  use physcons, only: &
       eps =>  con_eps & 
       ,cp_a => con_cp &          ! spec heat air @p    (j/kg/k)
       , epsm1 => con_epsm1 & 
       , hvap => con_hvap &       ! lat heat h2o cond   (j/kg)
       ,sigma_r => con_sbc  &     ! stefan-boltzmann    (w/m2/k4)
       ,grav => con_g         &   ! acceleration due to gravity (kg/m/s^2)
       ,omega => con_omega    &    ! ang vel of earth    (1/s)
       ,rvrdm1 => con_fvirt &
       ,rd => con_rd &
       ,rocp => con_rocp  &        ! r/cp          
       ,pi => con_pi
  !
  ! note: take timestep from here later
  public 
  integer :: &
       niter_conv = 5, &
       niter_z_w  = 5, &
       niter_sfs  = 5
  real (kind=kind_phys), parameter :: & 
       !
       ! general constants
        sec_in_day=86400.       &
       ,sec_in_hour=3600.       &
       ,solar_time_6am=21600.0  &
       ,const_rot=0.000073      &          !< constant to calculate corioli force
       ,ri_c=0.65               & 
       ,ri_g=0.25               & 
       ,eps_z_w=0.01            &          !< criteria to finish iterations for z_w
       ,eps_conv=0.01           &          !< criteria to finish iterations for d_conv
       ,eps_sfs=0.01            &          !< criteria to finish iterations for d_sfs
       ,z_w_max=30.0            &          !< max warm layer thickness
       ,z_w_min=0.2             &          !< min warm layer thickness
       ,z_w_ini=0.2             &          !< initial warm layer thickness in dtl_onset
       ,z_c_max=0.01            &          !< maximum of sub-layer thickness (m)
       ,z_c_ini=0.001           &          !< initial value of z_c
       ,ustar_a_min=0.031       &          !< minimum of friction wind speed (m/s): 0.031 ~ 1m/s at 10 m hight
       ,tau_min=0.005           &          !< minimum of wind stress for dtm
       ,exp_const=9.5           &          !< coefficient in exponet profile
       ,delz=0.1                &          !< vertical increment for integral calculation   (m)
       ,von=0.4                 &          !< von karman's "constant"      
       ,t0k=273.16              &          !<  celsius to kelvin
       ,gray=0.97               &
       ,sst_max=308.16          &
       ,tw_max=5.0              &
       ,wd_max=2.0              &
       ,omg_m =1.0              &          !< trace factor to apply salinity effect
       ,omg_rot = 1.0           &          !< trace factor to apply rotation effect
       ,omg_sh = 1.0            &          !< trace factor to apply sensible heat due to rainfall effect
       ,visw=1.e-6 &                       !< m2/s kinematic viscosity water
       ,novalue=0 &
       ,smallnumber=1.e-6 & 
       ,timestep_oc=sec_in_day/8. &        !< time step in the ocean model (3 hours)
       ,radian=2.*pi/180.       & 
       ,rad2deg=180./pi       & 
       ,cp_w=4000.   &                     !< specific heat water (j/kg/k )
       ,rho0_w=1022.0 &                    !< density water (kg/m3 ) (or 1024.438)
       ,vis_w=1.e-6  &                     !< kinematic viscosity water (m2/s )
       ,tc_w=0.6    &                      !< thermal conductivity water (w/m/k )
       ,capa_w =3950.0 &                   !< heat capacity of sea water      !
       ,thref =1.0e-3                      !< reference value of specific volume (m**3/kg) 

!!$!============================================
!!$
!!$  ,lvapor=2.453e6 &        ! latent heat of vaporization note: make it function of t ????? note the same as hvap        
!!$       ,alpha=1 ! thermal expansion coefficient
!!$  ,beta ! saline contraction coefficient
!!$  ,cp=1 !=1 specific heat of sea water
!!$  ,g=1 ! acceleration due to gravity
!!$  ,kw=1 ! thermal conductivity of water
!!$  ,nu=1    !kinematic wiscosity
!!$  ,rho_w=1 !water density
!!$  ,rho_a=1 !air density
!!$  ,l_vapr=2.453e6
!!$  ,novalue=--1.0e+10
!!$
!!$c factors
!!$      beta=1.2     !given as 1.25 in fairall et al.(1996)
!!$      von=0.4      ! von karman's "constant"
!!$c      fdg=1.00     ! fairall's lkb rr to von karman adjustment
!!$      fdg=1.00     !based on results from flux workshop august 1995
!!$      tok=273.16   ! celsius to kelvin
!!$      twopi=3.14159*2.
!!$ 
!!$c air constants and coefficients
!!$      rgas=287.1                  !j/kg/k     gas const. dry air
!!$      xlv=(2.501-0.00237*ts)*1e+6  !j/kg  latent heat of vaporization at ts
!!$      cpa=1004.67                 !j/kg/k specific heat of dry air (businger 1982)
!!$      cpv=cpa*(1+0.84*q)          !moist air - currently not used (businger 1982)
!!$      rhoa=p*100./(rgas*(t+tok)*(1.+.61*q)) !kg/m3  moist air density ( " )
!!$      visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t)   !m2/s
!!$          !kinematic viscosity of dry air - andreas (1989) crrel rep. 89-11
!!$c 
!!$c cool skin constants
!!$      al=2.1e-5*(ts+3.2)**0.79     !water thermal expansion coefft.
!!$      be=0.026                     !salinity expansion coefft.
!!$      cpw=4000.                    !j/kg/k specific heat water
!!$      rhow=1022.                   !kg/m3  density water
!!$      visw=1.e-6                   !m2/s kinematic viscosity water
!!$      tcw=0.6                      !w/m/k   thermal conductivity water
!!$      bigc=16.*grav*cpw*(rhow*visw)**3/(tcw*tcw*rhoa*rhoa)
!!$      wetc=0.622*xlv*qs/(rgas*(ts+tok)**2) !correction for dq;slope of sat. vap.
!!$
!!$!
!!$! functions
!!$
!!$
!!$  real,    parameter :: timestep=86400.    !integration time step, second
!!$
!!$  real, parameter    :: grav =9.81         !gravity, kg/m/s^2
!!$  real, parameter    :: capa =3950.0       !heat capacity of sea water
!!$  real, parameter    :: rhoref = 1024.438  !sea water reference density, kg/m^3
!!$  real   , parameter :: hslab=50.0         !slab ocean depth
!!$  real   , parameter :: bad=-1.0e+10
!!$  real   , parameter :: tmin=2.68e+02 
!!$  real   , parameter :: tmax=3.11e+02
!!$
!!$  real, parameter :: grav =9.81           !gravity, kg/m/s^2
!!$  real, parameter :: capa =3950.0         !heat capacity of sea water 
!!$  real, parameter :: rhoref = 1024.438    !sea water reference density, kg/m^3
!!$  real, parameter :: tmin=2.68e+02        !normal minimal temp
!!$  real, parameter :: tmax=3.11e+02        !normal max temp
!!$  real, parameter :: smin=1.0             !normal minimal salt
!!$  real, parameter :: smax=50.             !normal maximum salt
!!$  real, parameter :: visct=1.e-5          !viscocity for temperature diffusion
!!$  real, parameter :: viscs=1.e-5          !viscocity for salt diffusion
!!$
!!$
end module module_nst_parameters
!! @}
