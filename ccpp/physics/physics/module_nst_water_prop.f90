!>\file module_nst_water_prop.f90
!>\defgroup nst_water_prop GFS NSST Water Prop Module
!! @{
!>\ingroup GFS_NSST
!!\brief
module module_nst_water_prop
  use machine, only : kind_phys
  use module_nst_parameters, only : t0k
  !
  private
  public :: rhocoef,density,sw_rad,sw_rad_aw,sw_rad_sum,sw_rad_upper,sw_rad_upper_aw,sw_rad_skin,grv,solar_time_from_julian,compjd, &
            sw_ps_9b,sw_ps_9b_aw,get_dtzm_point,get_dtzm_2d
      
  !
  interface sw_ps_9b
     module procedure sw_ps_9b
  end interface
  interface sw_ps_9b_aw
     module procedure sw_ps_9b_aw
  end interface
  !
  interface sw_rad
     module procedure sw_fairall_6exp_v1  ! sw_wick_v1
  end interface
  interface sw_rad_aw
     module procedure sw_fairall_6exp_v1_aw
  end interface
  interface sw_rad_sum
     module procedure sw_fairall_6exp_v1_sum
  end interface
  interface sw_rad_upper
     module procedure sw_soloviev_3exp_v2
  end interface
  interface sw_rad_upper_aw
     module procedure sw_soloviev_3exp_v2_aw
  end interface
  interface sw_rad_skin
     module procedure sw_ohlmann_v1
  end interface
contains
  ! ------------------------------------------------------
!> This subroutine computes thermal expansion coefficient (alpha)
!! and saline contraction coefficient (beta). 
  subroutine rhocoef(t, s, rhoref, alpha, beta)
    ! ------------------------------------------------------

    !  compute thermal expansion coefficient (alpha) 
    !  and saline contraction coefficient (beta) using 
    !  the international equation of state of sea water 
    !  (1980). ref: pond and pickard, introduction to 
    !  dynamical oceanography, pp310.  
    !  note: compression effects are not included

    implicit none
    real(kind=kind_phys), intent(in)  :: t, s, rhoref 
    real(kind=kind_phys), intent(out) :: alpha, beta  
    real(kind=kind_phys) :: tc

    tc = t - t0k

    alpha =                                                        & 
         6.793952e-2                                              & 
         - 2.0 * 9.095290e-3 * tc     +  3.0 * 1.001685e-4 * tc**2  & 
         - 4.0 * 1.120083e-6 * tc**3  +  5.0 * 6.536332e-9 * tc**4  & 
         - 4.0899e-3 * s                                            & 
         + 2.0 * 7.6438e-5 * tc * s  -  3.0 * 8.2467e-7 * tc**2 * s & 
         + 4.0 * 5.3875e-9 * tc**3 * s                              & 
         + 1.0227e-4 * s**1.5 -  2.0 * 1.6546e-6 * tc * s**1.5

    ! note: rhoref - specify 
    !
    alpha =  -alpha/rhoref

    beta  =                                             &
         8.24493e-1          -  4.0899e-3 * tc           &
         + 7.6438e-5 * tc**2 -  8.2467e-7 * tc**3        &
         + 5.3875e-9 * tc**4 -  1.5 * 5.72466e-3 * s**.5 &
         + 1.5 * 1.0227e-4 * tc * s**.5                  &
         -  1.5 * 1.6546e-6 * tc**2 * s**.5              &
         + 2.0 * 4.8314e-4 * s

    beta = beta / rhoref

  end subroutine rhocoef
  ! ----------------------------------------
!> This subroutine computes sea water density.
  subroutine density(t, s, rho)
    ! ----------------------------------------
    implicit none

    ! input
    real(kind=kind_phys), intent(in)  :: t     !unit, k
    real(kind=kind_phys), intent(in)  :: s     !unit, 1/1000
    ! output
    real(kind=kind_phys), intent(out) :: rho   !unit, kg/m^3 
    ! local
    real(kind=kind_phys) :: tc

    ! compute density using the international equation 
    ! of state of sea water 1980, (pond and pickard, 
    ! introduction to dynamical oceanography, pp310). 
    ! compression effects are not included

    rho = 0.0
    tc = t - t0k

    !  effect of temperature on density (lines 1-3)
    !  effect of temperature and salinity on density (lines 4-8)
    rho = &
         999.842594                 +  6.793952e-2 * tc     &
         - 9.095290e-3 * tc**2        +  1.001685e-4 * tc**3     &
         - 1.120083e-6 * tc**4        +  6.536332e-9 * tc**5     &
         + 8.24493e-1 * s          -  4.0899e-3 * tc * s         &
         + 7.6438e-5 * tc**2 * s   -  8.2467e-7 * tc**3 * s      &
         + 5.3875e-9 * tc**4 * s   -  5.72466e-3 * s**1.5        &
         + 1.0227e-4 * tc * s**1.5 -  1.6546e-6 * tc**2 * s**1.5 &
         + 4.8314e-4 * s**2

  end subroutine density
  !
  !======================
  !
!> This subroutine computes the fraction of the solar radiation absorbed 
!! by the depth z following Paulson and Simpson (1981) \cite paulson_and_simpson_1981 .
  elemental subroutine sw_ps_9b(z,fxp)
    !
    ! fraction of the solar radiation absorbed by the ocean at the depth z 
    ! following paulson and simpson, 1981
    !
    ! input:
    ! z:       depth (m)
    !
    ! output:
    ! fxp: fraction of the solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real,intent(in):: z
    real,intent(out):: fxp
    real, dimension(9), parameter :: f=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/) &
                                ,gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    !
    if(z>0) then
      fxp=1.0-(f(1)*exp(-z/gamma(1))+f(2)*exp(-z/gamma(2))+f(3)*exp(-z/gamma(3))+ &
               f(4)*exp(-z/gamma(4))+f(5)*exp(-z/gamma(5))+f(6)*exp(-z/gamma(6))+ &
               f(7)*exp(-z/gamma(7))+f(8)*exp(-z/gamma(8))+f(9)*exp(-z/gamma(9)))
    else
       fxp=0.
    endif
    !
  end subroutine sw_ps_9b
  !
  !======================
  !
  !
  !======================
  !
  elemental subroutine sw_ps_9b_aw(z,aw)
    !
    ! d(fw)/d(z) for 9-band 
    !
    ! input:
    ! z:       depth (m)
    !
    ! output:
    ! fxp: fraction of the solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real,intent(in):: z
    real,intent(out):: aw
    real, dimension(9), parameter :: f=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/) &
                                ,gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    !
    if(z>0) then
      aw=(f(1)/gamma(1))*exp(-z/gamma(1))+(f(2)/gamma(2))*exp(-z/gamma(2))+(f(3)/gamma(3))*exp(-z/gamma(3))+ &
         (f(1)/gamma(4))*exp(-z/gamma(4))+(f(2)/gamma(5))*exp(-z/gamma(5))+(f(6)/gamma(6))*exp(-z/gamma(6))+ &
         (f(1)/gamma(7))*exp(-z/gamma(7))+(f(2)/gamma(8))*exp(-z/gamma(8))+(f(9)/gamma(9))*exp(-z/gamma(9))
    else
       aw=0.
    endif
    !
  end subroutine sw_ps_9b_aw
  !
  !======================
  elemental subroutine sw_fairall_6exp_v1(z,fxp)
    !
    ! fraction of the solar radiation absorbed by the ocean at the depth z (fairall et all, 1996, p. 1298)
    ! following paulson and simpson, 1981
    !
    ! input:
    ! z:       depth (m)
    !
    ! output:
    ! fxp: fraction of the solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real(kind=kind_phys),intent(in):: z
    real(kind=kind_phys),intent(out):: fxp
    real(kind=kind_phys), dimension(9), parameter :: f=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/) &
         ,gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    real(kind=kind_phys),dimension(9) :: zgamma
    real(kind=kind_phys),dimension(9) :: f_c
    !
    if(z>0) then
       zgamma=z/gamma
       f_c=f*(1.-1./zgamma*(1-exp(-zgamma)))
       fxp=sum(f_c)
    else
       fxp=0.
    endif
    !
  end subroutine sw_fairall_6exp_v1
  !
  !======================
  !
  !
  elemental subroutine sw_fairall_6exp_v1_aw(z,aw)
    !
    ! fraction of the solar radiation absorbed by the ocean at the depth z (fairall et all, 1996, p. 1298)
    ! following paulson and simpson, 1981
    !
    ! input:
    ! z:       depth (m)
    !
    ! output:
    ! aw: d(fxp)/d(z)
    !
    ! fxp: fraction of the solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real(kind=kind_phys),intent(in):: z
    real(kind=kind_phys),intent(out):: aw
    real(kind=kind_phys) :: fxp
    real(kind=kind_phys), dimension(9), parameter :: f=(/0.237,0.36,0.179,0.087,0.08,0.0246,0.025,0.007,0.0004/) &
         ,gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    real(kind=kind_phys),dimension(9) :: zgamma
    real(kind=kind_phys),dimension(9) :: f_aw
    !
    if(z>0) then
       zgamma=z/gamma
       f_aw=(f/z)*((gamma/z)*(1-exp(-zgamma))-exp(-zgamma))
       aw=sum(f_aw)

!      write(*,'(a,f6.2,f12.6,9f10.4)') 'z,aw in sw_rad_aw : ',z,aw,f_aw

    else
       aw=0.
    endif
    !
  end subroutine sw_fairall_6exp_v1_aw
  !
  elemental subroutine sw_fairall_6exp_v1_sum(z,sum)
    !
    ! fraction of the solar radiation absorbed by the ocean at the depth z (fairall et all, 1996, p. 1298)
    ! following paulson and simpson, 1981
    !
    ! input:
    ! z:       depth (m)
    !
    ! output:
    ! sum: for convection depth calculation
    !
    !
    implicit none
    real(kind=kind_phys),intent(in):: z
    real(kind=kind_phys),intent(out):: sum
    real(kind=kind_phys), dimension(9), parameter :: gamma=(/34.8,2.27,3.15e-2,5.48e-3,8.32e-4,1.26e-4,3.13e-4,7.82e-5,1.44e-5/)
    real(kind=kind_phys),dimension(9) :: zgamma
    real(kind=kind_phys),dimension(9) :: f_sum
    !
!    zgamma=z/gamma
!    f_sum=(zgamma/z)*exp(-zgamma)
!    sum=sum(f_sum)

    sum=(1.0/gamma(1))*exp(-z/gamma(1))+(1.0/gamma(2))*exp(-z/gamma(2))+(1.0/gamma(3))*exp(-z/gamma(3))+ &
        (1.0/gamma(4))*exp(-z/gamma(4))+(1.0/gamma(5))*exp(-z/gamma(5))+(1.0/gamma(6))*exp(-z/gamma(6))+ &
        (1.0/gamma(7))*exp(-z/gamma(7))+(1.0/gamma(8))*exp(-z/gamma(8))+(1.0/gamma(9))*exp(-z/gamma(9))
    !
  end subroutine sw_fairall_6exp_v1_sum
  !
  !======================

  elemental subroutine sw_fairall_simple_v1(f_sol_0,z,df_sol_z)
    !
    ! solar radiation absorbed by the ocean at the depth z (fairall et all, 1996, p. 1298)
    ! 
    ! input: 
    ! f_sol_0: solar radiation at the ocean surface (w/m^2)
    ! z:       depth (m)
    !
    ! output:
    ! df_sol_z: solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real(kind=kind_phys),intent(in):: z,f_sol_0
    real(kind=kind_phys),intent(out):: df_sol_z
    !
    if(z>0) then
       df_sol_z=f_sol_0*(0.137+11.0*z-6.6e-6/z*(1.-exp(-z/8.e-4)))
    else
       df_sol_z=0.
    endif
    !
  end subroutine sw_fairall_simple_v1
  !
  !======================
  !
  elemental subroutine sw_wick_v1(f_sol_0,z,df_sol_z)
    !
    ! solar radiation absorbed by the ocean at the depth z (zeng and beljaars, 2005, p.5)
    ! 
    ! input: 
    ! f_sol_0: solar radiation at the ocean surface (w/m^2)
    ! z:       depth (m)
    !
    ! output:
    ! df_sol_z: solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real(kind=kind_phys),intent(in):: z,f_sol_0
    real(kind=kind_phys),intent(out):: df_sol_z
    !
    if(z>0) then
       df_sol_z=f_sol_0*(0.065+11.0*z-6.6e-5/z*(1.-exp(-z/8.e-4)))
    else
       df_sol_z=0.
    endif
    !
  end subroutine sw_wick_v1
  !
  !======================
  !
  elemental subroutine sw_soloviev_3exp_v1(f_sol_0,z,df_sol_z)
    !
    ! solar radiation absorbed by the ocean at the depth z (fairall et all, 1996, p. 1301)
    ! following soloviev, 1982
    ! 
    ! input: 
    ! f_sol_0: solar radiation at the ocean surface (w/m^2)
    ! z:       depth (m)
    !
    ! output:
    ! df_sol_z: solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real(kind=kind_phys),intent(in):: z,f_sol_0
    real(kind=kind_phys),intent(out):: df_sol_z
    real(kind=kind_phys),dimension(3) :: f_c
    real(kind=kind_phys), dimension(3), parameter :: f=(/0.45,0.27,0.28/) &
         ,gamma=(/12.8,0.357,0.014/)
    !
    if(z>0) then
       f_c      = f*gamma(int(1-exp(-z/gamma)))
       df_sol_z = f_sol_0*(1.0-sum(f_c)/z)
    else
       df_sol_z = 0.
    endif
    !
  end subroutine sw_soloviev_3exp_v1
  !
  !======================
  !
  elemental subroutine sw_soloviev_3exp_v2(f_sol_0,z,df_sol_z)
    !
    ! solar radiation absorbed by the ocean at the depth z (fairall et all, 1996, p. 1301)
    ! following soloviev, 1982
    ! 
    ! input: 
    ! f_sol_0: solar radiation at the ocean surface (w/m^2)
    ! z:       depth (m)
    !
    ! output:
    ! df_sol_z: solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real(kind=kind_phys),intent(in):: z,f_sol_0
    real(kind=kind_phys),intent(out):: df_sol_z
    !
    if(z>0) then
       df_sol_z=f_sol_0*(1.0 &
            -(0.28*0.014*(1.-exp(-z/0.014)) &
            +0.27*0.357*(1.-exp(-z/0.357)) &        
            +.45*12.82*(1.-exp(-z/12.82)))/z &
            )
    else
       df_sol_z=0.
    endif
    !
  end subroutine sw_soloviev_3exp_v2

  elemental subroutine sw_soloviev_3exp_v2_aw(z,aw)
    !
    ! aw = d(fxp)/d(z)
    ! following soloviev, 1982
    !
    ! input:
    ! z:       depth (m)
    !
    ! output:
    ! aw: d(fxp)/d(z)
    !
    implicit none
    real(kind=kind_phys),intent(in):: z
    real(kind=kind_phys),intent(out):: aw
    real(kind=kind_phys):: fxp
    !
    if(z>0) then
       fxp=(1.0 &
            -(0.28*0.014*(1.-exp(-z/0.014)) &
            + 0.27*0.357*(1.-exp(-z/0.357)) &
            + 0.45*12.82*(1.-exp(-z/12.82)))/z &
            )
       aw=1.0-fxp-(0.28*exp(-z/0.014)+0.27*exp(-z/0.357)+0.45*exp(-z/12.82))
    else
       aw=0.
    endif
  end subroutine sw_soloviev_3exp_v2_aw
  !
  !
  !======================
  !
  elemental subroutine sw_ohlmann_v1(z,fxp)
    !
    ! fraction of the solar radiation absorbed by the ocean at the depth z
    !
    ! input:
    ! z:       depth (m)
    !
    ! output:
    ! fxp: fraction of the solar radiation absorbed by the ocean at depth z (w/m^2)
    !
    implicit none
    real(kind=kind_phys),intent(in):: z
    real(kind=kind_phys),intent(out):: fxp
    !
    if(z>0) then
       fxp=.065+11.*z-6.6e-5/z*(1.-exp(-z/8.0e-4))
    else
       fxp=0.
    endif
    !
  end subroutine sw_ohlmann_v1
  !

function grv(lat)
  real(kind=kind_phys) :: lat
  real(kind=kind_phys) :: gamma,c1,c2,c3,c4,pi,phi,x
  gamma=9.7803267715
  c1=0.0052790414
  c2=0.0000232718
  c3=0.0000001262
  c4=0.0000000007
  pi=3.141593
                                                                                                                                                             
  phi=lat*pi/180
  x=sin(phi)
  grv=gamma*(1+(c1*x**2)+(c2*x**4)+(c3*x**6)+(c4*x**8))
  !print *,'grav=',grv,lat
end function grv

!>This subroutine computes solar time from the julian date.
subroutine solar_time_from_julian(jday,xlon,soltim)
  !
  ! calculate solar time from the julian date
  !
  implicit none
  real(kind=kind_phys), intent(in)  :: jday
  real(kind=kind_phys), intent(in)  :: xlon
  real(kind=kind_phys), intent(out) :: soltim
  real(kind=kind_phys)                            :: fjd,xhr,xmin,xsec,intime
  integer                                        :: nn
  !
  fjd=jday-floor(jday)
  fjd=jday
  xhr=floor(fjd*24.0)-sign(12.0,fjd-0.5)
  xmin=nint(fjd*1440.0)-(xhr+sign(12.0,fjd-0.5))*60
  xsec=0
  intime=xhr+xmin/60.0+xsec/3600.0+24.0
  soltim=mod(xlon/15.0+intime,24.0)*3600.0
end subroutine solar_time_from_julian

!
!***********************************************************************
!
!> This subroutine computes julian day and fraction from year,
!! month, day and time UTC.
      subroutine compjd(jyr,jmnth,jday,jhr,jmn,jd,fjd)
!fpp$ noconcur r
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    compjd      computes julian day and fraction
!   prgmmr: kenneth campana  org: w/nmc23    date: 89-07-07
!
! abstract: computes julian day and fraction
!   from year, month, day and time utc.
!
! program history log:
!   77-05-06  ray orzol,gfdl
!   98-05-15  iredell   y2k compliance
!
! usage:    call compjd(jyr,jmnth,jday,jhr,jmn,jd,fjd)
!   input argument list:
!     jyr      - year (4 digits)
!     jmnth    - month
!     jday     - day
!     jhr      - hour
!     jmn      - minutes 
!   output argument list:
!     jd       - julian day.
!     fjd      - fraction of the julian day.
!
! subprograms called:
!   iw3jdn     compute julian day number
!
! attributes:
!   language: fortran.
!
!$$$
      use machine , only :kind_phys
      implicit none
!
      integer jyr,jmnth,jday,jhr,jmn,jd
      integer iw3jdn
      real (kind=kind_phys) fjd
      jd=iw3jdn(jyr,jmnth,jday)
      if(jhr.lt.12) then
        jd=jd-1
        fjd=0.5+jhr/24.+jmn/1440.
      else
        fjd=(jhr-12)/24.+jmn/1440.
      endif
      end subroutine compjd

!>This subroutine computes dtm (the mean of \f$dT(z)\f$). 
 subroutine get_dtzm_point(xt,xz,dt_cool,zc,z1,z2,dtm)
! ===================================================================== !
!                                                                       !
!  description:  get dtm = mean of dT(z) (z1 - z2) with NSST dT(z)      !
!                dT(z) = (1-z/xz)*dt_warm - (1-z/zc)*dt_cool            !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call get_dtm12                                                     !
!                                                                       !
!       inputs:                                                         !
!          (xt,xz,dt_cool,zc,z1,z2,                                     !
!       outputs:                                                        !
!          dtm)                                                         !
!                                                                       !
!  program history log:                                                 !
!                                                                       !
!         2015  -- xu li       createad original code                   !
!  inputs:                                                              !
!     xt      - real, heat content in dtl                            1  !
!     xz      - real, dtl thickness                                  1  !
!     dt_cool - real, sub-layer cooling amount                       1  !
!     zc      - sub-layer cooling thickness                          1  !
!     z1      - lower bound of depth of sea temperature              1  !
!     z2      - upper bound of depth of sea temperature              1  !
!  outputs:                                                             !
!     dtm   - mean of dT(z)  (z1 to z2)                              1  !
!
  use machine , only : kind_phys

  implicit none

  real (kind=kind_phys), intent(in)  :: xt,xz,dt_cool,zc,z1,z2
  real (kind=kind_phys), intent(out) :: dtm
! Local variables
  real (kind=kind_phys) :: dt_warm,dtw,dtc

!
! get the mean warming in the range of z=z1 to z=z2
!
  dtw = 0.0
  if ( xt > 0.0 ) then
    dt_warm = (xt+xt)/xz      ! Tw(0)
    if ( z1 < z2) then
      if ( z2 < xz ) then
        dtw = dt_warm*(1.0-(z1+z2)/(xz+xz))
      elseif ( z1 < xz .and. z2 >= xz ) then
        dtw = 0.5*(1.0-z1/xz)*dt_warm*(xz-z1)/(z2-z1)
      endif
    elseif ( z1 == z2 ) then
      if ( z1 < xz ) then
        dtw = dt_warm*(1.0-z1/xz)
      endif
    endif
  endif
!
! get the mean cooling in the range of z=z1 to z=z2
!
  dtc = 0.0
  if ( zc > 0.0 ) then
    if ( z1 < z2) then
      if ( z2 < zc ) then
        dtc = dt_cool*(1.0-(z1+z2)/(zc+zc))
      elseif ( z1 < zc .and. z2 >= zc ) then
        dtc = 0.5*(1.0-z1/zc)*dt_cool*(zc-z1)/(z2-z1)
      endif
    elseif ( z1 == z2 ) then
      if ( z1 < zc ) then
        dtc = dt_cool*(1.0-z1/zc)
      endif
    endif
  endif

!
! get the mean T departure from Tf in the range of z=z1 to z=z2
!
  dtm = dtw - dtc

 end subroutine get_dtzm_point

 subroutine get_dtzm_2d(xt,xz,dt_cool,zc,slmsk,z1,z2,nx,ny,dtm)
! ===================================================================== !
!                                                                       !
!  description:  get dtm = mean of dT(z) (z1 - z2) with NSST dT(z)      !
!                dT(z) = (1-z/xz)*dt_warm - (1-z/zc)*dt_cool            !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call get_dtzm_2d                                                     !
!                                                                       !
!       inputs:                                                         !
!          (xt,xz,dt_cool,zc,z1,z2,                                     !
!       outputs:                                                        !
!          dtm)                                                         !
!                                                                       !
!  program history log:                                                 !
!                                                                       !
!         2015  -- xu li       createad original code                   !
!  inputs:                                                              !
!     xt      - real, heat content in dtl                            1  !
!     xz      - real, dtl thickness                                  1  !
!     dt_cool - real, sub-layer cooling amount                       1  !
!     zc      - sub-layer cooling thickness                          1  !
!     nx      - integer, dimension in x-direction (zonal)            1  !
!     ny      - integer, dimension in y-direction (meridional)       1  !
!     z1      - lower bound of depth of sea temperature              1  !
!     z2      - upper bound of depth of sea temperature              1  !
!  outputs:                                                             !
!     dtm   - mean of dT(z)  (z1 to z2)                              1  !
!
  use machine , only : kind_phys

  implicit none

  integer, intent(in) :: nx,ny
  real (kind=kind_phys), dimension(nx,ny), intent(in)  :: xt,xz,dt_cool,zc,slmsk
  real (kind=kind_phys), intent(in)  :: z1,z2
  real (kind=kind_phys), dimension(nx,ny), intent(out) :: dtm                    
! Local variables
  integer :: i,j
  real (kind=kind_phys), dimension(nx,ny) :: dtw,dtc
  real (kind=kind_phys) :: dt_warm


!$omp parallel do private(j,i)
  do j = 1, ny
    do i= 1, nx
!
!     initialize dtw & dtc as zeros
!
      dtw(i,j) = 0.0      
      dtc(i,j) = 0.0      
      if ( slmsk(i,j) == 0.0 ) then
!
!       get the mean warming in the range of z=z1 to z=z2
!
        if ( xt(i,j) > 0.0 ) then
          dt_warm = (xt(i,j)+xt(i,j))/xz(i,j)      ! Tw(0)
          if ( z1 < z2) then
            if ( z2 < xz(i,j) ) then
              dtw(i,j) = dt_warm*(1.0-(z1+z2)/(xz(i,j)+xz(i,j)))
            elseif ( z1 < xz(i,j) .and. z2 >= xz(i,j) ) then
              dtw(i,j) = 0.5*(1.0-z1/xz(i,j))*dt_warm*(xz(i,j)-z1)/(z2-z1)
            endif
          elseif ( z1 == z2 ) then
            if ( z1 < xz(i,j) ) then
              dtw(i,j) = dt_warm*(1.0-z1/xz(i,j))
            endif
          endif
        endif
!
!       get the mean cooling in the range of z=0 to z=zsea
!
        if ( zc(i,j) > 0.0 ) then
          if ( z1 < z2) then
            if ( z2 < zc(i,j) ) then
              dtc(i,j) = dt_cool(i,j)*(1.0-(z1+z2)/(zc(i,j)+zc(i,j)))
            elseif ( z1 < zc(i,j) .and. z2 >= zc(i,j) ) then
              dtc(i,j) = 0.5*(1.0-z1/zc(i,j))*dt_cool(i,j)*(zc(i,j)-z1)/(z2-z1)
            endif
          elseif ( z1 == z2 ) then
            if ( z1 < zc(i,j) ) then
              dtc(i,j) = dt_cool(i,j)*(1.0-z1/zc(i,j))
            endif
          endif
        endif
      endif        ! if ( slmsk(i,j) == 0 ) then
    enddo
  enddo 
!
! get the mean T departure from Tf in the range of z=z1 to z=z2

!$omp parallel do private(j,i)
  do j = 1, ny
    do i= 1, nx
      if ( slmsk(i,j) == 0.0 ) then
        dtm(i,j) = dtw(i,j) - dtc(i,j)
      endif
    enddo
  enddo

 end subroutine get_dtzm_2d

end module module_nst_water_prop
!! @}
