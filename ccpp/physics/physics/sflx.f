!\file sflx.f
!! This file is the entity of GFS Noah LSM Model(Version 2.7).

!>\defgroup Noah_LSM GFS Noah LSM Model
!!\brief This is the entity of GFS Noah LSM model of physics subroutines.
!! It is a soil/veg/snowpack land-surface model to update soil moisture, soil
!! ice, soil temperature, skin temperature, snowpack water content, snowdepth,
!! and all terms of the surface energy balance and surface water balance
!! (excluding input atmospheric forcings of downward radiation and 
!! precipitation ).
!!
!! The land-surface model component was substantially upgraded from the Oregon
!! State University (OSU) land surface model to EMC's new Noah Land Surface Model
!! (Noah LSM) during the major implementation in the NCEP Global Forecast System
!! (GFS) on May 31, 2005. Forecast System (GFS). The Noah LSM embodies about 10
!! years of upgrades (see \cite chen_et_al_1996; \cite koren_et_al_1999; 
!! \cite ek_et_al_2003) to its ancestor, the OSU LSM.  The Noah LSM upgrade includes:
!!  - An increase from two (10, 190 cm thick) to four soil layers (10, 30, 60, 100 cm thick)
!!  - Addition of frozen soil physics
!!  - Add glacial ice treatment
!!  - Two snowpack states (SWE, density)
!!  - New  formulations for infiltration and runoff account for sub-grid variability in precipitation and soil moisture
!!  - Revised physics of the snowpack and its influence on surface heat fluxes and albedo
!!  - Higher canopy resistance
!!  - Spatially  varying root depth
!!  - Surface fluxes weighted by snow cover fraction
!!  - Improved thermal conduction in soil/snow
!!  - Improved seasonality of green vegetation cover.
!!  - Improved evaporation treatment over bare soil and snowpack
!!
!!\param[in] nsoil     integer, number of soil layers (>=2 but <=nsold) 
!!\param[in] couple    integer, =0:uncoupled (land model only),    
!! =1:coupled with parent atmos model  
!!\param[in] icein     integer, sea-ice flag (=1: sea-ice, =0: land) 
!!\param[in] ffrozp    real, flag for snow-rain detection (1.=snow, 0.=rain)                                     
!!\param[in] dt        real, time step (<3600 sec)                
!!\param[in] zlvl      real, height abv atmos ground forcing vars (\f$m\f$)
!!\param[in] sldpth    real, thickness of each soil layer (\f$m\f$), nsoil 
!!\param[in] swdn      real, downward SW radiation flux (\f$W/m^2\f$)
!!\param[in] swnet     real, downward SW net (dn-up) flux (\f$W/m^2\f$)
!!\param[in] lwdn      real, downward LW radiation flux (\f$W/m^2\f$)
!!\param[in] sfcems    real, sfc LW emissivity (fractional) 
!!\param[in] sfcprs    real, pressure at height zlvl above ground(\f$Pa\f$) 
!!\param[in] sfctmp    real, air temp at height zlvl above ground (\f$K\f$)   
!!\param[in] sfcspd    real, wind speed at height zlvl above ground (\f$m s^{-1}\f$)  
!!\param[in] prcp      real, precipitation rate (\f$kgm^{-2}s^{-1}\f$)                    
!!\param[in] q2        real, mixing ratio at hght zlvl above ground (\f$kgkg^{-1}\f$) 
!!\param[in] q2sat     real, sat mixing ratio at zlvl above ground (\f$kgkg^{-1}\f$) 
!!\param[in] dqsdt2    real, slope of sat specific humidity curve at t=sfctmp (\f$kgkg^{-1}k^{-1}\f$) 
!!\param[in] th2       real, air potential temperature at zlvl above ground (\f$K\f$) 
!!\param[in] ivegsrc   integer, sfc veg type data source UMD or IGBP   
!!\param[in] vegtyp    integer, vegetation type (integer index) 
!!\param[in] soiltyp   integer, soil type (integer index)    
!!\param[in] slopetyp  integer, class of sfc slope (integer index) 
!!\param[in] shdmin    real, min areal coverage of green veg (fraction)  
!!\param[in] alb       real, background snow-free sfc albedo (fraction) 
!!\param[in] snoalb    real, max albedo over deep snow (fraction)  
!!\param[in,out] tbot     real, bottom soil temp (\f$K\f$) (local yearly-mean sfc air temp)  
!!\param[in,out] cmc      real, canopy moisture content (\f$m\f$)
!!\param[in,out] t1       real, ground/canopy/snowpack eff skin temp (\f$K\f$) 
!!\param[in,out] stc      real, soil temp (\f$K\f$)     
!!\param[in,out] smc      real, total soil moisture (vol fraction) 
!!\param[in,out] sh2o     real, unfrozen soil moisture (vol fraction), note: frozen part = smc-sh2o 
!!\param[in,out] sneqv    real, water-equivalent snow depth (\f$m\f$), note: snow density = snwqv/snowh 
!!\param[in,out] ch       real, sfc exchange coeff for heat & moisture (\f$ms^{-1}\f$), note: conductance since it's been mult by wind   
!!\param[in,out] cm       real, sfc exchange coeff for momentum (\f$ms^{-1}\f$), note: conductance since it's been mult by wind
!!\param[in,out] z0       real, roughness length (\f$m\f$)
!!\param[out] nroot    integer, number of root layers          
!!\param[out] shdfac   real, aeral coverage of green veg (fraction) 
!!\param[out] snowh    real, snow depth (\f$m\f$)         
!!\param[out] albedo   real, sfc albedo incl snow effect (fraction) 
!!\param[out] eta      real, downward latent heat flux (\f$W/m^2\f$) 
!!\param[out] sheat    real, downward sensible heat flux (\f$W/m^2\f$)  
!!\param[out] ec       real, canopy water evaporation (\f$W/m^2\f$)  
!!\param[out] edir     real, direct soil evaporation (\f$W/m^2\f$)
!!\param[out] et       real, plant transpiration (\f$W/m^2\f$)
!!\param[out] ett      real, total plant transpiration (\f$W/m^2\f$) 
!!\param[out] esnow    real, sublimation from snowpack (\f$W/m^2\f$)
!!\param[out] drip     real, through-fall of precip and/or dew in excess of canopy water-holding capacity (\f$m\f$)             
!!\param[out] dew      real, dewfall (or frostfall for t<273.15) (\f$m\f$) 
!!\param[out] beta     real, ratio of actual/potential evap     
!!\param[out] etp      real, potential evaporation (\f$W/m^2\f$)  
!!\param[out] ssoil    real, upward soil heat flux (\f$W/m^2\f$)
!!\param[out] flx1     real, precip-snow sfc flux  (\f$W/m^2\f$) 
!!\param[out] flx2     real, freezing rain latent heat flux (\f$W/m^2\f$) 
!!\param[out] flx3     real, phase-change heat flux from snowmelt (\f$W/m^2\f$) 
!!\param[out] runoff1  real, surface runoff (\f$ms^{-1}\f$) not infiltrating sfc 
!!\param[out] runoff2  real, sub sfc runoff (\f$ms^{-1}\f$) (baseflow)   
!!\param[out] runoff3  real, excess of porosity for a given soil layer 
!!\param[out] snomlt   real, snow melt (\f$m\f$) (water equivalent)
!!\param[out] sncovr   real, fractional snow cover
!!\param[out] rc       real, canopy resistance (s/m) 
!!\param[out] pc       real, plant coeff (fraction) where pc*etp=transpi 
!!\param[out] rsmin    real, minimum canopy resistance (s/m) 
!!\param[out] xlai     real, leaf area index  (dimensionless) 
!!\param[out] rcs      real, incoming solar rc factor (dimensionless) 
!!\param[out] rct      real, air temperature rc factor (dimensionless) 
!!\param[out] rcq      real, atoms vapor press deficit rc factor   
!!\param[out] rcsoil   real, soil moisture rc factor (dimensionless) 
!!\param[out] soilw    real, available soil moisture in root zone  
!!\param[out] soilm    real, total soil column moisture (frozen+unfrozen) (\f$m\f$)
!!\param[out] smcwlt   real, wilting point (volumetric)       
!!\param[out] smcdry   real, dry soil moisture threshold (volumetric) 
!!\param[out] smcref   real, soil moisture threshold (volumetric) 
!!\param[out] smcmax   real, porosity (sat val of soil mois) 
!!\section general_sflx GFS Noah LSM General Algorithm
!! @{
!      subroutine sflx                                                   & !  ---  inputs:
!ccppdox: avoid to connect to sflx in mpbl
      subroutine gfssflx                                                & !  ---  inputs:
     &     ( nsoil, couple, icein, ffrozp, dt, zlvl, sldpth,            &
     &       swdn, swnet, lwdn, sfcems, sfcprs, sfctmp,                 &
     &       sfcspd, prcp, q2, q2sat, dqsdt2, th2, ivegsrc,             &
     &       vegtyp, soiltyp, slopetyp, shdmin, alb, snoalb,            & !  ---  input/outputs:
     &       tbot, cmc, t1, stc, smc, sh2o, sneqv, ch, cm,z0,           & !  ---  outputs:
     &       nroot, shdfac, snowh, albedo, eta, sheat, ec,              &
     &       edir, et, ett, esnow, drip, dew, beta, etp, ssoil,         &
     &       flx1, flx2, flx3, runoff1, runoff2, runoff3,               &
     &       snomlt, sncovr, rc, pc, rsmin, xlai, rcs, rct, rcq,        &
     &       rcsoil, soilw, soilm, smcwlt, smcdry, smcref, smcmax)

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine sflx - version 2.7:                                       !
!  sub-driver for "noah/osu lsm" family of physics subroutines for a    !
!  soil/veg/snowpack land-surface model to update soil moisture, soil   !
!  ice, soil temperature, skin temperature, snowpack water content,     !
!  snowdepth, and all terms of the surface energy balance and surface   !
!  water balance (excluding input atmospheric forcings of downward      !
!  radiation and precip)                                                !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!      call sflx                                                        !
!  ---  inputs:                                                         !
!          ( nsoil, couple, icein, ffrozp, dt, zlvl, sldpth,            !
!            swdn, swnet, lwdn, sfcems, sfcprs, sfctmp,                 !
!            sfcspd, prcp, q2, q2sat, dqsdt2, th2,ivegsrc,              !
!            vegtyp, soiltyp, slopetyp, shdmin, alb, snoalb,            !
!  ---  input/outputs:                                                  !
!            tbot, cmc, t1, stc, smc, sh2o, sneqv, ch, cm,              !
!  ---  outputs:                                                        !
!            nroot, shdfac, snowh, albedo, eta, sheat, ec,              !
!            edir, et, ett, esnow, drip, dew, beta, etp, ssoil,         !
!            flx1, flx2, flx3, runoff1, runoff2, runoff3,               !
!            snomlt, sncovr, rc, pc, rsmin, xlai, rcs, rct, rcq,        !
!            rcsoil, soilw, soilm, smcwlt, smcdry, smcref, smcmax )     !
!                                                                       !
!                                                                       !
!  subprograms called:  redprm, snow_new, csnow, snfrac, alcalc,        !
!            tdfcnd, snowz0, sfcdif, penman, canres, nopac, snopac.     !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!    jun  2003  -- k. mitchell et. al -- created version 2.7            !
!         200x  -- sarah lu    modified the code including:             !
!                       added passing argument, couple; replaced soldn  !
!                       and solnet by radflx; call sfcdif if couple=0;  !
!                       apply time filter to stc and tskin; and the     !
!                       way of namelist inport.                         !
!    feb  2004 -- m. ek noah v2.7.1 non-linear weighting of snow vs     !
!                       non-snow covered portions of gridbox            !
!    apr  2009  -- y.-t. hou   added lw surface emissivity effect,      !
!                       streamlined and reformatted the code, and       !
!                       consolidated constents/parameters by using      !
!                       module physcons, and added program documentation!               !
!    sep  2009 -- s. moorthi minor fixes
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers (>=2 but <=nsold)  1    !
!     couple   - integer, =0:uncoupled (land model only)           1    !
!                         =1:coupled with parent atmos model            !
!     icein    - integer, sea-ice flag (=1: sea-ice, =0: land)     1    !
!     ffrozp   - real,                                             1    !
!     dt       - real, time step (<3600 sec)                       1    !
!     zlvl     - real, height abv atmos ground forcing vars (m)    1    !
!     sldpth   - real, thickness of each soil layer (m)          nsoil  !
!     swdn     - real, downward sw radiation flux (w/m**2)         1    !
!     swnet    - real, downward sw net (dn-up) flux (w/m**2)       1    !
!     lwdn     - real, downward lw radiation flux (w/m**2)         1    !
!     sfcems   - real, sfc lw emissivity (fractional)              1    !
!     sfcprs   - real, pressure at height zlvl abv ground(pascals) 1    !
!     sfctmp   - real, air temp at height zlvl abv ground (k)      1    !
!     sfcspd   - real, wind speed at height zlvl abv ground (m/s)  1    !
!     prcp     - real, precip rate (kg m-2 s-1)                    1    !
!     q2       - real, mixing ratio at hght zlvl abv grnd (kg/kg)  1    !
!     q2sat    - real, sat mixing ratio at zlvl abv grnd (kg/kg)   1    !
!     dqsdt2   - real, slope of sat specific humidity curve at     1    !
!                      t=sfctmp (kg kg-1 k-1)                           !
!     th2      - real, air potential temp at zlvl abv grnd (k)     1    !
!     ivegsrc  - integer, sfc veg type data source umd or igbp          !
!     vegtyp   - integer, vegetation type (integer index)          1    !
!     soiltyp  - integer, soil type (integer index)                1    !
!     slopetyp - integer, class of sfc slope (integer index)       1    !
!     shdmin   - real, min areal coverage of green veg (fraction)  1    !
!     alb      - real, bkground snow-free sfc albedo (fraction)    1    !
!     snoalb   - real, max albedo over deep snow     (fraction)    1    !
!                                                                       !
!  input/outputs:                                                       !
!     tbot     - real, bottom soil temp (k)                        1    !
!                      (local yearly-mean sfc air temp)                 !
!     cmc      - real, canopy moisture content (m)                 1    !
!     t1       - real, ground/canopy/snowpack eff skin temp (k)    1    !
!     stc      - real, soil temp (k)                             nsoil  !
!     smc      - real, total soil moisture (vol fraction)        nsoil  !
!     sh2o     - real, unfrozen soil moisture (vol fraction)     nsoil  !
!                      note: frozen part = smc-sh2o                     !
!     sneqv    - real, water-equivalent snow depth (m)             1    !
!                      note: snow density = snwqv/snowh                 !
!     ch       - real, sfc exchange coeff for heat & moisture (m/s)1    !
!                      note: conductance since it's been mult by wind   !
!     cm       - real, sfc exchange coeff for momentum (m/s)       1    !
!                      note: conductance since it's been mult by wind   !
!                                                                       !
!  outputs:                                                             !
!     nroot    - integer, number of root layers                    1    !
!     shdfac   - real, aeral coverage of green veg (fraction)      1    !
!     snowh    - real, snow depth (m)                              1    !
!     albedo   - real, sfc albedo incl snow effect (fraction)      1    !
!     eta      - real, downward latent heat flux (w/m2)            1    !
!     sheat    - real, downward sensible heat flux (w/m2)          1    !
!     ec       - real, canopy water evaporation (w/m2)             1    !
!     edir     - real, direct soil evaporation (w/m2)              1    !
!     et       - real, plant transpiration     (w/m2)            nsoil  !
!     ett      - real, total plant transpiration (w/m2)            1    !
!     esnow    - real, sublimation from snowpack (w/m2)            1    !
!     drip     - real, through-fall of precip and/or dew in excess 1    !
!                      of canopy water-holding capacity (m)             !
!     dew      - real, dewfall (or frostfall for t<273.15) (m)     1    !
!     beta     - real, ratio of actual/potential evap              1    !
!     etp      - real, potential evaporation (w/m2)                1    !
!     ssoil    - real, upward soil heat flux (w/m2)                1    !
!     flx1     - real, precip-snow sfc flux  (w/m2)                1    !
!     flx2     - real, freezing rain latent heat flux (w/m2)       1    !
!     flx3     - real, phase-change heat flux from snowmelt (w/m2) 1    !
!     snomlt   - real, snow melt (m) (water equivalent)            1    !
!     sncovr   - real, fractional snow cover                       1    !
!     runoff1  - real, surface runoff (m/s) not infiltrating sfc   1    !
!     runoff2  - real, sub sfc runoff (m/s) (baseflow)             1    !
!     runoff3  - real, excess of porosity for a given soil layer   1    !
!     rc       - real, canopy resistance (s/m)                     1    !
!     pc       - real, plant coeff (fraction) where pc*etp=transpi 1    !
!     rsmin    - real, minimum canopy resistance (s/m)             1    !
!     xlai     - real, leaf area index  (dimensionless)            1    !
!     rcs      - real, incoming solar rc factor  (dimensionless)   1    !
!     rct      - real, air temp rc factor        (dimensionless)   1    !
!     rcq      - real, atoms vapor press deficit rc factor         1    !
!     rcsoil   - real, soil moisture rc factor   (dimensionless)   1    !
!     soilw    - real, available soil mois in root zone            1    !
!     soilm    - real, total soil column mois (frozen+unfrozen) (m)1    !
!     smcwlt   - real, wilting point (volumetric)                  1    !
!     smcdry   - real, dry soil mois threshold (volumetric)        1    !
!     smcref   - real, soil mois threshold     (volumetric)        1    !
!     smcmax   - real, porosity (sat val of soil mois)             1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
      use machine ,   only : kind_phys
!
      use physcons,   only : con_cp, con_rd, con_t0c, con_g, con_pi,    &
     &                       con_cliq, con_csol, con_hvap, con_hfus,    &
     &                       con_sbc
!
      implicit none

!  ---  constant parameters:
!      *** note: some of the constants are different in subprograms and need to
!          be consolidated with the standard def in module physcons at sometime
!          at the present time, those diverse values are kept temperately to 
!          provide the same result as the original codes.  -- y.t.h.  may09

      integer,               parameter :: nsold   = 4           ! max soil layers

!     real (kind=kind_phys), parameter :: gs      = con_g       ! con_g   =9.80665
      real (kind=kind_phys), parameter :: gs1     = 9.8         ! con_g in sfcdif
      real (kind=kind_phys), parameter :: gs2     = 9.81        ! con_g in snowpack, frh2o
      real (kind=kind_phys), parameter :: tfreez  = con_t0c     ! con_t0c =275.15
      real (kind=kind_phys), parameter :: lsubc   = 2.501e+6    ! con_hvap=2.5000e+6
      real (kind=kind_phys), parameter :: lsubf   = 3.335e5     ! con_hfus=3.3358e+5
      real (kind=kind_phys), parameter :: lsubs   = 2.83e+6     ! ? in sflx, snopac
      real (kind=kind_phys), parameter :: elcp    = 2.4888e+3   ! ? in penman
!     real (kind=kind_phys), parameter :: rd      = con_rd      ! con_rd  =287.05
      real (kind=kind_phys), parameter :: rd1     = 287.04      ! con_rd in sflx, penman, canres
      real (kind=kind_phys), parameter :: cp      = con_cp      ! con_cp  =1004.6
      real (kind=kind_phys), parameter :: cp1     = 1004.5      ! con_cp in sflx, canres
      real (kind=kind_phys), parameter :: cp2     = 1004.0      ! con_cp in htr
!     real (kind=kind_phys), parameter :: cph2o   = con_cliq    ! con_cliq=4.1855e+3
      real (kind=kind_phys), parameter :: cph2o1  = 4.218e+3    ! con_cliq in penman, snopac
      real (kind=kind_phys), parameter :: cph2o2  = 4.2e6       ! con_cliq in hrt *unit diff!
      real (kind=kind_phys), parameter :: cpice   = con_csol    ! con_csol=2.106e+3
      real (kind=kind_phys), parameter :: cpice1  = 2.106e6     ! con_csol in hrt *unit diff!
!     real (kind=kind_phys), parameter :: sigma   = con_sbc     ! con_sbc=5.6704e-8
      real (kind=kind_phys), parameter :: sigma1  = 5.67e-8     ! con_sbc in penman, nopac, snopac

!  ---  inputs:
      integer, intent(in) :: nsoil, couple, icein, vegtyp, soiltyp,     &
     &       slopetyp, ivegsrc

      real (kind=kind_phys), intent(in) :: ffrozp, dt, zlvl, lwdn,      &
     &       sldpth(nsoil), swdn, swnet, sfcems, sfcprs, sfctmp,        &
     &       sfcspd, prcp, q2, q2sat, dqsdt2, th2, shdmin, alb, snoalb

!  ---  input/outputs:
      real (kind=kind_phys), intent(inout) :: tbot, cmc, t1, sneqv,     &
     &       stc(nsoil), smc(nsoil), sh2o(nsoil), ch, cm

!  ---  outputs:
      integer, intent(out) :: nroot

      real (kind=kind_phys), intent(out) :: shdfac, snowh, albedo,      &
     &       eta, sheat, ec, edir, et(nsoil), ett, esnow, drip, dew,    &
     &       beta, etp, ssoil, flx1, flx2, flx3, snomlt, sncovr,        &
     &       runoff1, runoff2, runoff3, rc, pc, rsmin, xlai, rcs,       &
     &       rct, rcq, rcsoil, soilw, soilm, smcwlt, smcdry, smcref,    &
     &       smcmax

!  ---  locals:
!     real (kind=kind_phys) ::  df1h,
      real (kind=kind_phys) ::  bexp, cfactr, cmcmax, csoil, czil,      &
     &       df1, df1a, dksat, dwsat, dsoil, dtot, frcsno,              &
     &       frcsoi, epsca, fdown, f1, fxexp, frzx, hs, kdt, prcp1,     &
     &       psisat, quartz, rch, refkdt, rr, rgl, rsmax, sndens,       &
     &       sncond, sbeta, sn_new, slope, snup, salp, soilwm, soilww,  &
     &       t1v, t24, t2v, th2v, topt, tsnow, zbot, z0

      real (kind=kind_phys), dimension(nsold) :: rtdis, zsoil

      logical :: frzgra, snowng

      integer :: ice, k, kz

!
!===> ...  begin here
!
!  --- ...  initialization

      runoff1 = 0.0
      runoff2 = 0.0
      runoff3 = 0.0
      snomlt  = 0.0

!  --- ...  define local variable ice to achieve:
!             sea-ice case,          ice =  1
!             non-glacial land,      ice =  0
!             glacial-ice land,      ice = -1
!             if vegtype=15 (glacial-ice), re-set ice flag = -1 (glacial-ice)
!    note - for open-sea, sflx should *not* have been called. set green
!           vegetation fraction (shdfac) = 0.
!> - For open-sea, sea-ice and glacial-ice cases, sflx() should not have
!! been called (set green vegetation fraction (shdfac) =0.)  
      ice = icein

      if(ivegsrc == 2) then
       if (vegtyp == 13) then
        ice = -1
        shdfac = 0.0
       endif
      endif

      if(ivegsrc == 1) then
       if (vegtyp == 15) then
        ice = -1
        shdfac = 0.0
       endif
      endif

!> - Calculate soil layer depth below ground (sigin of \a zsoil is negative).
      if (ice == 1) then

        shdfac = 0.0

!  --- ...  set green vegetation fraction (shdfac) = 0.
!           set sea-ice layers of equal thickness and sum to 3 meters

        do kz = 1, nsoil
          zsoil(kz) = -3.0 * float(kz) / float(nsoil)
        enddo

      else

!  --- ...  calculate depth (negative) below ground from top skin sfc to 
!           bottom of each soil layer.
!    note - sign of zsoil is negative (denoting below ground)

        zsoil(1) = -sldpth(1)
        do kz = 2, nsoil
          zsoil(kz) = -sldpth(kz) + zsoil(kz-1)
        end do

      endif   ! end if_ice_block
         
!  --- ...  next is crucial call to set the land-surface parameters, 
!           including soil-type and veg-type dependent parameters.
!           set shdfac=0.0 for bare soil surfaces

!> - Call redprm() to set the land-surface paramters,
!! including soil-type and veg-type dependent parameters.
      call redprm
        if(ivegsrc == 1) then
!only igbp type has urban
!urban
         if(vegtyp == 13)then
              shdfac=0.05
              rsmin=400.0
              smcmax = 0.45
              smcref = 0.42
              smcwlt = 0.40
              smcdry = 0.40
         endif
        endif

!  ---  inputs:                                                            !
!          ( nsoil, vegtyp, soiltyp, slopetyp, sldpth, zsoil,              !
!  ---  outputs:                                                           !
!            cfactr, cmcmax, rsmin, rsmax, topt, refkdt, kdt,              !
!            sbeta, shdfac, rgl, hs, zbot, frzx, psisat, slope,            !
!            snup, salp, bexp, dksat, dwsat, smcmax, smcwlt,               !
!            smcref, smcdry, f1, quartz, fxexp, rtdis, nroot,              !
!            z0, czil, xlai, csoil )                                       !

!  --- ...  initialize precipitation logicals.

      snowng = .false.
      frzgra = .false.

!> - Over sea-ice or glacial-ice, if water-equivalent snow depth (\a sneqv) below threshold
!! lower bound (0.01 m for sea-ice, 0.10 m for glacial-ice), then
!! set at lower bound and store the source increment in subsurface
!! runoff/baseflow (runoff2).
!    note - runoff2 is then a negative value (as a flag) over sea-ice or
!           glacial-ice, in order to achieve water balance.

      if (ice == 1) then

        if (sneqv < 0.01) then
          sneqv = 0.01
          snowh = 0.10
!         snowh = sneqv / sndens
        endif

      elseif (ice == -1) then

        if (sneqv < 0.10) then
!         sndens = sneqv / snowh
!         runoff2 = -(0.10 - sneqv) / dt
          sneqv = 0.10
          snowh = 1.00
!         snowh = sneqv / sndens
        endif

      endif   ! end if_ice_block

!> - For sea-ice and glacial-ice cases, set smc and sh2o values = 1.0
!! as a flag for non-soil medium.

      if (ice /= 0) then
        do kz = 1, nsoil
          smc (kz) = 1.0
          sh2o(kz) = 1.0
        enddo
      endif

!> - If input snowpack (\a sneqv) is nonzero, then call csnow() to compute 
!! snow density (\a sndens) and snow thermal conductivity (\a sncond). 
! (note that csnow is a function subroutine)

      if (sneqv .eq. 0.0) then
        sndens = 0.0
        snowh = 0.0
        sncond = 1.0
      else
        sndens = sneqv / snowh
        sndens = max( 0.0, min( 1.0, sndens ))   ! added by moorthi

        call csnow
!  ---  inputs:                                                         !
!          ( sndens,                                                    !
!  ---  outputs:                                                        !
!            sncond )                                                   !

      endif

!> - Determine if it's precipitating and what kind of precipitation it is.
!! if it's precipitating and the air temperature is colder than \f$0^oC\f$, 
!! it's snowing! if it's precipitating and the air temperature is warmer than 
!! \f$0^oC\f$, but the ground temperature is colder than \f$0^oC\f$, freezing 
!! rain is presumed to be falling.

      if (prcp > 0.0) then
        if (ffrozp > 0.5) then
          snowng = .true.
        else
          if (t1 <= tfreez) frzgra = .true.
        endif
      endif

!> - If either precipitation flag (\a snowng, \a frzgra) is set as true:
! determine new snowfall (converting precipitation rate from 
! \f$kg m^{-2} s^{-1}\f$ to a liquid equiv snow depth in meters)
!  and add it to the existing snowpack.
!>  - Since all precip is added to snowpack, no precip infiltrates
!! into the soil so that \a prcp1 is set to zero.

      if (snowng .or. frzgra) then

        sn_new = prcp * dt * 0.001
        sneqv = sneqv + sn_new
        prcp1 = 0.0

!>  - Call snow_new() to update snow density based on new snowfall, 
!! using old and new snow. 
        call snow_new
!  ---  inputs:                                                         !
!          ( sfctmp, sn_new,                                            !
!  ---  input/outputs:                                                  !
!            snowh, sndens )                                            !

!>  - Call csnow() to update snow thermal conductivity.
        call csnow
!  ---  inputs:                                                         !
!          ( sndens,                                                    !
!  ---  outputs:                                                        !
!            sncond )                                                   !

      else

!> - If precipitation is liquid (rain), hence save in the precip variable
!! that later can wholely or partially infiltrate the soil (along
!! with any canopy "drip" added to this later).

        prcp1 = prcp

      endif   ! end if_snowng_block

!> - Determine snowcover fraction and albedo fraction over sea-ice, 
!! glacial-ice, and land. For nonzero snow depth over land case:

      if (ice /= 0) then

!  --- ...  snow cover, albedo over sea-ice, glacial-ice

        sncovr = 1.0
        albedo = 0.65

      else

!  --- ...  non-glacial land
!           if snow depth=0, set snowcover fraction=0, albedo=snow free albedo.

        if (sneqv == 0.0) then

          sncovr = 0.0
          albedo = alb

        else

!  --- ...  determine snow fraction cover.
!           determine surface albedo modification due to snowdepth state.
!>  - Call snfrac() to calculate snow fraction cover.
          call snfrac
!  ---  inputs:                                                         !
!          ( sneqv, snup, salp, snowh,                                  !
!  ---  outputs:                                                        !
!            sncovr )                                                   !

!>  - Call alcalc() to calculate surface albedo modification due to snowdepth
!! state.
          call alcalc
!  ---  inputs:                                                         !
!          ( alb, snoalb, shdfac, shdmin, sncovr, tsnow,                !
!  ---  outputs:                                                        !
!            albedo )                                                   !

        endif   ! end if_sneqv_block

      endif   ! end if_ice_block

!  --- ...  thermal conductivity for sea-ice case, glacial-ice case
!> - Calculate thermal diffusivity (\a df1):
!>  - For sea-ice case and glacial-ice case, this is constant(\f$df1=2.2\f$).

      if (ice /= 0) then

        df1 = 2.2

      else
!>  - For non-glacial land case, call tdfcnd() to calculate the thermal
!! diffusivity of top soil layer (\cite peters-lidard_et_al_1998).

!  --- ...  next calculate the subsurface heat flux, which first requires
!           calculation of the thermal diffusivity.  treatment of the
!           latter follows that on pages 148-149 from "heat transfer in 
!           cold climates", by v. j. lunardini (published in 1981 
!           by van nostrand reinhold co.) i.e. treatment of two contiguous 
!           "plane parallel" mediums (namely here the first soil layer 
!           and the snowpack layer, if any). this diffusivity treatment 
!           behaves well for both zero and nonzero snowpack, including the 
!           limit of very thin snowpack.  this treatment also eliminates
!           the need to impose an arbitrary upper bound on subsurface 
!           heat flux when the snowpack becomes extremely thin.

!  --- ...   first calculate thermal diffusivity of top soil layer, using
!            both the frozen and liquid soil moisture, following the 
!            soil thermal diffusivity function of peters-lidard et al.
!            (1998,jas, vol 55, 1209-1224), which requires the specifying
!            the quartz content of the given soil class (see routine redprm)

        call tdfcnd                                                     &
!  ---  inputs:
     &     ( smc(1), quartz, smcmax, sh2o(1),                           &
!  ---  outputs:
     &       df1                                                        &
     &     )
!>   - For IGBP/urban, \f$df1=3.24\f$.
        if(ivegsrc == 1) then
!only igbp type has urban
!urban
            if ( vegtyp == 13 ) df1=3.24
        endif

!>   - Add subsurface heat flux reduction effect from the 
!!  overlying green canopy, adapted from section 2.1.2 of 
!!  \cite peters-lidard_et_al_1997.

        df1 = df1 * exp( sbeta*shdfac )

      endif   ! end if_ice_block

!  --- ...  finally "plane parallel" snowpack effect following 
!           v.j. linardini reference cited above. note that dtot is
!           combined depth of snowdepth and thickness of first soil layer

      dsoil = -0.5 * zsoil(1)

      if (sneqv == 0.0) then

        ssoil = df1 * (t1 - stc(1)) / dsoil

      else

        dtot = snowh + dsoil
        frcsno = snowh / dtot
        frcsoi = dsoil / dtot

!  --- ...  1. harmonic mean (series flow)

!       df1  = (sncond*df1) / (frcsoi*sncond + frcsno*df1)
!       df1h = (sncond*df1) / (frcsoi*sncond + frcsno*df1)

!  --- ...  2. arithmetic mean (parallel flow)

!       df1  = frcsno*sncond + frcsoi*df1
        df1a = frcsno*sncond + frcsoi*df1

!  --- ...  3. geometric mean (intermediate between harmonic and arithmetic mean)

!       df1 = (sncond**frcsno) * (df1**frcsoi)
!       df1 = df1h*sncovr + df1a*(1.0-sncovr)
!       df1 = df1h*sncovr + df1 *(1.0-sncovr)
        df1 = df1a*sncovr + df1 *(1.0-sncovr)

!> - Calculate subsurface heat flux, \a ssoil, from final thermal
!! diffusivity of surface mediums,\a df1 above, and skin
!! temperature and top mid-layer soil temperature.

        ssoil = df1 * (t1 - stc(1)) / dtot

      endif   ! end if_sneqv_block

!> - For uncoupled mode, call snowz0() to calculate surface roughness 
!! (\a z0) over snowpack using snow condition from the previous timestep.

!     if (couple == 0) then            ! uncoupled mode
        if (sncovr > 0.0) then

          call snowz0
!  ---  inputs:                                                         !
!          ( sncovr,                                                    !
!  ---  input/outputs:                                                  !
!            z0 )                                                       !

        endif
!     endif

!> - Calculate virtual temps and virtual potential temps needed by
!!           subroutines sfcdif and penman.

      t2v = sfctmp * (1.0 + 0.61*q2)

!  --- ...  next call routine sfcdif to calculate the sfc exchange coef (ch)
!           for heat and moisture.
!    note - comment out call sfcdif, if sfcdif already called in calling 
!           program (such as in coupled atmospheric model).
!         - do not call sfcdif until after above call to redprm, in case
!           alternative values of roughness length (z0) and zilintinkevich
!           coef (czil) are set there via namelist i/o.
!         - routine sfcdif returns a ch that represents the wind spd times
!           the "original" nondimensional "ch" typical in literature.  hence
!           the ch returned from sfcdif has units of m/s.  the important 
!           companion coefficient of ch, carried here as "rch", is the ch
!           from sfcdif times air density and parameter "cp".  "rch" is
!           computed in "call penman". rch rather than ch is the coeff 
!           usually invoked later in eqns.
!         - sfcdif also returns the surface exchange coefficient for momentum,
!           cm, also known as the surface drage coefficient, but cm is not
!           used here.

!  --- ...  key required radiation term is the total downward radiation
!           (fdown) = net solar (swnet) + downward longwave (lwdn),
!           for use in penman ep calculation (penman) and other surface
!           energy budget calcuations.  also need downward solar (swdn)
!           for canopy resistance routine (canres).
!    note - fdown, swdn are derived differently in the uncoupled and
!           coupled modes.

!> - Calculate the total downward radiation (\a fdown) = net solar (\a swnet) +
!!  downward longwave (\a lwdn) as input of penman() and other surface
!! energy budget calculations.

      if (couple == 0) then                      !......uncoupled mode

!  --- ...  uncoupled mode:
!           compute surface exchange coefficients

        t1v  = t1  * (1.0 + 0.61 * q2)
        th2v = th2 * (1.0 + 0.61 * q2)

        call sfcdif
!  ---  inputs:                                                         !
!          ( zlvl, z0, t1v, th2v, sfcspd, czil,                         !
!  ---  input/outputs:                                                  !
!            cm, ch )                                                   !

!     swnet = net solar radiation into the ground (w/m2; dn-up) from input
!     fdown  = net solar + downward lw flux at sfc (w/m2)

        fdown = swnet + lwdn

      else                                       !......coupled mode

!  --- ...  coupled mode (couple .ne. 0):
!           surface exchange coefficients computed externally and passed in,
!           hence subroutine sfcdif not called.

!     swnet = net solar radiation into the ground (w/m2; dn-up) from input
!     fdown  = net solar + downward lw flux at sfc (w/m2)

        fdown = swnet + lwdn

      endif   ! end if_couple_block

!> - Call penman() to calculate potential evaporation (\a etp),
!! and other partial products and sums for later
!! calculations.

      call penman
!  ---  inputs:                                                         !
!          ( sfctmp, sfcprs, sfcems, ch, t2v, th2, prcp, fdown,         !
!            ssoil, q2, q2sat, dqsdt2, snowng, frzgra,                  !
!  ---  outputs:                                                        !
!            t24, etp, rch, epsca, rr, flx2 )                           !

!> - Call canres() to calculate the canopy resistance and convert it
!! into pc if nonzero greenness fraction.

      if (shdfac > 0.) then

!  --- ...  frozen ground extension: total soil water "smc" was replaced 
!           by unfrozen soil water "sh2o" in call to canres below

        call canres
!  ---  inputs:                                                         !
!          ( nsoil, nroot, swdn, ch, q2, q2sat, dqsdt2, sfctmp,         !
!            sfcprs, sfcems, sh2o, smcwlt, smcref, zsoil, rsmin,        !
!            rsmax, topt, rgl, hs, xlai,                                !
!  ---  outputs:                                                        !
!            rc, pc, rcs, rct, rcq, rcsoil )                            !

      endif

!> - Now decide major pathway branch to take depending on whether
!!           snowpack exists or not:

      esnow = 0.0

      if (sneqv .eq. 0.0) then
!>  - For no snowpack is present, call nopac() to calculate soil moisture
!! and heat flux values and update soil moisture contant and soil heat
!! content values. 
        call nopac
!  ---  inputs:                                                         !
!          ( nsoil, nroot, etp, prcp, smcmax, smcwlt, smcref,           !
!            smcdry, cmcmax, dt, shdfac, sbeta, sfctmp, sfcems,         !
!            t24, th2, fdown, epsca, bexp, pc, rch, rr, cfactr,         !
!            slope, kdt, frzx, psisat, zsoil, dksat, dwsat,             !
!            zbot, ice, rtdis, quartz, fxexp, csoil,                    !
!  ---  input/outputs:                                                  !
!            cmc, t1, stc, sh2o, tbot,                                  !
!  ---  outputs:                                                        !
!            eta, smc, ssoil, runoff1, runoff2, runoff3, edir,          !
!            ec, et, ett, beta, drip, dew, flx1, flx3 )                 !

      else

!>  - For a snowpack is present, call snopac().
        call snopac
!  ---  inputs:                                                         !
!          ( nsoil, nroot, etp, prcp, smcmax, smcwlt, smcref, smcdry,   !
!            cmcmax, dt, df1, sfcems, sfctmp, t24, th2, fdown, epsca,   !
!            bexp, pc, rch, rr, cfactr, slope, kdt, frzx, psisat,       !
!            zsoil, dwsat, dksat, zbot, shdfac, ice, rtdis, quartz,     !
!            fxexp, csoil, flx2, snowng,                                !
!  ---  input/outputs:                                                  !
!            prcp1, cmc, t1, stc, sncovr, sneqv, sndens, snowh,         !
!            sh2o, tbot, beta,                                          !
!  ---  outputs:                                                        !
!            smc, ssoil, runoff1, runoff2, runoff3, edir, ec, et,       !
!            ett, snomlt, drip, dew, flx1, flx3, esnow )                !

      endif
!> - Noah LSM post-processing: 
!>  - Calculate sensible heat (h) for return to parent model.

      sheat = -(ch*cp1*sfcprs) / (rd1*t2v) * (th2 - t1)

!>  - Convert units and/or sign of total evap (eta), potential evap (etp),
!!  subsurface heat flux (s), and runoffs for what parent model expects.
!   convert eta from kg m-2 s-1 to w m-2
!     eta = eta * lsubc
!     etp = etp * lsubc

      edir = edir * lsubc
      ec = ec * lsubc

      do k = 1, 4
        et(k) = et(k) * lsubc
      enddo

      ett = ett * lsubc
      esnow = esnow * lsubs
      etp = etp * ((1.0 - sncovr)*lsubc + sncovr*lsubs)

      if (etp > 0.) then
        eta = edir + ec + ett + esnow
      else
        eta = etp
      endif

      beta = eta / etp

!>  - Convert the sign of soil heat flux so that:
!!   -  ssoil>0: warm the surface  (night time)
!!   -  ssoil<0: cool the surface  (day time)

      ssoil = -1.0 * ssoil      

      if (ice == 0) then

!>  - For the case of land (but not glacial-ice):
!!  convert runoff3 (internal layer runoff from supersat) from \f$m\f$ 
!!  to \f$ms^-1\f$ and add to subsurface runoff/baseflow (runoff2).
!!  runoff2 is already a rate at this point.

        runoff3 = runoff3 / dt
        runoff2 = runoff2 + runoff3

      else

!>  - For the case of sea-ice (ice=1) or glacial-ice (ice=-1), add any
!! snowmelt directly to surface runoff (runoff1) since there is no
!! soil medium, and thus no call to subroutine smflx (for soil
!! moisture tendency).

        runoff1 = snomlt / dt
      endif

!>  - Calculate total column soil moisture in meters (soilm) and root-zone 
!! soil moisture availability (fraction) relative to porosity/saturation.

      soilm = -1.0 * smc(1) * zsoil(1)
      do k = 2, nsoil
        soilm = soilm + smc(k)*(zsoil(k-1) - zsoil(k))
      enddo

      soilwm = -1.0 * (smcmax - smcwlt) * zsoil(1)
      soilww = -1.0 * (smc(1) - smcwlt) * zsoil(1)
      do k = 2, nroot
        soilwm = soilwm + (smcmax - smcwlt) * (zsoil(k-1) - zsoil(k))
        soilww = soilww + (smc(k) - smcwlt) * (zsoil(k-1) - zsoil(k))
      enddo

      soilw = soilww / soilwm
!
      return


! =================
      contains
! =================

!*************************************!
!  section-1  1st level subprograms   !
!*************************************!

!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates albedo including snow effect (0 -> 1).
!!\ingroup Noah_LSM
      subroutine alcalc
!...................................
!  ---  inputs:
!    &     ( alb, snoalb, shdfac, shdmin, sncovr, tsnow,                &
!  ---  outputs:
!    &       albedo                                                     &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine alcalc calculates albedo including snow effect (0 -> 1)   !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from calling program:                                  size   !
!     alb      - real, snowfree albedo                             1    !
!     snoalb   - real, maximum (deep) snow albedo                  1    !
!     shdfac   - real, areal fractional coverage of green veg.     1    !
!     shdmin   - real, minimum areal coverage of green veg.        1    !
!     sncovr   - real, fractional snow cover                       1    !
!     tsnow    - real, snow surface temperature (k)                1    !
!                                                                       !
!  outputs to calling program:                                          !
!     albedo   - real, surface albedo including snow effect        1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
!     real (kind=kind_phys), intent(in) :: alb, snoalb, shdfac,         &
!    &       shdmin, sncovr, tsnow

!  ---  outputs:
!     real (kind=kind_phys), intent(out) :: albedo

!  ---  locals: (none)

!
!===> ...  begin here
!
!  --- ...  snoalb is argument representing maximum albedo over deep snow,
!           as passed into sflx, and adapted from the satellite-based
!           maximum snow albedo fields provided by d. robinson and g. kukla
!           (1985, jcam, vol 24, 402-411)

!         albedo = alb + (1.0-(shdfac-shdmin))*sncovr*(snoalb-alb)
          albedo = alb + sncovr*(snoalb - alb)

          if (albedo > snoalb) albedo = snoalb

!  --- ...  base formulation (dickinson et al., 1986, cogley et al., 1990)

!     if (tsnow <= 263.16) then
!       albedo = snoalb
!     else
!       if (tsnow < 273.16) then
!         tm = 0.1 * (tsnow - 263.16)
!         albedo = 0.5 * ((0.9 - 0.2*(tm**3)) + (0.8 - 0.16*(tm**3)))
!       else
!         albedo = 0.67
!       endif
!     endif

!  --- ...  isba formulation (verseghy, 1991; baker et al., 1990)

!     if (tsnow < 273.16) then
!       albedo = snoalb - 0.008*dt/86400
!     else
!       albedo = (snoalb - 0.5) * exp( -0.24*dt/86400 ) + 0.5
!     endif

!
      return
!...................................
      end subroutine alcalc
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates canopy resistance which depends on incoming
!! solar radiation, air temperature, atmospheric water vapor pressure
!! deficit at the lowest model level, and soil moisture (preferably unfrozen
!! soil moisture rather than total).
      subroutine canres
!  ---  inputs:
!    &     ( nsoil, nroot, swdn, ch, q2, q2sat, dqsdt2, sfctmp,         &
!    &       sfcprs, sfcems, sh2o, smcwlt, smcref, zsoil, rsmin,        &
!    &       rsmax, topt, rgl, hs, xlai,                                &
!  ---  outputs:
!    &       rc, pc, rcs, rct, rcq, rcsoil                              &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine canres calculates canopy resistance which depends on      !
!  incoming solar radiation, air temperature, atmospheric water vapor   !
!  pressure deficit at the lowest model level, and soil moisture        !
!  (preferably unfrozen soil moisture rather than total)                !
!                                                                       !
!  source:   jarvis (1976), noilhan and planton (1989, mwr), jacquemin  !
!            and noilhan (1990, blm)                                    !
!  see also: chen et al (1996, jgr, vol 101(d3), 7251-7268), eqns       !
!            12-14 and table 2 of sec. 3.1.2                            !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from calling program:                                  size   !
!     nsoil    - integer, no. of soil layers                       1    !
!     nroot    - integer, no. of soil layers in root zone (<nsoil) 1    !
!     swdn     - real, incoming solar radiation                    1    !
!     ch       - real, sfc exchange coeff for heat and moisture    1    !
!     q2       - real, air humidity at 1st level above ground      1    !
!     q2sat    - real, sat. air humidity at 1st level abv ground   1    !
!     dqsdt2   - real, slope of sat. humidity function wrt temp    1    !
!     sfctmp   - real, sfc temperature at 1st level above ground   1    !
!     sfcprs   - real, sfc pressure                                1    !
!     sfcems   - real, sfc emissivity for lw radiation             1    !
!     sh2o     - real, volumetric soil moisture                  nsoil  !
!     smcwlt   - real, wilting point                               1    !
!     smcref   - real, reference soil moisture                     1    !
!     zsoil    - real, soil depth (negative sign, as below grd)  nsoil  !
!     rsmin    - real, mimimum stomatal resistance                 1    !
!     rsmax    - real, maximum stomatal resistance                 1    !
!     topt     - real, optimum transpiration air temperature       1    !
!     rgl      - real, canopy resistance func (in solar rad term)  1    !
!     hs       - real, canopy resistance func (vapor deficit term) 1    !
!     xlai     - real, leaf area index                             1    !
!                                                                       !
!  outputs to calling program:                                          !
!     rc       - real, canopy resistance                           1    !
!     pc       - real, plant coefficient                           1    !
!     rcs      - real, incoming solar rc factor                    1    !
!     rct      - real, air temp rc factor                          1    !
!     rcq      - real, atoms vapor press deficit rc factor         1    !
!     rcsoil   - real, soil moisture rc factor                     1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
!     integer, intent(in) :: nsoil, nroot

!     real (kind=kind_phys), intent(in) :: swdn, ch, q2, q2sat,         &
!    &       dqsdt2, sfctmp, sfcprs, sfcems, smcwlt, smcref, rsmin,     &
!    &       rsmax, topt, rgl, hs, xlai, sh2o(nsoil), zsoil(nsoil)

!  ---  outputs:
!     real (kind=kind_phys), intent(out) :: rc, pc, rcs, rct, rcq,      &
!    &       rcsoil

!  ---  locals:
      real (kind=kind_phys) :: delta, ff, gx, rr, part(nsold)

      integer :: k

!
!===> ...  begin here
!
!  --- ...  initialize canopy resistance multiplier terms.

      rcs = 0.0
      rct = 0.0
      rcq = 0.0
      rcsoil = 0.0
      rc = 0.0

!  --- ...  contribution due to incoming solar radiation

      ff = 0.55 * 2.0 * swdn / (rgl*xlai)
      rcs = (ff + rsmin/rsmax) / (1.0 + ff)
      rcs = max( rcs, 0.0001 )

!  --- ...  contribution due to air temperature at first model level above ground
!           rct expression from noilhan and planton (1989, mwr).

      rct = 1.0 - 0.0016 * (topt - sfctmp)**2.0
      rct = max( rct, 0.0001 )

!  --- ...  contribution due to vapor pressure deficit at first model level.
!           rcq expression from ssib

      rcq = 1.0 / (1.0 + hs*(q2sat-q2))
      rcq = max( rcq, 0.01 )

!  --- ...  contribution due to soil moisture availability.
!           determine contribution from each soil layer, then add them up.

      gx = (sh2o(1) - smcwlt) / (smcref - smcwlt)
      gx = max( 0.0, min( 1.0, gx ) )

!  --- ...  use soil depth as weighting factor
      part(1) = (zsoil(1)/zsoil(nroot)) * gx

!  --- ...  use root distribution as weighting factor
!     part(1) = rtdis(1) * gx

      do k = 2, nroot

        gx = (sh2o(k) - smcwlt) / (smcref - smcwlt)
        gx = max( 0.0, min( 1.0, gx ) )

!  --- ...  use soil depth as weighting factor
        part(k) = ((zsoil(k) - zsoil(k-1)) / zsoil(nroot)) * gx

!  --- ...  use root distribution as weighting factor
!       part(k) = rtdis(k) * gx

      enddo

      do k = 1, nroot
        rcsoil = rcsoil + part(k)
      enddo
      rcsoil = max( rcsoil, 0.0001 )

!  --- ...  determine canopy resistance due to all factors.  convert canopy
!           resistance (rc) to plant coefficient (pc) to be used with
!           potential evap in determining actual evap.  pc is determined by:
!           pc * linerized penman potential evap = penman-monteith actual
!           evaporation (containing rc term).

      rc = rsmin / (xlai*rcs*rct*rcq*rcsoil)
      rr = (4.0*sfcems*sigma1*rd1/cp1) * (sfctmp**4.0)/(sfcprs*ch) + 1.0
      delta = (lsubc/cp1) * dqsdt2

      pc = (rr + delta) / (rr*(1.0 + rc*ch) + delta)
!
      return
!...................................
      end subroutine canres
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates snow termal conductivity
      subroutine csnow
!...................................
!  ---  inputs:
!    &     ( sndens,                                                    &
!  ---  outputs:
!    &       sncond                                                     &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine csnow calculates snow termal conductivity                 !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from the calling program:                              size   !
!     sndens   - real, snow density                                1    !
!                                                                       !
!  outputs to the calling program:                                      !
!     sncond   - real, snow termal conductivity                    1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  constant parameters:
      real (kind=kind_phys), parameter :: unit = 0.11631

!  ---  inputs:
!     real (kind=kind_phys), intent(in) :: sndens

!  ---  outputs:
!     real (kind=kind_phys), intent(out) :: sncond

!  ---  locals:
      real (kind=kind_phys) :: c

!
!===> ...  begin here
!
!  --- ...  sncond in units of cal/(cm*hr*c), returned in w/(m*c)
!           basic version is dyachkova equation (1960), for range 0.1-0.4

      c = 0.328 * 10**(2.25*sndens)
      sncond = unit * c

!  --- ...  de vaux equation (1933), in range 0.1-0.6

!      sncond = 0.0293 * (1.0 + 100.0*sndens**2)

!  --- ...  e. andersen from flerchinger

!      sncond = 0.021 + 2.51 * sndens**2
!
      return
!...................................
      end subroutine csnow
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates soil moisture and heat flux values and 
!! update soil moisture content and soil heat content values for the
!! case when no snow pack is present.
      subroutine nopac
!...................................
!  ---  inputs:
!    &     ( nsoil, nroot, etp, prcp, smcmax, smcwlt, smcref,           &
!    &       smcdry, cmcmax, dt, shdfac, sbeta, sfctmp, sfcems,         &
!    &       t24, th2, fdown, epsca, bexp, pc, rch, rr, cfactr,         &
!    &       slope, kdt, frzx, psisat, zsoil, dksat, dwsat,             &
!    &       zbot, ice, rtdis, quartz, fxexp, csoil,                    &
!  ---  input/outputs:
!    &       cmc, t1, stc, sh2o, tbot,                                  &
!  ---  outputs:
!    &       eta, smc, ssoil, runoff1, runoff2, runoff3, edir,          &
!    &       ec, et, ett, beta, drip, dew, flx1, flx3                   &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine nopac calculates soil moisture and heat flux values and   !
!  update soil moisture content and soil heat content values for the    !
!  case when no snow pack is present.                                   !
!                                                                       !
!                                                                       !
!  subprograms called:  evapo, smflx, tdfcnd, shflx                     !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from calling program:                                  size   !
!     nsoil    - integer, number of soil layers                    1    !
!     nroot    - integer, number of root layers                    1    !
!     etp      - real, potential evaporation                       1    !
!     prcp     - real, precip rate                                 1    !
!     smcmax   - real, porosity (sat val of soil mois)             1    !
!     smcwlt   - real, wilting point                               1    !
!     smcref   - real, soil mois threshold                         1    !
!     smcdry   - real, dry soil mois threshold                     1    !
!     cmcmax   - real, maximum canopy water parameters             1    !
!     dt       - real, time step                                   1    !
!     shdfac   - real, aeral coverage of green veg                 1    !
!     sbeta    - real, param to cal veg effect on soil heat flux   1    !
!     sfctmp   - real, air temp at height zlvl abv ground          1    !
!     sfcems   - real, sfc lw emissivity                           1    !
!     t24      - real, sfctmp**4                                   1    !
!     th2      - real, air potential temp at zlvl abv grnd         1    !
!     fdown    - real, net solar + downward lw flux at sfc         1    !
!     epsca    - real,                                             1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     pc       - real, plant coeff                                 1    !
!     rch      - real, companion coefficient of ch                 1    !
!     rr       - real,                                             1    !
!     cfactr   - real, canopy water parameters                     1    !
!     slope    - real, linear reservoir coefficient                1    !
!     kdt      - real,                                             1    !
!     frzx     - real, frozen ground parameter                     1    !
!     psisat   - real, saturated soil potential                    1    !
!     zsoil    - real, soil layer depth below ground (negative)  nsoil  !
!     dksat    - real, saturated soil hydraulic conductivity       1    !
!     dwsat    - real, saturated soil diffusivity                  1    !
!     zbot     - real, specify depth of lower bd soil              1    !
!     ice      - integer, sea-ice flag (=1: sea-ice, =0: land)     1    !
!     rtdis    - real, root distribution                         nsoil  !
!     quartz   - real, soil quartz content                         1    !
!     fxexp    - real, bare soil evaporation exponent              1    !
!     csoil    - real, soil heat capacity                          1    !
!                                                                       !
!  input/outputs from and to the calling program:                       !
!     cmc      - real, canopy moisture content                     1    !
!     t1       - real, ground/canopy/snowpack eff skin temp        1    !
!     stc      - real, soil temp                                 nsoil  !
!     sh2o     - real, unfrozen soil moisture                    nsoil  !
!     tbot     - real, bottom soil temp                            1    !
!                                                                       !
!  outputs to the calling program:                                      !
!     eta      - real, downward latent heat flux                   1    !
!     smc      - real, total soil moisture                       nsoil  !
!     ssoil    - real, upward soil heat flux                       1    !
!     runoff1  - real, surface runoff not infiltrating sfc         1    !
!     runoff2  - real, sub surface runoff (baseflow)               1    !
!     runoff3  - real, excess of porosity                          1    !
!     edir     - real, direct soil evaporation                     1    !
!     ec       - real, canopy water evaporation                    1    !
!     et       - real, plant transpiration                       nsoil  !
!     ett      - real, total plant transpiration                   1    !
!     beta     - real, ratio of actual/potential evap              1    !
!     drip     - real, through-fall of precip and/or dew           1    !
!     dew      - real, dewfall (or frostfall)                      1    !
!     flx1     - real, precip-snow sfc flux                        1    !
!     flx3     - real, phase-change heat flux from snowmelt        1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
!     integer, intent(in) :: nsoil, nroot, ice

!     real (kind=kind_phys), intent(in) :: etp, prcp, smcmax,           &
!    &       smcwlt, smcref, smcdry, cmcmax, dt, shdfac, sbeta,         &
!    &       sfctmp, sfcems, t24, th2, fdown, epsca, bexp, pc,          &
!    &       rch, rr, cfactr, slope, kdt, frzx, psisat,                 &
!    &       zsoil(nsoil), dksat, dwsat, zbot, rtdis(nsoil),            &
!    &       quartz, fxexp, csoil

!  ---  input/outputs:
!     real (kind=kind_phys), intent(inout) :: cmc, t1, stc(nsoil),      &
!    &       sh2o(nsoil), tbot

!  ---  outputs:
!     real (kind=kind_phys), intent(out) :: eta, smc(nsoil), ssoil,     &
!    &       runoff1, runoff2, runoff3, edir, ec, et(nsoil), ett,       &
!    &       beta, drip, dew, flx1, flx3

!  ---  locals:
      real (kind=kind_phys) :: df1, eta1, etp1, prcp1, yy, yynum,       &
     &       zz1, ec1, edir1, et1(nsoil), ett1

      integer :: k

!
!===> ...  begin here
!
!  --- ...  convert etp from kg m-2 s-1 to ms-1 and initialize dew.

      prcp1= prcp * 0.001
      etp1 = etp  * 0.001
      dew  = 0.0
      edir = 0.0
      edir1= 0.0
      ec   = 0.0
      ec1  = 0.0

      do k = 1, nsoil
        et (k) = 0.0
        et1(k) = 0.0
      enddo

      ett  = 0.0
      ett1 = 0.0

      if (etp > 0.0) then

!  --- ...  convert prcp from 'kg m-2 s-1' to 'm s-1'.

        call evapo                                                      &
!  ---  inputs:
     &     ( nsoil, nroot, cmc, cmcmax, etp1, dt, zsoil,                &
     &       sh2o, smcmax, smcwlt, smcref, smcdry, pc,                  &
     &       shdfac, cfactr, rtdis, fxexp,                              &
!  ---  outputs:
     &       eta1, edir1, ec1, et1, ett1                                &
     &     )

        call smflx                                                      &
!  ---  inputs:
     &     ( nsoil, dt, kdt, smcmax, smcwlt, cmcmax, prcp1,             &
     &       zsoil, slope, frzx, bexp, dksat, dwsat, shdfac,            &
     &       edir1, ec1, et1,                                           &
!  ---  input/outputs:
     &       cmc, sh2o,                                                 &
!  ---  outputs:
     &       smc, runoff1, runoff2, runoff3, drip                       &
     &     )

      else

!  --- ...  if etp < 0, assume dew forms (transform etp1 into dew and
!           reinitialize etp1 to zero).

        eta1 = 0.0
        dew  = -etp1

!  --- ...  convert prcp from 'kg m-2 s-1' to 'm s-1' and add dew amount.

        prcp1 = prcp1 + dew

        call smflx                                                      &
!  ---  inputs:
     &     ( nsoil, dt, kdt, smcmax, smcwlt, cmcmax, prcp1,             &
     &       zsoil, slope, frzx, bexp, dksat, dwsat, shdfac,            &
     &       edir1, ec1, et1,                                           &
!  ---  input/outputs:
     &       cmc, sh2o,                                                 &
!  ---  outputs:
     &       smc, runoff1, runoff2, runoff3, drip                       &
     &     )

      endif   ! end if_etp_block

!  --- ...  convert modeled evapotranspiration fm  m s-1  to  kg m-2 s-1

      eta  = eta1 * 1000.0
      edir = edir1 * 1000.0
      ec   = ec1 * 1000.0

      do k = 1, nsoil
        et(k) = et1(k) * 1000.0
      enddo

      ett = ett1 * 1000.0

!  --- ...  based on etp and e values, determine beta

      if ( etp <= 0.0 ) then
        beta = 0.0
        if ( etp < 0.0 ) then
          beta = 1.0
        endif
      else
        beta = eta / etp
      endif

!  --- ...  get soil thermal diffuxivity/conductivity for top soil lyr,
!           calc. adjusted top lyr soil temp and adjusted soil flux, then
!           call shflx to compute/update soil heat flux and soil temps.

      call tdfcnd                                                       &
!  ---  inputs:
     &     ( smc(1), quartz, smcmax, sh2o(1),                           &
!  ---  outputs:
     &       df1                                                        &
     &     )
       if(ivegsrc == 1) then
!urban
         if ( vegtyp == 13 ) df1=3.24
       endif

!  --- ... vegetation greenness fraction reduction in subsurface heat
!          flux via reduction factor, which is convenient to apply here
!          to thermal diffusivity that is later used in hrt to compute
!          sub sfc heat flux (see additional comments on veg effect
!          sub-sfc heat flx in routine sflx)

      df1 = df1 * exp( sbeta*shdfac )

!  --- ...  compute intermediate terms passed to routine hrt (via routine
!           shflx below) for use in computing subsurface heat flux in hrt

      yynum = fdown - sfcems*sigma1*t24
      yy = sfctmp + (yynum/rch + th2 - sfctmp - beta*epsca)/rr
      zz1 = df1/(-0.5*zsoil(1)*rch*rr) + 1.0

      call shflx                                                        &
!  ---  inputs:
     &     ( nsoil, smc, smcmax, dt, yy, zz1, zsoil, zbot,              &
     &       psisat, bexp, df1, ice, quartz, csoil, vegtyp,             &
!  ---  input/outputs:
     &       stc, t1, tbot, sh2o,                                       &
!  ---  outputs:
     &       ssoil                                                      &
     &     )

!  --- ...  set flx1 and flx3 (snopack phase change heat fluxes) to zero since
!           they are not used here in snopac.  flx2 (freezing rain heat flux)
!           was similarly initialized in the penman routine.

      flx1 = 0.0
      flx3 = 0.0
!
      return
!...................................
      end subroutine nopac
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates potential evaporation for the current point.
!! various partial sums/products are also calculated and passed back
!! to the calling routine for later use
      subroutine penman
!...................................
!  ---  inputs:
!    &     ( sfctmp, sfcprs, sfcems, ch, t2v, th2, prcp, fdown,         &
!    &       ssoil, q2, q2sat, dqsdt2, snowng, frzgra,                  &
!  ---  outputs:
!    &       t24, etp, rch, epsca, rr, flx2                             &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine penman calculates potential evaporation for the current   !
!  point.  various partial sums/products are also calculated and passed !
!  back to the calling routine for later use.                           !
!                                                                       !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     sfctmp   - real, sfc temperature at 1st level above ground   1    !
!     sfcprs   - real, sfc pressure                                1    !
!     sfcems   - real, sfc emissivity for lw radiation             1    !
!     ch       - real, sfc exchange coeff for heat & moisture      1    !
!     t2v      - real, sfc virtual temperature                     1    !
!     th2      - real, air potential temp at zlvl abv grnd         1    !
!     prcp     - real, precip rate                                 1    !
!     fdown    - real, net solar + downward lw flux at sfc         1    !
!     ssoil    - real, upward soil heat flux                       1    !
!     q2       - real, mixing ratio at hght zlvl abv ground        1    !
!     q2sat    - real, sat mixing ratio at zlvl abv ground         1    !
!     dqsdt2   - real, slope of sat specific humidity curve        1    !
!     snowng   - logical, snow flag                                1    !
!     frzgra   - logical, freezing rain flag                       1    !
!                                                                       !
!  outputs:                                                             !
!     t24      - real, sfctmp**4                                   1    !
!     etp      - real, potential evaporation                       1    !
!     rch      - real, companion coefficient of ch                 1    !
!     epsca    - real,                                             1    !
!     rr       - real,                                             1    !
!     flx2     - real, freezing rain latent heat flux              1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
!     real (kind=kind_phys), intent(in) :: sfctmp, sfcprs, sfcems,      &
!    &       ch, t2v, th2, prcp, fdown, ssoil, q2, q2sat, dqsdt2

!     logical, intent(in) :: snowng, frzgra

!  ---  outputs:
!     real (kind=kind_phys), intent(out) :: t24, etp, rch, epsca,       &
!    &       rr, flx2

!  ---  locals:
      real (kind=kind_phys) :: a, delta, fnet, rad, rho

!
!===> ...  begin here
!
      flx2 = 0.0

!  --- ...  prepare partial quantities for penman equation.

      delta = elcp * dqsdt2
      t24 = sfctmp * sfctmp * sfctmp * sfctmp
      rr  = t24 * 6.48e-8 / (sfcprs*ch) + 1.0
      rho = sfcprs / (rd1*t2v)
      rch = rho * cp * ch

!  --- ...  adjust the partial sums / products with the latent heat
!           effects caused by falling precipitation.

      if (.not. snowng) then
        if (prcp > 0.0)  rr = rr + cph2o1*prcp/rch
      else
        rr = rr + cpice*prcp/rch
      endif

      fnet = fdown - sfcems*sigma1*t24 - ssoil

!  --- ...  include the latent heat effects of frzng rain converting to ice
!           on impact in the calculation of flx2 and fnet.

      if (frzgra) then
        flx2 = -lsubf * prcp
        fnet = fnet - flx2
      endif

!  --- ...  finish penman equation calculations.

      rad = fnet/rch + th2 - sfctmp
      a = elcp * (q2sat - q2)
      epsca = (a*rr + rad*delta) / (delta + rr)
      etp = epsca * rch / lsubc
!
      return
!...................................
      end subroutine penman
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine internally sets default values or optionally read-in 
!! via namelist i/o, all soil and vegetation parateters requied for the execusion
!! of the Noah LSM.
      subroutine redprm
!...................................
!  ---  inputs:
!    &     ( nsoil, vegtyp, soiltyp, slopetyp, sldpth, zsoil,              &
!  ---  outputs:
!    &       cfactr, cmcmax, rsmin, rsmax, topt, refkdt, kdt,              &
!    &       sbeta, shdfac, rgl, hs, zbot, frzx, psisat, slope,            &
!    &       snup, salp, bexp, dksat, dwsat, smcmax, smcwlt,               &
!    &       smcref, smcdry, f1, quartz, fxexp, rtdis, nroot,              &
!    &       z0, czil, xlai, csoil                                         &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine redprm internally sets(default valuess), or optionally    !
!  read-in via namelist i/o, all soil and vegetation parameters         !
!  required for the execusion of the noah lsm.                          !
!                                                                       !
!  optional non-default parameters can be read in, accommodating up to  !
!  30 soil, veg, or slope classes, if the default max number of soil,   !
!  veg, and/or slope types is reset.                                    !
!                                                                       !
!  future upgrades of routine redprm must expand to incorporate some    !
!  of the empirical parameters of the frozen soil and snowpack physics  !
!  (such as in routines frh2o, snowpack, and snow_new) not yet set in   !
!  this redprm routine, but rather set in lower level subroutines.      !
!                                                                       !
!  all soil, veg, slope, and universal parameters values are defined    !
!  externally (in subroutine "set_soilveg.f") and then accessed via     !
!  "use namelist_soilveg" (below) and then set here.                    !
!                                                                       !
!  soil types   zobler (1986)      cosby et al (1984) (quartz cont.(1)) !
!      1         coarse            loamy sand            (0.82)         !
!      2         medium            silty clay loam       (0.10)         !
!      3         fine              light clay            (0.25)         !
!      4         coarse-medium     sandy loam            (0.60)         !
!      5         coarse-fine       sandy clay            (0.52)         !
!      6         medium-fine       clay loam             (0.35)         !
!      7         coarse-med-fine   sandy clay loam       (0.60)         !
!      8         organic           loam                  (0.40)         !
!      9         glacial land ice  loamy sand            (na using 0.82)!
!     13: <old>- glacial land ice -<old>                                !
!     13:        glacial-ice (no longer use these parameters), now      !
!                treated as ice-only surface and sub-surface            !
!                (in subroutine hrtice)                                 !
!  upgraded to statsgo (19-type)
!     1: sand
!     2: loamy sand
!     3: sandy loam
!     4: silt loam
!     5: silt
!     6:loam
!     7:sandy clay loam
!     8:silty clay loam
!     9:clay loam
!     10:sandy clay
!     11: silty clay
!     12: clay
!     13: organic material
!     14: water
!     15: bedrock
!     16: other (land-ice)
!     17: playa
!     18: lava
!     19: white sand
!                                                                       !
!  ssib vegetation types (dorman and sellers, 1989; jam)                !
!      1:  broadleaf-evergreen trees  (tropical forest)                 !
!      2:  broadleaf-deciduous trees                                    !
!      3:  broadleaf and needleleaf trees (mixed forest)                !
!      4:  needleleaf-evergreen trees                                   !
!      5:  needleleaf-deciduous trees (larch)                           !
!      6:  broadleaf trees with groundcover (savanna)                   !
!      7:  groundcover only (perennial)                                 !
!      8:  broadleaf shrubs with perennial groundcover                  !
!      9:  broadleaf shrubs with bare soil                              !
!     10:  dwarf trees and shrubs with groundcover (tundra)             !
!     11:  bare soil                                                    !
!     12:  cultivations (the same parameters as for type 7)             !
!     13: <old>- glacial (the same parameters as for type 11) -<old>    !
!     13:  glacial-ice (no longer use these parameters), now treated as !
!          ice-only surface and sub-surface (in subroutine hrtice)      !
!  upgraded to IGBP (20-type)
!      1:Evergreen Needleleaf Forest
!      2:Evergreen Broadleaf Forest
!      3:Deciduous Needleleaf Forest
!      4:Deciduous Broadleaf Forest
!      5:Mixed Forests
!      6:Closed Shrublands
!      7:Open Shrublands
!      8:Woody Savannas
!      9:Savannas
!      10:Grasslands
!      11:Permanent wetlands
!      12:Croplands
!      13:Urban and Built-Up
!      14:Cropland/natural vegetation mosaic
!      15:Snow and Ice
!      16:Barren or Sparsely Vegetated
!      17:Water
!      18:Wooded Tundra
!      19:Mixed Tundra
!      20:Bare Ground Tundra
!                                                                       !
!  slopetyp is to estimate linear reservoir coefficient slope to the    !
!  baseflow runoff out of the bottom layer. lowest class (slopetyp=0)   !
!  means highest slope parameter = 1.                                   !
!                                                                       !
!  slope class       percent slope                                      !
!      1                0-8                                             !
!      2                8-30                                            !
!      3                > 30                                            !
!      4                0-30                                            !
!      5                0-8 & > 30                                      !
!      6                8-30 & > 30                                     !
!      7                0-8, 8-30, > 30                                 !
!      9                glacial ice                                     !
!    blank              ocean/sea                                       !
!                                                                       !
!  note: class 9 from zobler file should be replaced by 8 and 'blank' 9 !
!                                                                       !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from calling program:                                  size   !
!     nsoil    - integer, number of soil layers                    1    !
!     vegtyp   - integer, vegetation type (integer index)          1    !
!     soiltyp  - integer, soil type (integer index)                1    !
!     slopetyp - integer, class of sfc slope (integer index)       1    !
!     sldpth   - integer, thickness of each soil layer (m)       nsoil  !
!     zsoil    - integer, soil depth (negative sign) (m)         nsoil  !
!                                                                       !
!  outputs to the calling program:                                      !
!     cfactr   - real, canopy water parameters                     1    !
!     cmcmax   - real, maximum canopy water parameters             1    !
!     rsmin    - real, mimimum stomatal resistance                 1    !
!     rsmax    - real, maximum stomatal resistance                 1    !
!     topt     - real, optimum transpiration air temperature       1    !
!     refkdt   - real, =2.e-6 the sat. dk. val for soil type 2     1    !
!     kdt      - real,                                             1    !
!     sbeta    - real, param to cal veg effect on soil heat flux   1    !
!     shdfac   - real, vegetation greenness fraction               1    !
!     rgl      - real, canopy resistance func (in solar rad term)  1    !
!     hs       - real, canopy resistance func (vapor deficit term) 1    !
!     zbot     - real, specify depth of lower bd soil temp (m)     1    !
!     frzx     - real, frozen ground parameter, ice content        1    !
!                      threshold above which frozen soil is impermeable !
!     psisat   - real, saturated soil potential                    1    !
!     slope    - real, linear reservoir coefficient                1    !
!     snup     - real, threshold snow depth (water equi m)         1    !
!     salp     - real, snow cover shape parameter                  1    !
!                      from anderson's hydro-17 best fit salp = 2.6     !
!     bexp     - real, the 'b' parameter                           1    !
!     dksat    - real, saturated soil hydraulic conductivity       1    !
!     dwsat    - real, saturated soil diffusivity                  1    !
!     smcmax   - real, max soil moisture content (porosity)        1    !
!     smcwlt   - real, wilting pt soil moisture contents           1    !
!     smcref   - real, reference soil moisture (onset stress)      1    !
!     smcdry   - real, air dry soil moist content limits           1    !
!     f1       - real, used to comp soil diffusivity/conductivity  1    !
!     quartz   - real, soil quartz content                         1    !
!     fxexp    - real, bare soil evaporation exponent              1    !
!     rtdis    - real, root distribution                         nsoil  !
!     nroot    - integer, number of root layers                    1    !
!     z0       - real, roughness length (m)                        1    !
!     czil     - real, param to cal roughness length of heat       1    !
!     xlai     - real, leaf area index                             1    !
!     csoil    - real, soil heat capacity (j m-3 k-1)              1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
      use namelist_soilveg

!  ---  input:
!     integer, intent(in) :: nsoil, vegtyp, soiltyp, slopetyp

!     real (kind=kind_phys), intent(in) :: sldpth(nsoil), zsoil(nsoil)

!  ---  outputs:
!     real (kind=kind_phys), intent(out) :: cfactr, cmcmax, rsmin,      &
!    &       rsmax, topt, refkdt, kdt, sbeta, shdfac, rgl, hs, zbot,    &
!    &       frzx, psisat, slope, snup, salp, bexp, dksat, dwsat,       &
!    &       smcmax, smcwlt, smcref, smcdry, f1, quartz, fxexp, z0,     &
!    &       czil, xlai, csoil, rtdis(nsoil)

!     integer, intent(out) :: nroot

!  ---  locals:
      real (kind=kind_phys) :: frzfact, frzk, refdk

      integer :: i

!
!===> ...  begin here
!
      if (soiltyp > defined_soil) then
        write(*,*) 'warning: too many soil types,soiltyp=',soiltyp,     &
     &   'defined_soil=',defined_soil
        stop 333
      endif

      if (vegtyp > defined_veg) then
        write(*,*) 'warning: too many veg types'
        stop 333
      endif

      if (slopetyp > defined_slope) then
        write(*,*) 'warning: too many slope types'
        stop 333
      endif

!  --- ...  set-up universal parameters (not dependent on soiltyp, vegtyp
!           or slopetyp)

      zbot   = zbot_data
      salp   = salp_data
      cfactr = cfactr_data
      cmcmax = cmcmax_data
      sbeta  = sbeta_data
      rsmax  = rsmax_data
      topt   = topt_data
      refdk  = refdk_data
      frzk   = frzk_data
      fxexp  = fxexp_data
      refkdt = refkdt_data
      czil   = czil_data
      csoil  = csoil_data

!  --- ...  set-up soil parameters

      bexp  = bb   (soiltyp)
      dksat = satdk(soiltyp)
      dwsat = satdw(soiltyp)
      f1    = f11  (soiltyp)
      kdt   = refkdt * dksat / refdk

      psisat = satpsi(soiltyp)
      quartz = qtz   (soiltyp)
      smcdry = drysmc(soiltyp)
      smcmax = maxsmc(soiltyp)
      smcref = refsmc(soiltyp)
      smcwlt = wltsmc(soiltyp)

      frzfact = (smcmax / smcref) * (0.412 / 0.468)

!  --- ...  to adjust frzk parameter to actual soil type: frzk * frzfact

      frzx = frzk * frzfact

!  --- ...  set-up vegetation parameters

      nroot = nroot_data(vegtyp)
      snup  = snupx(vegtyp)
      rsmin = rsmtbl(vegtyp)

      rgl = rgltbl(vegtyp)
      hs  = hstbl(vegtyp)
! roughness lengthe is defined in sfcsub
!     z0  = z0_data(vegtyp)
      xlai= lai_data(vegtyp)

      if (vegtyp == bare) shdfac = 0.0

      if (nroot > nsoil) then
        write(*,*) 'warning: too many root layers'
        stop 333
      endif

!  --- ...  calculate root distribution.  present version assumes uniform
!           distribution based on soil layer depths.

      do i = 1, nroot
        rtdis(i) = -sldpth(i) / zsoil(nroot)
      enddo

!  --- ...  set-up slope parameter

      slope = slope_data(slopetyp)
!
      return
!...................................
      end subroutine redprm
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates surface layer exchange coefficients
!! via iterative process(see \cite chen_et_al_1997).
      subroutine sfcdif
!...................................
!  ---  inputs:
!    &     ( zlvl, z0, t1v, th2v, sfcspd, czil,                         &
!  ---  input/outputs:
!    &       cm, ch                                                     &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine sfcdif calculates surface layer exchange coefficients     !
!  via iterative process. see chen et al (1997, blm)                    !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from the calling program:                              size   !
!     zlvl     - real, height abv atmos ground forcing vars (m)    1    !
!     z0       - real, roughness length (m)                        1    !
!     t1v      - real, surface exchange coefficient                1    !
!     th2v     - real, surface exchange coefficient                1    !
!     sfcspd   - real, wind speed at height zlvl abv ground (m/s)  1    !
!     czil     - real, param to cal roughness length of heat       1    !
!                                                                       !
!  input/outputs from and to the calling program:                       !
!     cm       - real, sfc exchange coeff for momentum (m/s)       1    !
!     ch       - real, sfc exchange coeff for heat & moisture (m/s)1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  --- constant parameters:
      integer,               parameter :: itrmx  = 5
      real (kind=kind_phys), parameter :: wwst   = 1.2
      real (kind=kind_phys), parameter :: wwst2  = wwst*wwst
      real (kind=kind_phys), parameter :: vkrm   = 0.40
      real (kind=kind_phys), parameter :: excm   = 0.001
      real (kind=kind_phys), parameter :: beta   = 1.0/270.0
      real (kind=kind_phys), parameter :: btg    = beta*gs1
      real (kind=kind_phys), parameter :: elfc   = vkrm*btg
      real (kind=kind_phys), parameter :: wold   = 0.15
      real (kind=kind_phys), parameter :: wnew   = 1.0-wold
      real (kind=kind_phys), parameter :: pihf   = 3.14159265/2.0  ! con_pi/2.0

      real (kind=kind_phys), parameter :: epsu2  = 1.e-4
      real (kind=kind_phys), parameter :: epsust = 0.07
      real (kind=kind_phys), parameter :: ztmin  = -5.0
      real (kind=kind_phys), parameter :: ztmax  = 1.0
      real (kind=kind_phys), parameter :: hpbl   = 1000.0
      real (kind=kind_phys), parameter :: sqvisc = 258.2

      real (kind=kind_phys), parameter :: ric    = 0.183
      real (kind=kind_phys), parameter :: rric   = 1.0/ric
      real (kind=kind_phys), parameter :: fhneu  = 0.8
      real (kind=kind_phys), parameter :: rfc    = 0.191
      real (kind=kind_phys), parameter :: rfac   = ric/(fhneu*rfc*rfc)

!  ---  inputs:
!     real (kind=kind_phys),  intent(in) :: zlvl, z0, t1v, th2v,        &
!    &       sfcspd, czil

!  ---  input/outputs:
!     real (kind=kind_phys),  intent(inout) :: cm, ch

!  ---  locals:
      real (kind=kind_phys) :: zilfc, zu, zt, rdz, cxch, dthv, du2,     &
     &       btgh, wstar2, ustar, zslu, zslt, rlogu, rlogt, rlmo,       &
     &       zetalt, zetalu, zetau, zetat, xlu4, xlt4, xu4, xt4,        &
     &       xlu, xlt, xu, xt, psmz, simm, pshz, simh, ustark,          &
     &       rlmn, rlma

      integer :: ilech, itr

!  ---  define local in-line functions:

      real (kind=kind_phys) :: pslmu, pslms, pslhu, pslhs, zz
      real (kind=kind_phys) :: pspmu, pspms, psphu, psphs, xx, yy

!  ...  1) lech's surface functions

      pslmu( zz ) = -0.96 * log( 1.0-4.5*zz )
      pslms( zz ) = zz*rric - 2.076*(1.0 - 1.0/(zz + 1.0))
      pslhu( zz ) = -0.96 * log( 1.0-4.5*zz )
      pslhs( zz ) = zz*rfac - 2.076*(1.0 - 1.0/(zz + 1.0))

!  ...  2) paulson's surface functions

      pspmu( xx ) = -2.0 * log( (xx + 1.0)*0.5 )                        &
     &            - log( (xx*xx + 1.0)*0.5 ) + 2.0*atan(xx) - pihf
      pspms( yy ) = 5.0 * yy
      psphu( xx ) = -2.0 * log( (xx*xx + 1.0)*0.5 )
      psphs( yy ) = 5.0 * yy

!
!===> ...  begin here
!
!  --- ...  this routine sfcdif can handle both over open water (sea, ocean) and
!           over solid surface (land, sea-ice).

      ilech = 0

!   --- ...  ztfc: ratio of zoh/zom  less or equal than 1
!            czil: constant c in zilitinkevich, s. s.1995,:note about zt

      zilfc = -czil * vkrm * sqvisc

      zu = z0

      rdz = 1.0 / zlvl
      cxch = excm * rdz
      dthv = th2v - t1v
      du2 = max( sfcspd*sfcspd, epsu2 )

!  --- ...  beljars correction of ustar

      btgh = btg * hpbl

!  --- ...  if statements to avoid tangent linear problems near zero
      if (btgh*ch*dthv /= 0.0) then
        wstar2 = wwst2 * abs( btgh*ch*dthv )**(2.0/3.0)
      else
        wstar2 = 0.0
      endif

      ustar = max( sqrt( cm*sqrt( du2+wstar2 ) ), epsust )

!  --- ...  zilitinkevitch approach for zt

      zt = exp( zilfc*sqrt( ustar*z0 ) ) * z0

      zslu = zlvl + zu
      zslt = zlvl + zt

!     print*,'zslt=',zslt
!     print*,'zlvl=',zvll
!     print*,'zt=',zt

      rlogu = log( zslu/zu )
      rlogt = log( zslt/zt )

      rlmo = elfc*ch*dthv / ustar**3

!     print*,'rlmo=',rlmo
!     print*,'elfc=',elfc
!     print*,'ch=',ch
!     print*,'dthv=',dthv
!     print*,'ustar=',ustar

      do itr = 1, itrmx

!  --- ...  1./ monin-obukkhov length-scale

        zetalt = max( zslt*rlmo, ztmin )
        rlmo   = zetalt / zslt
        zetalu = zslu * rlmo
        zetau  = zu * rlmo
        zetat  = zt * rlmo

        if (ilech == 0) then

          if (rlmo < 0.0) then
            xlu4 = 1.0 - 16.0 * zetalu
            xlt4 = 1.0 - 16.0 * zetalt
            xu4  = 1.0 - 16.0 * zetau
            xt4  = 1.0 - 16.0* zetat

            xlu = sqrt( sqrt( xlu4 ) )
            xlt = sqrt( sqrt( xlt4 ) )
            xu  = sqrt( sqrt( xu4  ) )
            xt  = sqrt( sqrt( xt4  ) )

            psmz = pspmu(xu)

!           print*,'-----------1------------'
!           print*,'psmz=',psmz
!           print*,'pspmu(zetau)=',pspmu( zetau )
!           print*,'xu=',xu
!           print*,'------------------------'

            simm = pspmu( xlu ) - psmz + rlogu
            pshz = psphu( xt  )
            simh = psphu( xlt ) - pshz + rlogt
          else
            zetalu = min( zetalu, ztmax )
            zetalt = min( zetalt, ztmax )
            psmz = pspms( zetau )

!           print*,'-----------2------------'
!           print*,'psmz=',psmz
!           print*,'pspms(zetau)=',pspms( zetau )
!           print*,'zetau=',zetau
!           print*,'------------------------'

            simm = pspms( zetalu ) - psmz + rlogu
            pshz = psphs( zetat  )
            simh = psphs( zetalt ) - pshz + rlogt
          endif   ! end if_rlmo_block

        else

!  --- ...  lech's functions

          if (rlmo < 0.0) then
            psmz = pslmu( zetau )

!           print*,'-----------3------------'
!           print*,'psmz=',psmz
!           print*,'pslmu(zetau)=',pslmu( zetau )
!           print*,'zetau=',zetau
!           print*,'------------------------'

            simm = pslmu( zetalu ) - psmz + rlogu
            pshz = pslhu( zetat  )
            simh = pslhu( zetalt ) - pshz + rlogt
          else
            zetalu = min( zetalu, ztmax )
            zetalt = min( zetalt, ztmax )

            psmz = pslms( zetau  )

!           print*,'-----------4------------'
!           print*,'psmz=',psmz
!           print*,'pslms(zetau)=',pslms( zetau )
!           print*,'zetau=',zetau
!           print*,'------------------------'

            simm = pslms( zetalu ) - psmz + rlogu
            pshz = pslhs( zetat  )
            simh = pslhs( zetalt ) - pshz + rlogt
          endif   ! end if_rlmo_block

        endif   ! end if_ilech_block

!  --- ...  beljaars correction for ustar

        ustar = max( sqrt( cm*sqrt( du2+wstar2 ) ), epsust )

!  --- ...  zilitinkevitch fix for zt

        zt = exp( zilfc*sqrt( ustar*z0 ) ) * z0

        zslt = zlvl + zt
        rlogt = log( zslt/zt )

        ustark = ustar * vkrm
        cm = max( ustark/simm, cxch )
        ch = max( ustark/simh, cxch )

!  --- ...  if statements to avoid tangent linear problems near zero

        if (btgh*ch*dthv /= 0.0) then
          wstar2 = wwst2 * abs(btgh*ch*dthv) ** (2.0/3.0)
        else
          wstar2 = 0.0
        endif

        rlmn = elfc*ch*dthv / ustar**3
        rlma = rlmo*wold + rlmn*wnew

        rlmo = rlma

      enddo   ! end do_itr_loop

!     print*,'----------------------------'
!     print*,'sfcdif output !  ! ! ! ! ! ! ! !  !   !    !'
!
!     print*,'zlvl=',zlvl
!     print*,'z0=',z0
!     print*,'t1v=',t1v
!     print*,'th2v=',th2v
!     print*,'sfcspd=',sfcspd
!     print*,'czil=',czil
!     print*,'cm=',cm
!     print*,'ch=',ch
!     print*,'----------------------------'
!
      return
!...................................
      end subroutine sfcdif
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates snow fraction (0->1).
      subroutine snfrac
!...................................
!  ---  inputs:
!    &     ( sneqv, snup, salp, snowh,                                  &
!  ---  outputs:
!    &       sncovr                                                     &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine snfrac calculatexsnow fraction (0 -> 1)                   !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from the calling program:                              size   !
!     sneqv    - real, snow water equivalent (m)                   1    !
!     snup     - real, threshold sneqv depth above which sncovr=1  1    !
!     salp     - real, tuning parameter                            1    !
!     snowh    - real, snow depth (m)                              1    !
!                                                                       !
!  outputs to the calling program:                                      !
!     sncovr   - real, fractional snow cover                       1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
!     real (kind=kind_phys),  intent(in) :: sneqv, snup, salp, snowh

!  ---  outputs:
!     real (kind=kind_phys),  intent(out) :: sncovr

!  ---  locals:
      real (kind=kind_phys) :: rsnow, z0n

!
!===> ...  begin here
!
!  --- ...  snup is veg-class dependent snowdepth threshhold (set in routine
!           redprm) above which snocvr=1.

          if (sneqv < snup) then
            rsnow = sneqv / snup
            sncovr = 1.0 - (exp(-salp*rsnow) - rsnow*exp(-salp))
          else
            sncovr = 1.0
          endif

          z0n = 0.035

!  --- ...  formulation of dickinson et al. 1986

!       sncovr = snowh / (snowh + 5.0*z0n)

!  --- ...  formulation of marshall et al. 1994

!       sncovr = sneqv / (sneqv + 2.0*z0n)

!
      return
!...................................
      end subroutine snfrac
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates soil moisture and heat flux values and 
!! update soil moisture content and soil heat content values for the
!! case when a snow pack is present.
      subroutine snopac
!...................................
!  ---  inputs:
!    &     ( nsoil, nroot, etp, prcp, smcmax, smcwlt, smcref, smcdry,   &
!    &       cmcmax, dt, df1, sfcems, sfctmp, t24, th2, fdown, epsca,   &
!    &       bexp, pc, rch, rr, cfactr, slope, kdt, frzx, psisat,       &
!    &       zsoil, dwsat, dksat, zbot, shdfac, ice, rtdis, quartz,     &
!    &       fxexp, csoil, flx2, snowng,                                &
!  ---  input/outputs:
!    &       prcp1, cmc, t1, stc, sncovr, sneqv, sndens, snowh,         &
!    &       sh2o, tbot, beta,                                          &
!  ---  outputs:
!    &       smc, ssoil, runoff1, runoff2, runoff3, edir, ec, et,       &
!    &       ett, snomlt, drip, dew, flx1, flx3, esnow                  &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine snopac calculates soil moisture and heat flux values and  !
!  update soil moisture content and soil heat content values for the    !
!  case when a snow pack is present.                                    !
!                                                                       !
!                                                                       !
!  subprograms called:  evapo, smflx, shflx, snowpack
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from the calling program:                              size   !
!     nsoil    - integer, number of soil layers                    1    !
!     nroot    - integer, number of root layers                    1    !
!     etp      - real, potential evaporation                       1    !
!     prcp     - real, precip rate                                 1    !
!     smcmax   - real, porosity                                    1    !
!     smcwlt   - real, wilting point                               1    !
!     smcref   - real, soil mois threshold                         1    !
!     smcdry   - real, dry soil mois threshold                     1    !
!     cmcmax   - real, maximum canopy water parameters             1    !
!     dt       - real, time step                                   1    !
!     df1      - real, thermal diffusivity                         m    !
!     sfcems   - real, lw surface emissivity                       1    !
!     sfctmp   - real, sfc temperature                             1    !
!     t24      - real, sfctmp**4                                   1    !
!     th2      - real, sfc air potential temperature               1    !
!     fdown    - real, net solar + downward lw flux at sfc         1    !
!     epsca    - real,                                             1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     pc       - real, plant coeff                                 1    !
!     rch      - real, companion coefficient of ch                 1    !
!     rr       - real,                                             1    !
!     cfactr   - real, canopy water parameters                     1    !
!     slope    - real, linear reservoir coefficient                1    !
!     kdt      - real,                                             1    !
!     frzx     - real, frozen ground parameter                     1    !
!     psisat   - real, saturated soil potential                    1    !
!     zsoil    - real, soil layer depth below ground (negative)  nsoil  !
!     dwsat    - real, saturated soil diffusivity                  1    !
!     dksat    - real, saturated soil hydraulic conductivity       1    !
!     zbot     - real, specify depth of lower bd soil              1    !
!     shdfac   - real, aeral coverage of green vegetation          1    !
!     ice      - integer, sea-ice flag (=1: sea-ice, =0: land)     1    !
!     rtdis    - real, root distribution                         nsoil  !
!     quartz   - real, soil quartz content                         1    !
!     fxexp    - real, bare soil evaporation exponent              1    !
!     csoil    - real, soil heat capacity                          1    !
!     flx2     - real, freezing rain latent heat flux              1    !
!     snowng   - logical, snow flag                                1    !
!                                                                       !
!  input/outputs from and to the calling program:                       !
!     prcp1    - real, effective precip                            1    !
!     cmc      - real, canopy moisture content                     1    !
!     t1       - real, ground/canopy/snowpack eff skin temp        1    !
!     stc      - real, soil temperature                          nsoil  !
!     sncovr   - real, snow cover                                  1    !
!     sneqv    - real, water-equivalent snow depth                 1    !
!     sndens   - real, snow density                                1    !
!     snowh    - real, snow depth                                  1    !
!     sh2o     - real, unfrozen soil moisture                    nsoil  !
!     tbot     - real, bottom soil temperature                     1    !
!     beta     - real, ratio of actual/potential evap              1    !
!                                                                       !
!  outputs to the calling program:                                      !
!     smc      - real, total soil moisture                       nsoil  !
!     ssoil    - real, upward soil heat flux                       1    !
!     runoff1  - real, surface runoff not infiltrating sfc         1    !
!     runoff2  - real, sub surface runoff                          1    !
!     runoff3  - real, excess of porosity for a given soil layer   1    !
!     edir     - real, direct soil evaporation                     1    !
!     ec       - real, canopy water evaporation                    1    !
!     et       - real, plant transpiration                       nsoil  !
!     ett      - real, total plant transpiration                   1    !
!     snomlt   - real, snow melt water equivalent                  1    !
!     drip     - real, through-fall of precip                      1    !
!     dew      - real, dewfall (or frostfall)                      1    !
!     flx1     - real, precip-snow sfc flux                        1    !
!     flx3     - real, phase-change heat flux from snowmelt        1    !
!     esnow    - real, sublimation from snowpack                   1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  constant parameters:
      real, parameter :: esdmin = 1.e-6

!  ---  inputs:
!     integer, intent(in) :: nsoil, nroot, ice

!     real (kind=kind_phys), intent(in) :: etp, prcp, smcmax, smcref,   &
!    &       smcwlt, smcdry, cmcmax, dt, df1, sfcems, sfctmp, t24,      &
!    &       th2, fdown, epsca, bexp, pc, rch, rr, cfactr, slope, kdt,  &
!    &       frzx, psisat, dwsat, dksat, zbot, shdfac, quartz,          &
!    &       csoil, fxexp, flx2, zsoil(nsoil), rtdis(nsoil)

!     logical, intent(in) :: snowng

!  ---  input/outputs:
!     real (kind=kind_phys), intent(inout) :: prcp1, t1, sncovr, sneqv, &
!    &       sndens, snowh, cmc, tbot, beta, sh2o(nsoil), stc(nsoil)

!  ---  outputs:
!     real (kind=kind_phys), intent(out) :: ssoil, runoff1, runoff2,    &
!    &       runoff3, edir, ec, et(nsoil), ett, snomlt, drip, dew,      &
!    &       flx1, flx3, esnow, smc(nsoil)

!  ---  locals:
      real (kind=kind_phys):: denom, dsoil, dtot, etp1, ssoil1,         &
     &       snoexp, ex, t11, t12, t12a, t12b, yy, zz1, seh, t14,       &
     &       ec1, edir1, ett1, etns, etns1, esnow1, esnow2, etanrg,     &
     &       et1(nsoil)

      integer k

!     data snoexp /1.0/    !!! <----- for noah v2.7
      data snoexp /2.0/    !!! <----- for noah v2.7.1

!  --- ...  convert potential evap (etp) from kg m-2 s-1 to m s-1 and then to an
!           amount (m) given timestep (dt) and call it an effective snowpack
!           reduction amount, esnow2 (m) for a snowcover fraction = 1.0.  this is
!           the amount the snowpack would be reduced due to sublimation from the
!           snow sfc during the timestep.  sublimation will proceed at the
!           potential rate unless the snow depth is less than the expected
!           snowpack reduction.  for snowcover fraction = 1.0, 0=edir=et=ec, and
!           hence total evap = esnow = sublimation (potential evap rate)

!  --- ...  if sea-ice (ice=1) or glacial-ice (ice=-1), snowcover fraction = 1.0,
!           and sublimation is at the potential rate.
!           for non-glacial land (ice=0), if snowcover fraction < 1.0, total
!           evaporation < potential due to non-potential contribution from
!           non-snow covered fraction.

      prcp1 = prcp1 * 0.001

      edir  = 0.0
      edir1 = 0.0

      ec  = 0.0
      ec1 = 0.0

      do k = 1, nsoil
        et (k) = 0.0
        et1(k) = 0.0
      enddo

      ett   = 0.0
      ett1  = 0.0
      etns  = 0.0
      etns1 = 0.0
      esnow = 0.0
      esnow1= 0.0
      esnow2= 0.0

      dew = 0.0
      etp1 = etp * 0.001

      if (etp < 0.0) then

!  --- ...  if etp<0 (downward) then dewfall (=frostfall in this case).

        dew = -etp1
        esnow2 = etp1 * dt
        etanrg = etp * ((1.0-sncovr)*lsubc + sncovr*lsubs)

      else

!  --- ...  etp >= 0, upward moisture flux

        if (ice /= 0) then           ! for sea-ice and glacial-ice case

          esnow = etp
          esnow1 = esnow * 0.001
          esnow2 = esnow1 * dt
          etanrg = esnow * lsubs

        else                         ! for non-glacial land case

          if (sncovr < 1.0) then

            call evapo                                                  &
!  ---  inputs:
     &     ( nsoil, nroot, cmc, cmcmax, etp1, dt, zsoil,                &
     &       sh2o, smcmax, smcwlt, smcref, smcdry, pc,                  &
     &       shdfac, cfactr, rtdis, fxexp,                              &
!  ---  outputs:
     &       etns1, edir1, ec1, et1, ett1                               &
     &     )

            edir1 = edir1 * (1.0 - sncovr)
            ec1 = ec1 * (1.0 - sncovr)

            do k = 1, nsoil
              et1(k) = et1(k) * (1.0 - sncovr)
            enddo

            ett1  = ett1  * (1.0 - sncovr)
            etns1 = etns1 * (1.0 - sncovr)

            edir = edir1 * 1000.0
            ec = ec1 * 1000.0

            do k = 1, nsoil
              et(k) = et1(k) * 1000.0
            enddo

            ett = ett1 * 1000.0
            etns = etns1 * 1000.0

          endif   ! end if_sncovr_block

          esnow = etp * sncovr
!         esnow1 = etp * 0.001
          esnow1 = esnow * 0.001
          esnow2 = esnow1 * dt
          etanrg = esnow*lsubs + etns*lsubc

        endif   ! end if_ice_block

      endif   ! end if_etp_block

!  --- ...  if precip is falling, calculate heat flux from snow sfc to newly
!           accumulating precip.  note that this reflects the flux appropriate for
!           the not-yet-updated skin temperature (t1).  assumes temperature of the
!           snowfall striking the gound is =sfctmp (lowest model level air temp).

      flx1 = 0.0
      if ( snowng ) then
        flx1 = cpice * prcp * (t1 - sfctmp)
      else
        if (prcp > 0.0) flx1 = cph2o1 * prcp * (t1 - sfctmp)
      endif

!  --- ...  calculate an 'effective snow-grnd sfc temp' (t12) based on heat
!           fluxes between the snow pack and the soil and on net radiation.
!           include flx1 (precip-snow sfc) and flx2 (freezing rain latent
!           heat) fluxes.
!           flx2 reflects freezing rain latent heat flux using t1 calculated
!           in penman.

      dsoil = -0.5 * zsoil(1)
      dtot = snowh + dsoil
      denom = 1.0 + df1 / (dtot * rr * rch)

!     t12a = ( (fdown - flx1 - flx2 - sigma1*t24) / rch                 &
!    &     + th2 - sfctmp - beta*epsca ) / rr
      t12a = ( (fdown - flx1 - flx2 - sfcems*sigma1*t24) / rch          &
     &     + th2 - sfctmp - etanrg/rch ) / rr

      t12b = df1 * stc(1) / (dtot * rr * rch)
      t12 = (sfctmp + t12a + t12b) / denom

!  --- ...  if the 'effective snow-grnd sfc temp' is at or below freezing, no snow
!           melt will occur.  set the skin temp to this effective temp.  reduce
!           (by sublimination ) or increase (by frost) the depth of the snowpack,
!           depending on sign of etp.
!           update soil heat flux (ssoil) using new skin temperature (t1)
!           since no snowmelt, set accumulated snowmelt to zero, set 'effective'
!           precip from snowmelt to zero, set phase-change heat flux from snowmelt
!           to zero.

      if (t12 <= tfreez) then

        t1 = t12
!       ssoil = df1 * (t1 - stc(1)) / dtot
        ssoil = (t1 - stc (1)) * max(7.0, df1/dtot)
        sneqv = max(0.0, sneqv-esnow2)
        flx3 = 0.0
        ex = 0.0
        snomlt = 0.0

      else

!  --- ...  if the 'effective snow-grnd sfc temp' is above freezing, snow melt
!           will occur.  call the snow melt rate,ex and amt, snomlt.  revise the
!           effective snow depth.  revise the skin temp because it would have chgd
!           due to the latent heat released by the melting. calc the latent heat
!           released, flx3. set the effective precip, prcp1 to the snow melt rate,
!           ex for use in smflx.  adjustment to t1 to account for snow patches.
!           calculate qsat valid at freezing point.  note that esat (saturation
!           vapor pressure) value of 6.11e+2 used here is that valid at frzzing
!           point.  note that etp from call penman in sflx is ignored here in
!           favor of bulk etp over 'open water' at freezing temp.
!           update soil heat flux (s) using new skin temperature (t1)

!  --- ...  noah v2.7.1   mek feb2004
!           non-linear weighting of snow vs non-snow covered portions of gridbox
!           so with snoexp = 2.0 (>1), surface skin temperature is higher than
!           for the linear case (snoexp = 1).

        t1 = tfreez * sncovr**snoexp + t12 * (1.0 - sncovr**snoexp)

        beta = 1.0
        ssoil = df1 * (t1 - stc(1)) / dtot

!  --- ...  if potential evap (sublimation) greater than depth of snowpack.
!           beta<1
!           snowpack has sublimated away, set depth to zero.

        if (sneqv-esnow2 <= esdmin) then

          sneqv = 0.0
          ex = 0.0
          snomlt = 0.0
          flx3 = 0.0

        else

!  --- ...  potential evap (sublimation) less than depth of snowpack, retain
!           beta=1.

          sneqv = sneqv - esnow2
          seh = rch * (t1 - th2)

          t14 = t1 * t1
          t14 = t14 * t14

          flx3 = fdown - flx1 - flx2 - sfcems*sigma1*t14                &
     &         - ssoil - seh - etanrg
          if (flx3 <= 0.0) flx3 = 0.0

          ex = flx3 * 0.001 / lsubf

!  --- ...  snowmelt reduction depending on snow cover
!           if snow cover less than 5% no snowmelt reduction
!     note: does 'if' below fail to match the melt water with the melt
!           energy?

!         if (sncovr > 0.05) ex = ex * sncovr
          snomlt = ex * dt

!  --- ...  esdmin represents a snowpack depth threshold value below which we
!           choose not to retain any snowpack, and instead include it in snowmelt.

          if (sneqv-snomlt >= esdmin) then

            sneqv = sneqv - snomlt

          else

!  --- ...  snowmelt exceeds snow depth

            ex = sneqv / dt
            flx3 = ex * 1000.0 * lsubf
            snomlt = sneqv 
            sneqv = 0.0

          endif   ! end if_sneqv-snomlt_block

        endif   ! end if_sneqv-esnow2_block

!       prcp1 = prcp1 + ex

!  --- ...  if non-glacial land, add snowmelt rate (ex) to precip rate to be used
!           in subroutine smflx (soil moisture evolution) via infiltration.

!  --- ...  for sea-ice and glacial-ice, the snowmelt will be added to subsurface
!           runoff/baseflow later near the end of sflx (after return from call to
!           subroutine snopac)

        if (ice == 0) prcp1 = prcp1 + ex

      endif   ! end if_t12<=tfreez_block

!  --- ...  final beta now in hand, so compute evaporation.  evap equals etp
!           unless beta<1.

!      eta = beta * etp

!  --- ...  smflx returns updated soil moisture values for non-glacial land.
!           if sea-ice (ice=1) or glacial-ice (ice=-1), skip call to smflx, since
!           no soil medium for sea-ice or glacial-ice

      if (ice == 0) then

        call smflx                                                        &
!  ---  inputs:
     &     ( nsoil, dt, kdt, smcmax, smcwlt, cmcmax, prcp1,               &
     &       zsoil, slope, frzx, bexp, dksat, dwsat, shdfac,              &
     &       edir1, ec1, et1,                                             &
!  ---  input/outputs:
     &       cmc, sh2o,                                                   &
!  ---  outputs:
     &       smc, runoff1, runoff2, runoff3, drip                         &
     &     )

      endif

!  --- ...  before call shflx in this snowpack case, set zz1 and yy arguments to
!           special values that ensure that ground heat flux calculated in shflx
!           matches that already computer for below the snowpack, thus the sfc
!           heat flux to be computed in shflx will effectively be the flux at the
!           snow top surface.  t11 is a dummy arguement so we will not use the
!           skin temp value as revised by shflx.

      zz1 = 1.0
      yy = stc(1) - 0.5*ssoil*zsoil(1)*zz1 / df1
      t11 = t1

!  --- ...  shflx will calc/update the soil temps.  note:  the sub-sfc heat flux
!           (ssoil1) and the skin temp (t11) output from this shflx call are not
!           used  in any subsequent calculations. rather, they are dummy variables
!           here in the snopac case, since the skin temp and sub-sfc heat flux are
!           updated instead near the beginning of the call to snopac.

      call shflx                                                        &
!  ---  inputs:
     &     ( nsoil, smc, smcmax, dt, yy, zz1, zsoil, zbot,              &
     &       psisat, bexp, df1, ice, quartz, csoil, vegtyp,             &
!  ---  input/outputs:
     &       stc, t11, tbot, sh2o,                                      &
!  ---  outputs:
     &       ssoil1                                                     &
     &     )

!  --- ...  snow depth and density adjustment based on snow compaction.  yy is
!           assumed to be the soil temperture at the top of the soil column.

      if (ice == 0) then              ! for non-glacial land

        if (sneqv > 0.0) then

          call snowpack                                                 &
!  ---  inputs:
     &     ( sneqv, dt, t1, yy,                                         &
!  ---  input/outputs:
     &       snowh, sndens                                              &
     &     )

        else

          sneqv = 0.0
          snowh = 0.0
          sndens = 0.0
!         sncond = 1.0
          sncovr = 0.0

        endif   ! end if_sneqv_block

!  --- ...  over sea-ice or glacial-ice, if s.w.e. (sneqv) below threshold lower
!           bound (0.01 m for sea-ice, 0.10 m for glacial-ice), then set at
!           lower bound and store the source increment in subsurface runoff/
!           baseflow (runoff2).  note:  runoff2 is then a negative value (as
!           a flag) over sea-ice or glacial-ice, in order to achieve water balance.

      elseif (ice == 1) then          ! for sea-ice

        if (sneqv >= 0.01) then

          call snowpack                                                 &
!  ---  inputs:
     &     ( sneqv, dt, t1, yy,                                         &
!  ---  input/outputs:
     &       snowh, sndens                                              &
     &     )

        else

!         sndens = sneqv / snowh
!         runoff2 = -(0.01 - sneqv) / dt
          sneqv = 0.01
          snowh = 0.05
          sncovr = 1.0
!         snowh = sneqv / sndens

        endif   ! end if_sneqv_block

      else                            ! for glacial-ice

        if (sneqv >= 0.10) then

          call snowpack                                                 &
!  ---  inputs:
     &     ( sneqv, dt, t1, yy,                                         &
!  ---  input/outputs:
     &       snowh, sndens                                              &
     &     )

        else

!         sndens = sneqv / snowh
!         runoff2 = -(0.10 - sneqv) / dt
          sneqv = 0.10
          snowh = 0.50
          sncovr = 1.0
!         snowh = sneqv / sndens

        endif   ! end if_sneqv_block

      endif   ! end if_ice_block

!
      return
!...................................
      end subroutine snopac
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates snow depth and densitity to account
!! for the new snowfall. new values of snow depth & density returned.
      subroutine snow_new
!...................................
!  ---  inputs:
!    &     ( sfctmp, sn_new,                                            &
!  ---  input/outputs:
!    &       snowh, sndens                                              &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!    subroutine snow_new calculates snow depth and densitity to account !
!    for the new snowfall. new values of snow depth & density returned. !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from the calling program:                              size   !
!     sfctmp   - real, surface air temperature (k)                 1    !
!     sn_new   - real, new snowfall (m)                            1    !
!                                                                       !
!  input/outputs from and to the calling program:                       !
!     snowh    - real, snow depth (m)                              1    !
!     sndens   - real, snow density                                1    !
!                      (g/cm3=dimensionless fraction of h2o density)    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
!     real(kind=kind_phys), intent(in) :: sfctmp, sn_new

!  ---  input/outputs:
!     real(kind=kind_phys), intent(inout) :: snowh, sndens

!  ---  locals:
      real(kind=kind_phys) :: dsnew, snowhc, hnewc, newsnc, tempc

!
!===> ...  begin here
!
!  --- ...  conversion into simulation units

      snowhc = snowh * 100.0
      newsnc = sn_new * 100.0
      tempc  = sfctmp - tfreez

!  --- ...  calculating new snowfall density depending on temperature
!           equation from gottlib l. 'a general runoff model for
!           snowcovered and glacierized basin', 6th nordic hydrological
!           conference, vemadolen, sweden, 1980, 172-177pp.

      if (tempc <= -15.0) then
        dsnew = 0.05
      else
        dsnew = 0.05 + 0.0017*(tempc + 15.0)**1.5
      endif

!  --- ...  adjustment of snow density depending on new snowfall

      hnewc  = newsnc / dsnew
      sndens = (snowhc*sndens + hnewc*dsnew) / (snowhc + hnewc)
      snowhc = snowhc + hnewc
      snowh  = snowhc * 0.01
!
      return
!...................................
      end subroutine snow_new
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates total roughness length over snow.
      subroutine snowz0
!...................................
!  ---  inputs:
!    &     ( sncovr,                                                    &
!  ---  input/outputs:
!    &       z0                                                         &
!    &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!    subroutine snowz0 calculates total roughness length over snow      !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from the calling program:                              size   !
!     sncovr   - real, fractional snow cover                       1    !
!                                                                       !
!  input/outputs from and to the calling program:                       !
!     z0       - real, roughness length (m)                        1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
!     real(kind=kind_phys), intent(in) :: sncovr

!  ---  input/outputs:
!     real(kind=kind_phys), intent(inout) :: z0

!  ---  locals:
      real(kind=kind_phys) :: z0s
!
!===> ...  begin here
!
!     z0s = 0.001                     ! snow roughness length:=0.001 (m)
!  --- ...  current noah lsm condition - mbek, 09-oct-2001
      z0s = z0

      z0 = (1.0 - sncovr)*z0 + sncovr*z0s

!
      return
!...................................
      end subroutine snowz0
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates thermal diffusivity and conductivity
!! of the soil for a given point and time.
      subroutine tdfcnd                                                 &
!...................................
!  ---  inputs:
     &     ( smc, qz, smcmax, sh2o,                                     &
!  ---  outputs:
     &       df                                                         &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!    subroutine tdfcnd calculates thermal diffusivity and conductivity  !
!    of the soil for a given point and time.                            !
!                                                                       !
!    peters-lidard approach (peters-lidard et al., 1998)                !
!    june 2001 changes: frozen soil condition.                          !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!  use as in peters-lidard, 1998 (modif. from johansen, 1975).          !
!                                 pablo grunmann, 08/17/98              !
!  refs.:                                                               !
!    farouki, o.t.,1986: thermal properties of soils. series on rock    !
!             and soil mechanics, vol. 11, trans tech, 136 pp.          !
!    johansen, o., 1975: thermal conductivity of soils. ph.d. thesis,   !
!             university of trondheim,                                  !
!    peters-lidard, c. d., et al., 1998: the effect of soil thermal     !
!             conductivity parameterization on surface energy fluxes    !
!             and temperatures. journal of the atmospheric sciences,    !
!             vol. 55, pp. 1209-1224.                                   !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     smc      - real, top layer total soil moisture               1    !
!     qz       - real, quartz content (soil type dependent)        1    !
!     smcmax   - real, porosity                                    1    !
!     sh2o     - real, top layer unfrozen soil moisture            1    !
!                                                                       !
!  outputs:                                                             !
!     df       - real, soil thermal diffusivity and conductivity   1    !
!                                                                       !
!  locals:                                                              !
!     thkw     - water thermal conductivity                        1    !
!     thkqtz   - thermal conductivity for quartz                   1    !
!     thko     - thermal conductivity for other soil components    1    !
!     thkqtz   - thermal conductivity for the solids combined      1    !
!     thkice   - ice thermal conductivity                          1    !
!                                                                       !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  input:
      real (kind=kind_phys), intent(in) :: smc, qz, smcmax, sh2o

!  ---  output:
      real (kind=kind_phys), intent(out) :: df

!  ---  locals:
      real (kind=kind_phys) :: gammd, thkdry, ake, thkice, thko,        &
     &       thkqtz, thksat, thks, thkw, satratio, xu, xunfroz
!
!===> ...  begin here
!
!  --- ...  if the soil has any moisture content compute a partial sum/product
!           otherwise use a constant value which works well with most soils

!  --- ...  saturation ratio:

      satratio = smc / smcmax

!  --- ...  parameters  w/(m.k)
      thkice = 2.2
      thkw   = 0.57
      thko   = 2.0
!     if (qz <= 0.2) thko = 3.0
      thkqtz = 7.7

!  --- ...  solids' conductivity

      thks = (thkqtz**qz) * (thko**(1.0-qz))

!  --- ...  unfrozen fraction (from 1., i.e., 100%liquid, to 0. (100% frozen))

      xunfroz = (sh2o + 1.e-9) / (smc + 1.e-9)

!  --- ...  unfrozen volume for saturation (porosity*xunfroz)

      xu=xunfroz*smcmax

!  --- ...  saturated thermal conductivity

      thksat = thks**(1.-smcmax) * thkice**(smcmax-xu) * thkw**(xu)

!  --- ...  dry density in kg/m3

      gammd = (1.0 - smcmax) * 2700.0

!  --- ...  dry thermal conductivity in w.m-1.k-1

      thkdry = (0.135*gammd + 64.7) / (2700.0 - 0.947*gammd)

      if ( sh2o+0.0005 < smc ) then         ! frozen

        ake = satratio

      else                                  ! unfrozen

!  --- ...  range of validity for the kersten number (ake)
        if ( satratio > 0.1 ) then

!  --- ...  kersten number (using "fine" formula, valid for soils containing
!           at least 5% of particles with diameter less than 2.e-6 meters.)
!           (for "coarse" formula, see peters-lidard et al., 1998).

          ake = log10( satratio ) + 1.0

        else

!  --- ...  use k = kdry
          ake = 0.0

        endif   ! end if_satratio_block

      endif   ! end if_sh2o+0.0005_block

!  --- ...  thermal conductivity

      df = ake * (thksat - thkdry) + thkdry
!
      return
!...................................
      end subroutine tdfcnd
!-----------------------------------


!*********************************************!
!  section-2  2nd level subprograms           !
!*********************************************!


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates soil moisture flux. The soil moisture
!! content (smc - a per unit volume measurement) is a dependent variable
!! that is updated with prognostic equations. The canopy moisture content
!! (cmc) is also updated. Frozen ground version: new states added: sh2o,
!! and frozen ground correction factor, frzfact and paramter slope.
      subroutine evapo                                                  &
!...................................
!  ---  inputs:
     &     ( nsoil, nroot, cmc, cmcmax, etp1, dt, zsoil,                &
     &       sh2o, smcmax, smcwlt, smcref, smcdry, pc,                  &
     &       shdfac, cfactr, rtdis, fxexp,                              &
!  ---  outputs:
     &       eta1, edir1, ec1, et1, ett1                                &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine evapo calculates soil moisture flux.  the soil moisture   !
!  content (smc - a per unit volume measurement) is a dependent variable!
!  that is updated with prognostic eqns. the canopy moisture content    !
!  (cmc) is also updated. frozen ground version:  new states added:     !
!  sh2o, and frozen ground correction factor, frzfact and parameter     !
!  slope.                                                               !
!                                                                       !
!                                                                       !
!  subprogram called:  devap, transp                                    !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs from calling program:                                  size   !
!     nsoil    - integer, number of soil layers                    1    !
!     nroot    - integer, number of root layers                    1    !
!     cmc      - real, canopy moisture content                     1    !
!     cmcmax   - real, maximum canopy water parameters             1    !
!     etp1     - real, potential evaporation                       1    !
!     dt       - real, time step                                   1    !
!     zsoil    - real, soil layer depth below ground             nsoil  !
!     sh2o     - real, unfrozen soil moisture                    nsoil  !
!     smcmax   - real, porosity                                    1    !
!     smcwlt   - real, wilting point                               1    !
!     smcref   - real, soil mois threshold                         1    !
!     smcdry   - real, dry soil mois threshold                     1    !
!     pc       - real, plant coeff                                 1    !
!     cfactr   - real, canopy water parameters                     1    !
!     rtdis    - real, root distribution                         nsoil  !
!     fxexp    - real, bare soil evaporation exponent              1    !
!                                                                       !
!  outputs to calling program:                                          !
!     eta1     - real, latent heat flux                            1    !
!     edir1    - real, direct soil evaporation                     1    !
!     ec1      - real, canopy water evaporation                    1    !
!     et1      - real, plant transpiration                       nsoil  !
!     ett1     - real, total plant transpiration                   1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
      integer, intent(in) :: nsoil, nroot

      real (kind=kind_phys),  intent(in) :: cmc, cmcmax, etp1, dt, pc,  &
     &       smcmax, smcwlt, smcref, smcdry, shdfac, cfactr, fxexp,     &
     &       zsoil(nsoil), sh2o(nsoil), rtdis(nsoil)

!  ---  outputs:
      real (kind=kind_phys),  intent(out) :: eta1, edir1, ec1, ett1,    &
     &       et1(nsoil)

!  ---  locals:
      real (kind=kind_phys) :: cmc2ms

      integer :: i, k

!
!===> ...  begin here
!
!  --- ...  executable code begins here if the potential evapotranspiration
!           is greater than zero.

      edir1 = 0.0
      ec1 = 0.0

      do k = 1, nsoil
        et1(k) = 0.0
      enddo
      ett1 = 0.0

      if (etp1 > 0.0) then

!  --- ...  retrieve direct evaporation from soil surface.  call this function
!           only if veg cover not complete.
!           frozen ground version:  sh2o states replace smc states.

        if (shdfac < 1.0) then

          call devap                                                    &
!  ---  inputs:
     &     ( etp1, sh2o(1), shdfac, smcmax, smcdry, fxexp,              &
!  ---  outputs:
     &       edir1                                                      &
     &     )

        endif

!  --- ...  initialize plant total transpiration, retrieve plant transpiration,
!           and accumulate it for all soil layers.

        if (shdfac > 0.0) then

          call transp                                                   &
!  ---  inputs:
     &     ( nsoil, nroot, etp1, sh2o, smcwlt, smcref,                   &
     &       cmc, cmcmax, zsoil, shdfac, pc, cfactr, rtdis,             &
!  ---  outputs:
     &       et1                                                        &
     &     )

          do k = 1, nsoil
            ett1 = ett1 + et1(k)
          enddo

!  --- ...  calculate canopy evaporation.
!           if statements to avoid tangent linear problems near cmc=0.0.

          if (cmc > 0.0) then
            ec1 = shdfac * ( (cmc/cmcmax)**cfactr ) * etp1
          else
            ec1 = 0.0
          endif

!  --- ...  ec should be limited by the total amount of available water
!           on the canopy.  -f.chen, 18-oct-1994

          cmc2ms = cmc / dt
          ec1 = min ( cmc2ms, ec1 )
        endif

      endif   ! end if_etp1_block

!  --- ...  total up evap and transp types to obtain actual evapotransp

      eta1 = edir1 + ett1 + ec1

!
      return
!...................................
      end subroutine evapo
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine updates the temperature state of the soil column
!! based on the thermal diffusion equation and update the frozen soil
!! moisture content based on the temperature.
      subroutine shflx                                                  &
!...................................
!  ---  inputs:
     &     ( nsoil, smc, smcmax, dt, yy, zz1, zsoil, zbot,              &
     &       psisat, bexp, df1, ice, quartz, csoil, vegtyp,             &
!  ---  input/outputs:
     &       stc, t1, tbot, sh2o,                                       &
!  ---  outputs:
     &       ssoil                                                      &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine shflx updates the temperature state of the soil column    !
!  based on the thermal diffusion equation and update the frozen soil   !
!  moisture content based on the temperature.                           !
!                                                                       !
!  subprogram called:  hstep, hrtice, hrt                               !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     smc      - real, total soil moisture                       nsoil  !
!     smcmax   - real, porosity (sat val of soil mois)             1    !
!     dt       - real, time step                                   1    !
!     yy       - real, soil temperature at the top of column       1    !
!     zz1      - real,                                             1    !
!     zsoil    - real, soil layer depth below ground (negative)  nsoil  !
!     zbot     - real, specify depth of lower bd soil              1    !
!     psisat   - real, saturated soil potential                    1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     df1      - real, thermal diffusivity and conductivity        1    !
!     ice      - integer, sea-ice flag (=1: sea-ice, =0: land)     1    !
!     quartz   - real, soil quartz content                         1    !
!     csoil    - real, soil heat capacity                          1    !
!     vegtyp   - integer, vegtation type                           1    !
!                                                                       !
!  input/outputs:                                                       !
!     stc      - real, soil temp                                 nsoil  !
!     t1       - real, ground/canopy/snowpack eff skin temp        1    !
!     tbot     - real, bottom soil temp                            1    !
!     sh2o     - real, unfrozen soil moisture                    nsoil  !
!                                                                       !
!  outputs:                                                             !
!     ssoil    - real, upward soil heat flux                       1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  parameter constants:
      real (kind=kind_phys), parameter :: ctfil1 = 0.5
      real (kind=kind_phys), parameter :: ctfil2 = 1.0 - ctfil1

!  ---  inputs:
      integer, intent(in) :: nsoil, ice, vegtyp

      real (kind=kind_phys), intent(in) :: smc(nsoil), smcmax, dt, yy,  &
     &       zz1, zsoil(nsoil), zbot, psisat, bexp, df1, quartz, csoil

!  ---  input/outputs:
      real (kind=kind_phys), intent(inout) :: stc(nsoil), t1, tbot,     &
     &       sh2o(nsoil)

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: ssoil

!  ---  locals:
      real (kind=kind_phys) :: ai(nsold), bi(nsold), ci(nsold), oldt1,  &
     &       rhsts(nsold), stcf(nsold), stsoil(nsoil)

      integer :: i

!
!===> ...  begin here
!
      oldt1 = t1
      do i = 1, nsoil
         stsoil(i) = stc(i)
      enddo

!  --- ...  hrt routine calcs the right hand side of the soil temp dif eqn

      if (ice /= 0) then

!  --- ...  sea-ice case, glacial-ice case

        call hrtice                                                     &
!  ---  inputs:
     &     ( nsoil, stc, zsoil, yy, zz1, df1, ice,                      &
!  ---  input/outputs:
     &       tbot,                                                      &
!  ---  outputs:
     &       rhsts, ai, bi, ci                                          &
     &     )

        call hstep                                                      &
!  ---  inputs:
     &     ( nsoil, stc, dt,                                            &
!  ---  input/outputs:
     &       rhsts, ai, bi, ci,                                         &
!  ---  outputs:
     &       stcf                                                       &
     &     )

      else

!  --- ...  land-mass case

        call hrt                                                        &
!  ---  inputs:
     &     ( nsoil, stc, smc, smcmax, zsoil, yy, zz1, tbot,             &
     &       zbot, psisat, dt, bexp, df1, quartz, csoil,vegtyp,         &
!  ---  input/outputs:
     &       sh2o,                                                      &
!  ---  outputs:
     &       rhsts, ai, bi, ci                                          &
     &     )

        call hstep                                                      &
!  ---  inputs:
     &     ( nsoil, stc, dt,                                            &
!  ---  input/outputs:
     &       rhsts, ai, bi, ci,                                         &
!  ---  outputs:
     &       stcf                                                       &
     &     )

      endif

      do i = 1, nsoil
         stc(i) = stcf(i)
      enddo

!  --- ...  in the no snowpack case (via routine nopac branch,) update the grnd
!           (skin) temperature here in response to the updated soil temperature
!           profile above.  (note: inspection of routine snopac shows that t1
!           below is a dummy variable only, as skin temperature is updated
!           differently in routine snopac)

      t1 = (yy + (zz1 - 1.0)*stc(1)) / zz1
      t1 = ctfil1*t1 + ctfil2*oldt1

      do i = 1, nsoil
        stc(i) = ctfil1*stc(i) + ctfil2*stsoil(i)
      enddo

!  --- ...  calculate surface soil heat flux

      ssoil = df1*(stc(1) - t1) / (0.5*zsoil(1))

!
      return
!...................................
      end subroutine shflx
!-----------------------------------



!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates soil moisture flux. The soil moisture
!! content (smc - a per unit vulume measurement) is a dependent variable
!! that is updated with prognostic equations. The canopy moisture content
!! (cmc) is also updated. Frozen ground version: new states added: sh2o and
!! frozen ground correction factor, frzx and parameter slope.
      subroutine smflx                                                    &
!...................................
!  ---  inputs:
     &     ( nsoil, dt, kdt, smcmax, smcwlt, cmcmax, prcp1,               &
     &       zsoil, slope, frzx, bexp, dksat, dwsat, shdfac,              &
     &       edir1, ec1, et1,                                             &
!  ---  input/outputs:
     &       cmc, sh2o,                                                   &
!  ---  outputs:
     &       smc, runoff1, runoff2, runoff3, drip                         &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine smflx calculates soil moisture flux.  the soil moisture   !
!  content (smc - a per unit volume measurement) is a dependent variable!
!  that is updated with prognostic eqns. the canopy moisture content    !
!  (cmc) is also updated. frozen ground version:  new states added: sh2o!
!  and frozen ground correction factor, frzx and parameter slope.       !
!                                                                       !
!                                                                       !
!  subprogram called:  srt, sstep                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     dt       - real, time step                                   1    !
!     kdt      - real,                                             1    !
!     smcmax   - real, porosity                                    1    !
!     smcwlt   - real, wilting point                               1    !
!     cmcmax   - real, maximum canopy water parameters             1    !
!     prcp1    - real, effective precip                            1    !
!     zsoil    - real, soil layer depth below ground (negative)  nsoil  !
!     slope    - real, linear reservoir coefficient                1    !
!     frzx     - real, frozen ground parameter                     1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     dksat    - real, saturated soil hydraulic conductivity       1    !
!     dwsat    - real, saturated soil diffusivity                  1    !
!     shdfac   - real, aeral coverage of green veg                 1    !
!     edir1    - real, direct soil evaporation                     1    !
!     ec1      - real, canopy water evaporation                    1    !
!     et1      - real, plant transpiration                       nsoil  !
!                                                                       !
!  input/outputs:                                                       !
!     cmc      - real, canopy moisture content                     1    !
!     sh2o     - real, unfrozen soil moisture                    nsoil  !
!                                                                       !
!  outputs:                                                             !
!     smc      - real, total soil moisture                       nsoil  !
!     runoff1  - real, surface runoff not infiltrating sfc         1    !
!     runoff2  - real, sub surface runoff (baseflow)               1    !
!     runoff3  - real, excess of porosity                          1    !
!     drip     - real, through-fall of precip and/or dew           1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
      integer, intent(in) :: nsoil

      real (kind=kind_phys),  intent(in) :: dt, kdt, smcmax, smcwlt,    &
     &       cmcmax, prcp1, slope, frzx, bexp, dksat, dwsat, shdfac,    &
     &       edir1, ec1, et1(nsoil), zsoil(nsoil)

!  ---  input/outputs:
      real (kind=kind_phys),  intent(inout) :: cmc, sh2o(nsoil)

!  ---  outputs:
      real (kind=kind_phys),  intent(out) :: smc(nsoil), runoff1,       &
     &       runoff2, runoff3, drip

!  ---  locals:
      real (kind=kind_phys) :: dummy, excess, pcpdrp, rhsct, trhsct,    &
     &       rhstt(nsold), sice(nsold), sh2oa(nsold), sh2ofg(nsold),    &
     &       ai(nsold), bi(nsold), ci(nsold)

      integer :: i, k
!
!===> ...  begin here
!
!  --- ...  executable code begins here.

      dummy = 0.0

!  --- ...  compute the right hand side of the canopy eqn term ( rhsct )

      rhsct = shdfac*prcp1 - ec1

!  --- ...  convert rhsct (a rate) to trhsct (an amount) and add it to
!           existing cmc.  if resulting amt exceeds max capacity, it becomes
!           drip and will fall to the grnd.

      drip = 0.0
      trhsct = dt * rhsct
      excess = cmc + trhsct

      if (excess > cmcmax) drip = excess - cmcmax

!  --- ...  pcpdrp is the combined prcp1 and drip (from cmc) that goes into
!           the soil

      pcpdrp = (1.0 - shdfac)*prcp1 + drip/dt

!  --- ...  store ice content at each soil layer before calling srt & sstep

      do i = 1, nsoil
        sice(i) = smc(i) - sh2o(i)
      enddo

!  --- ...  call subroutines srt and sstep to solve the soil moisture
!           tendency equations.

!  ---  if the infiltrating precip rate is nontrivial,
!         (we consider nontrivial to be a precip total over the time step
!         exceeding one one-thousandth of the water holding capacity of
!         the first soil layer)
!       then call the srt/sstep subroutine pair twice in the manner of
!         time scheme "f" (implicit state, averaged coefficient)
!         of section 2 of kalnay and kanamitsu (1988, mwr, vol 116,
!         pages 1945-1958)to minimize 2-delta-t oscillations in the
!         soil moisture value of the top soil layer that can arise because
!         of the extreme nonlinear dependence of the soil hydraulic
!         diffusivity coefficient and the hydraulic conductivity on the
!         soil moisture state
!       otherwise call the srt/sstep subroutine pair once in the manner of
!         time scheme "d" (implicit state, explicit coefficient)
!         of section 2 of kalnay and kanamitsu
!       pcpdrp is units of kg/m**2/s or mm/s, zsoil is negative depth in m

!     if ( pcpdrp .gt. 0.0 ) then
      if ( (pcpdrp*dt) > (0.001*1000.0*(-zsoil(1))*smcmax) ) then

!  --- ...  frozen ground version:
!           smc states replaced by sh2o states in srt subr.  sh2o & sice states
!           included in sstep subr.  frozen ground correction factor, frzx
!           added.  all water balance calculations using unfrozen water

        call srt                                                        &
!  ---  inputs:
     &     ( nsoil, edir1, et1, sh2o, sh2o, pcpdrp, zsoil, dwsat,       &
     &       dksat, smcmax, bexp, dt, smcwlt, slope, kdt, frzx, sice,   &
!  ---  outputs:
     &       rhstt, runoff1, runoff2, ai, bi, ci                        &
     &     )

        call sstep                                                      &
!  ---  inputs:
     &     ( nsoil, sh2o, rhsct, dt, smcmax, cmcmax, zsoil, sice,       &
!  ---  input/outputs:
     &       dummy, rhstt, ai, bi, ci,                                  &
!  ---  outputs:
     &       sh2ofg, runoff3, smc                                       &
     &     )

        do k = 1, nsoil
          sh2oa(k) = (sh2o(k) + sh2ofg(k)) * 0.5
        enddo

        call srt                                                        &
!  ---  inputs:
     &     ( nsoil, edir1, et1, sh2o, sh2oa, pcpdrp, zsoil, dwsat,      &
     &       dksat, smcmax, bexp, dt, smcwlt, slope, kdt, frzx, sice,   &
!  ---  outputs:
     &       rhstt, runoff1, runoff2, ai, bi, ci                        &
     &     )

        call sstep                                                      &
!  ---  inputs:
     &     ( nsoil, sh2o, rhsct, dt, smcmax, cmcmax, zsoil, sice,       &
!  ---  input/outputs:
     &       cmc, rhstt, ai, bi, ci,                                    &
!  ---  outputs:
     &       sh2o, runoff3, smc                                         &
     &     )

      else

        call srt                                                        &
!  ---  inputs:
     &     ( nsoil, edir1, et1, sh2o, sh2o, pcpdrp, zsoil, dwsat,       &
     &       dksat, smcmax, bexp, dt, smcwlt, slope, kdt, frzx, sice,   &
!  ---  outputs:
     &       rhstt, runoff1, runoff2, ai, bi, ci                        &
     &     )

        call sstep                                                      &
!  ---  inputs:
     &     ( nsoil, sh2o, rhsct, dt, smcmax, cmcmax, zsoil, sice,       &
!  ---  input/outputs:
     &       cmc, rhstt, ai, bi, ci,                                    &
!  ---  outputs:
     &       sh2o, runoff3, smc                                         &
     &     )

      endif

!     runof = runoff
!
      return
!...................................
      end subroutine smflx
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates compaction of a snowpack under conditions of 
!! increasing snow density, as obtained from an approximate solution of 
!! E. Anderson's differential equation (3.29),NOAA technical report NWS 19,
!! by Victor Koren, 03/25/95. subroutine will return new values of \a snowh
!! and \a sndens .
      subroutine snowpack                                               &
!...................................
!  ---  inputs:
     &     ( esd, dtsec, tsnow, tsoil,                                  &
!  ---  input/outputs:
     &       snowh, sndens                                              &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!    subroutine snowpack calculates compaction of snowpack under        !
!    conditions of increasing snow density, as obtained from an         !
!    approximate solution of e. anderson's differential equation (3.29),!
!    noaa technical report nws 19, by victor koren, 03/25/95.           !
!    subroutine will return new values of snowh and sndens              !
!                                                                       !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     esd      - real, water equivalent of snow (m)                1    !
!     dtsec    - real, time step (sec)                             1    !
!     tsnow    - real, snow surface temperature (k)                1    !
!     tsoil    - real, soil surface temperature (k)                1    !
!                                                                       !
!  input/outputs:                                                       !
!     snowh    - real, snow depth (m)                              1    !
!     sndens   - real, snow density                                1    !
!                      (g/cm3=dimensionless fraction of h2o density)    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  parameter constants:
      real (kind=kind_phys), parameter :: c1 = 0.01
      real (kind=kind_phys), parameter :: c2 = 21.0

!  ---  inputs:
      real (kind=kind_phys), intent(in) :: esd, dtsec, tsnow, tsoil

!  ---  input/outputs:
      real (kind=kind_phys), intent(inout) :: snowh, sndens

!  ---  locals:
      real (kind=kind_phys) :: bfac, dsx, dthr, dw, snowhc, pexp,       &
     &       tavgc, tsnowc, tsoilc, esdc, esdcx

      integer :: ipol, j
!
!===> ...  begin here
!
!  --- ...  conversion into simulation units

      snowhc = snowh * 100.0
      esdc   = esd * 100.0
      dthr   = dtsec / 3600.0
      tsnowc = tsnow - tfreez
      tsoilc = tsoil - tfreez

!  --- ...  calculating of average temperature of snow pack

      tavgc = 0.5 * (tsnowc + tsoilc)

!  --- ...  calculating of snow depth and density as a result of compaction
!           sndens=ds0*(exp(bfac*esd)-1.)/(bfac*esd)
!           bfac=dthr*c1*exp(0.08*tavgc-c2*ds0)
!     note: bfac*esd in sndens eqn above has to be carefully treated
!           numerically below:
!        c1 is the fractional increase in density (1/(cm*hr))
!        c2 is a constant (cm3/g) kojima estimated as 21 cms/g

      if (esdc > 1.e-2) then
        esdcx = esdc
      else
        esdcx = 1.e-2
      endif

      bfac = dthr*c1 * exp(0.08*tavgc - c2*sndens)

!     dsx = sndens * ((dexp(bfac*esdc)-1.0) / (bfac*esdc))

!  --- ...  the function of the form (e**x-1)/x imbedded in above expression
!           for dsx was causing numerical difficulties when the denominator "x"
!           (i.e. bfac*esdc) became zero or approached zero (despite the fact
!           that the analytical function (e**x-1)/x has a well defined limit
!           as "x" approaches zero), hence below we replace the (e**x-1)/x
!           expression with an equivalent, numerically well-behaved
!           polynomial expansion.

!  --- ...  number of terms of polynomial expansion, and hence its accuracy,
!           is governed by iteration limit "ipol".
!           ipol greater than 9 only makes a difference on double
!           precision (relative errors given in percent %).
!       ipol=9, for rel.error <~ 1.6 e-6 % (8 significant digits)
!       ipol=8, for rel.error <~ 1.8 e-5 % (7 significant digits)
!       ipol=7, for rel.error <~ 1.8 e-4 % ...

      ipol = 4
      pexp = 0.0

      do j = ipol, 1, -1
!       pexp = (1.0 + pexp)*bfac*esdc /real(j+1)
        pexp = (1.0 + pexp)*bfac*esdcx/real(j+1)
      enddo
      pexp = pexp + 1.

      dsx = sndens * pexp

!  --- ...  above line ends polynomial substitution
!           end of koren formulation

!! --- ...  base formulation (cogley et al., 1990)
!           convert density from g/cm3 to kg/m3

!!      dsm = sndens * 1000.0

!!      dsx = dsm + dtsec*0.5*dsm*gs2*esd /                             &
!!   &        (1.e7*exp(-0.02*dsm + kn/(tavgc+273.16)-14.643))

!! --- ...  convert density from kg/m3 to g/cm3

!!      dsx = dsx / 1000.0

!! --- ...  end of cogley et al. formulation

!  --- ...  set upper/lower limit on snow density

      dsx = max( min( dsx, 0.40 ), 0.05 )
      sndens = dsx

!  --- ...  update of snow depth and density depending on liquid water
!           during snowmelt.  assumed that 13% of liquid water can be
!           stored in snow per day during snowmelt till snow density 0.40.

      if (tsnowc >= 0.0) then
        dw = 0.13 * dthr / 24.0
        sndens = sndens*(1.0 - dw) + dw
        if (sndens > 0.40) sndens = 0.40
      endif

!  --- ...  calculate snow depth (cm) from snow water equivalent and snow
!           density. change snow depth units to meters

      snowhc = esdc / sndens
      snowh  = snowhc * 0.01

!
      return
!...................................
      end subroutine snowpack
!-----------------------------------


!*********************************************!
!  section-3  3rd or lower level subprograms  !
!*********************************************!


!-----------------------------------
!>\ingroup Noah_LSM
!> This subrtouine calculates direct soil evaporation.
      subroutine devap                                                  &
!...................................
!  ---  inputs:
     &     ( etp1, smc, shdfac, smcmax, smcdry, fxexp,                  &
!  ---  outputs:
     &       edir1                                                      &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine devap calculates direct soil evaporation                  !
!                                                                       !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     etp1     - real, potential evaporation                       1    !
!     smc      - real, unfrozen soil moisture                      1    !
!     shdfac   - real, aeral coverage of green vegetation          1    !
!     smcmax   - real, porosity (sat val of soil mois)             1    !
!     smcdry   - real, dry soil mois threshold                     1    !
!     fxexp    - real, bare soil evaporation exponent              1    !
!                                                                       !
!  outputs:                                                             !
!     edir1    - real, direct soil evaporation                     1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
      real (kind=kind_phys), intent(in) :: etp1, smc, shdfac, smcmax,   &
     &       smcdry, fxexp

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: edir1

!  ---  locals:
      real (kind=kind_phys) :: fx, sratio
!
!===> ...  begin here
!
!  --- ...  direct evap a function of relative soil moisture availability,
!           linear when fxexp=1.
!           fx > 1 represents demand control
!           fx < 1 represents flux control

      sratio = (smc - smcdry) / (smcmax - smcdry)

      if (sratio > 0.0) then
        fx = sratio**fxexp
        fx = max ( min ( fx, 1.0 ), 0.0 )
      else
        fx = 0.0
      endif

!  --- ...  allow for the direct-evap-reducing effect of shade

      edir1 = fx * ( 1.0 - shdfac ) * etp1
!
      return
!...................................
      end subroutine devap
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates amount of supercooled liquid soil water
!! content if temperature is below 273.15K (t0). It requires Newton-type
!! iteration to solve the nonlinear implicit equation given in eqn 17
!! of \cite koren_et_al_1999.
!!
!! New version (June 2001): much faster and more accurate Newton iteration
!! achieved by first taking log of eqn cited above -- less than 4 (typically
!! 1 or 2) iterations achieves convergence. Also, explicit 1-step solution
!! option for special case of paramter ck=0, which reduces the orginal 
!! implicit equation to a simpler explicit form, known as the "flerchinger eqn".
!! Improved handling of solution in the limit of freezing point temperature t0.
      subroutine frh2o                                                  &
!...................................
!  ---  inputs:
     &     ( tkelv, smc, sh2o, smcmax, bexp, psis,                      &
!  ---  outputs:
     &       liqwat                                                     &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine frh2o calculates amount of supercooled liquid soil water  !
!  content if temperature is below 273.15k (t0).  requires newton-type  !
!  iteration to solve the nonlinear implicit equation given in eqn 17   !
!  of koren et al (1999, jgr, vol 104(d16), 19569-19585).               !
!                                                                       !
!  new version (june 2001): much faster and more accurate newton        !
!  iteration achieved by first taking log of eqn cited above -- less    !
!  than 4 (typically 1 or 2) iterations achieves convergence.  also,    !
!  explicit 1-step solution option for special case of parameter ck=0,  !
!  which reduces the original implicit equation to a simpler explicit   !
!  form, known as the "flerchinger eqn". improved handling of solution  !
!  in the limit of freezing point temperature t0.                       !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     tkelv    - real, temperature (k)                             1    !
!     smc      - real, total soil moisture content (volumetric)    1    !
!     sh2o     - real, liquid soil moisture content (volumetric)   1    !
!     smcmax   - real, saturation soil moisture content            1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     psis     - real, saturated soil matric potential             1    !
!                                                                       !
!  outputs:                                                             !
!     liqwat   - real, supercooled liquid water content            1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  constant parameters:
      real (kind=kind_phys), parameter :: ck    = 8.0
!     real (kind=kind_phys), parameter :: ck    = 0.0
      real (kind=kind_phys), parameter :: blim  = 5.5
      real (kind=kind_phys), parameter :: error = 0.005

!  ---  inputs:
      real (kind=kind_phys), intent(in) :: tkelv, smc, sh2o, smcmax,    &
     &       bexp, psis

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: liqwat

!  ---  locals:
      real (kind=kind_phys) :: bx, denom, df, dswl, fk, swl, swlk

      integer :: nlog, kcount
!
!===> ...  begin here
!
!  --- ...  limits on parameter b: b < 5.5  (use parameter blim)
!           simulations showed if b > 5.5 unfrozen water content is
!           non-realistically high at very low temperatures.

      bx = bexp
      if (bexp > blim)  bx = blim

!  --- ...  initializing iterations counter and iterative solution flag.

      nlog  = 0
      kcount= 0

!  --- ...  if temperature not significantly below freezing (t0), sh2o = smc

      if (tkelv > (tfreez-1.e-3)) then

        liqwat = smc

      else

        if (ck /= 0.0) then

!  --- ...  option 1: iterated solution for nonzero ck
!                     in koren et al, jgr, 1999, eqn 17

!  --- ...  initial guess for swl (frozen content)

          swl = smc - sh2o

!  --- ...  keep within bounds.

          swl = max( min( swl, smc-0.02 ), 0.0 )

!  --- ...  start of iterations

          do while ( (nlog < 10) .and. (kcount == 0) )
            nlog = nlog + 1

            df = alog( (psis*gs2/lsubf) * ( (1.0 + ck*swl)**2.0 )       &
     &         * (smcmax/(smc-swl))**bx ) - alog(-(tkelv-tfreez)/tkelv)

            denom = 2.0*ck/(1.0 + ck*swl) + bx/(smc - swl)
            swlk  = swl - df/denom

!  --- ...  bounds useful for mathematical solution.

            swlk = max( min( swlk, smc-0.02 ), 0.0 )

!  --- ...  mathematical solution bounds applied.

            dswl = abs(swlk - swl)
            swl = swlk

!  --- ...  if more than 10 iterations, use explicit method (ck=0 approx.)
!           when dswl less or eq. error, no more iterations required.

            if ( dswl <= error )  then
              kcount = kcount + 1
            endif
          enddo   !  end do_while_loop

!  --- ...  bounds applied within do-block are valid for physical solution.

          liqwat = smc - swl

        endif   ! end if_ck_block

!  --- ...  option 2: explicit solution for flerchinger eq. i.e. ck=0
!                     in koren et al., jgr, 1999, eqn 17
!           apply physical bounds to flerchinger solution

        if (kcount == 0) then
          fk = ( ( (lsubf/(gs2*(-psis)))                                &
     &       * ((tkelv-tfreez)/tkelv) )**(-1/bx) ) * smcmax

          fk = max( fk, 0.02 )

          liqwat = min( fk, smc )
        endif

      endif   ! end if_tkelv_block
!
      return
!...................................
      end subroutine frh2o
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates the right hand side of the time tendency 
!! term of the soil thermal diffusion equation. Also to compute (prepare)
!! the matrix coefficients for the tri-diagonal matrix of the implicit time
!! scheme.
      subroutine hrt                                                    &
!...................................
!  ---  inputs:
     &     ( nsoil, stc, smc, smcmax, zsoil, yy, zz1, tbot,             &
     &       zbot, psisat, dt, bexp, df1, quartz, csoil, vegtyp,        &
!  ---  input/outputs:
     &       sh2o,                                                      &
!  ---  outputs:
     &       rhsts, ai, bi, ci                                          &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine hrt calculates the right hand side of the time tendency   !
!  term of the soil thermal diffusion equation.  also to compute        !
!  (prepare) the matrix coefficients for the tri-diagonal matrix of     !
!  the implicit time scheme.                                            !
!                                                                       !
!  subprogram called:  tbnd, snksrc, tmpavg                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     stc      - real, soil temperature                          nsoil  !
!     smc      - real, total soil moisture                       nsoil  !
!     smcmax   - real, porosity                                    1    !
!     zsoil    - real, soil layer depth below ground (negative)  nsoil  !
!     yy       - real,                                             1    !
!     zz1      - real, soil temperture at the top soil column      1    !
!     tbot     - real, bottom soil temp                            1    !
!     zbot     - real, specify depth of lower bd soil              1    !
!     psisat   - real, saturated soil potential                    1    !
!     dt       - real, time step                                   1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     df1      - real, thermal diffusivity                         1    !
!     quartz   - real, soil quartz content                         1    !
!     csoil    - real, soil heat capacity                          1    !
!     vegtyp   - integer, vegetation type                          1    !
!                                                                       !
!  input/outputs:                                                       !
!     sh2o     - real, unfrozen soil moisture                    nsoil  !
!                                                                       !
!  outputs:                                                             !
!     rhsts    - real, time tendency of soil thermal diffusion   nsoil  !
!     ai       - real, matrix coefficients                       nsold  !
!     bi       - real, matrix coefficients                       nsold  !
!     ci       - real, matrix coefficients                       nsold  !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
      integer, intent(in) :: nsoil, vegtyp

      real (kind=kind_phys),  intent(in) :: stc(nsoil), smc(nsoil),     &
     &       smcmax, zsoil(nsoil), yy, zz1, tbot, zbot, psisat, dt,     &
     &       bexp, df1, quartz, csoil

!  ---  input/outputs:
      real (kind=kind_phys),  intent(inout) :: sh2o(nsoil)

!  ---  outputs:
      real (kind=kind_phys),  intent(out) :: rhsts(nsoil), ai(nsold),   &
     &       bi(nsold), ci(nsold)

!  ---  locals:
      real (kind=kind_phys) :: ddz, ddz2, denom, df1n, df1k, dtsdz,     &
     &       dtsdz2, hcpct, qtot, ssoil, sice, tavg, tbk, tbk1,         &
     &       tsnsr, tsurf, csoil_loc

      integer :: i, k

      logical :: itavg

!
!===> ...  begin here
!
        csoil_loc=csoil

       if (ivegsrc == 1)then
!urban
        if( vegtyp == 13 ) then
            csoil_loc=3.0e6
        endif
       endif

!  --- ...  initialize logical for soil layer temperature averaging.

      itavg = .true.
!     itavg = .false.

!  ===  begin section for top soil layer

!  --- ...  calc the heat capacity of the top soil layer

      hcpct = sh2o(1)*cph2o2 + (1.0 - smcmax)*csoil_loc                 &
     &      + (smcmax - smc(1))*cp2 + (smc(1) - sh2o(1))*cpice1

!  --- ...  calc the matrix coefficients ai, bi, and ci for the top layer

      ddz = 1.0 / ( -0.5*zsoil(2) )
      ai(1) = 0.0
      ci(1) = (df1*ddz) / ( zsoil(1)*hcpct )
      bi(1) = -ci(1) + df1 / ( 0.5*zsoil(1)*zsoil(1)*hcpct*zz1 )

!  --- ...  calculate the vertical soil temp gradient btwn the 1st and 2nd soil
!           layers.  then calculate the subsurface heat flux. use the temp
!           gradient and subsfc heat flux to calc "right-hand side tendency
!           terms", or "rhsts", for top soil layer.

      dtsdz = (stc(1) - stc(2)) / (-0.5*zsoil(2))
      ssoil = df1 * (stc(1) - yy) / (0.5*zsoil(1)*zz1)
      rhsts(1) = (df1*dtsdz - ssoil) / (zsoil(1)*hcpct)

!  --- ...  next capture the vertical difference of the heat flux at top and
!           bottom of first soil layer for use in heat flux constraint applied to
!           potential soil freezing/thawing in routine snksrc.

      qtot = ssoil - df1*dtsdz

!  --- ...  if temperature averaging invoked (itavg=true; else skip):
!           set temp "tsurf" at top of soil column (for use in freezing soil
!           physics later in subroutine snksrc).  if snowpack content is
!           zero, then tsurf expression below gives tsurf = skin temp.  if
!           snowpack is nonzero (hence argument zz1=1), then tsurf expression
!           below yields soil column top temperature under snowpack.  then
!           calculate temperature at bottom interface of 1st soil layer for use
!           later in subroutine snksrc

      if (itavg) then

        tsurf = (yy + (zz1-1)*stc(1)) / zz1

        call tbnd                                                       &
!  ---  inputs:
     &     ( stc(1), stc(2), zsoil, zbot, 1, nsoil,                     &
!  ---  outputs:
     &       tbk                                                        &
     &     )

      endif

!  --- ...  calculate frozen water content in 1st soil layer.

      sice = smc(1) - sh2o(1)

!  --- ...  if frozen water present or any of layer-1 mid-point or bounding
!           interface temperatures below freezing, then call snksrc to
!           compute heat source/sink (and change in frozen water content)
!           due to possible soil water phase change

      if ( (sice > 0.0) .or. (tsurf < tfreez) .or.                      &
     &     (stc(1) < tfreez) .or. (tbk < tfreez) ) then

        if (itavg) then

          call tmpavg                                                   &
!  ---  inputs:
     &     ( tsurf, stc(1), tbk, zsoil, nsoil, 1,                       &
!  ---  outputs:
     &       tavg                                                       &
     &     )

        else

          tavg = stc(1)

        endif   ! end if_itavg_block

        call snksrc                                                     &
!  ---  inputs:
     &     ( nsoil, 1, tavg, smc(1), smcmax, psisat, bexp, dt,          &
     &       qtot, zsoil,                                               &
!  ---  input/outputs:
     &       sh2o(1),                                                   &
!  ---  outputs:
     &       tsnsr                                                      &
     &     )


        rhsts(1) = rhsts(1) - tsnsr / ( zsoil(1)*hcpct )

      endif   ! end if_sice_block

!  ===  this ends section for top soil layer.

!  --- ...  initialize ddz2

      ddz2 = 0.0

!  --- ...  loop thru the remaining soil layers, repeating the above process
!           (except subsfc or "ground" heat flux not repeated in lower layers)

      df1k = df1

      do k = 2, nsoil

!  --- ...  calculate heat capacity for this soil layer.

        hcpct = sh2o(k)*cph2o2 + (1.0 - smcmax)*csoil_loc               &
     &        + (smcmax - smc(k))*cp2 + (smc(k) - sh2o(k))*cpice1

        if (k /= nsoil) then

!  --- ...  this section for layer 2 or greater, but not last layer.
!           calculate thermal diffusivity for this layer.

          call tdfcnd                                                   &
!  ---  inputs:
     &     ( smc(k), quartz, smcmax, sh2o(k),                           &
!  ---  outputs:
     &       df1n                                                       &
     &     )
!urban
      if (ivegsrc == 1)then
       if ( vegtyp == 13 ) df1n = 3.24
      endif

!  --- ...  calc the vertical soil temp gradient thru this layer

          denom = 0.5 * (zsoil(k-1) - zsoil(k+1))
          dtsdz2 = (stc(k) - stc(k+1)) / denom

!  --- ...  calc the matrix coef, ci, after calc'ng its partial product

          ddz2 = 2.0 / (zsoil(k-1) - zsoil(k+1))
          ci(k) = -df1n*ddz2 / ((zsoil(k-1) - zsoil(k)) * hcpct)

!  --- ...  if temperature averaging invoked (itavg=true; else skip):
!           calculate temp at bottom of layer.

          if (itavg) then

            call tbnd                                                   &
!  ---  inputs:
     &     ( stc(k), stc(k+1), zsoil, zbot, k, nsoil,                   &
!  ---  outputs:
     &       tbk1                                                       &
     &     )

          endif

        else

!  --- ...  special case of bottom soil layer:  calculate thermal diffusivity
!           for bottom layer.

          call tdfcnd                                                   &
!  ---  inputs:
     &     ( smc(k), quartz, smcmax, sh2o(k),                           &
!  ---  outputs:
     &       df1n                                                       &
     &     )
!urban
      if (ivegsrc == 1)then
       if ( vegtyp == 13 ) df1n = 3.24
      endif

!  --- ...  calc the vertical soil temp gradient thru bottom layer.

          denom = 0.5 * (zsoil(k-1) + zsoil(k)) - zbot
          dtsdz2 = (stc(k) - tbot) / denom

!  --- ...  set matrix coef, ci to zero if bottom layer.

          ci(k) = 0.0

!  --- ...  if temperature averaging invoked (itavg=true; else skip):
!           calculate temp at bottom of last layer.

          if (itavg) then

            call tbnd                                                   &
!  ---  inputs:
     &     ( stc(k), tbot, zsoil, zbot, k, nsoil,                       &
!  ---  outputs:
     &       tbk1                                                       &
     &     )

          endif

        endif   ! end if_k_block

!  --- ...  calculate rhsts for this layer after calc'ng a partial product.

        denom = (zsoil(k) - zsoil(k-1)) * hcpct
        rhsts(k) = ( df1n*dtsdz2 - df1k*dtsdz ) / denom

        qtot = -1.0 * denom * rhsts(k)
        sice = smc(k) - sh2o(k)

        if ( (sice > 0.0) .or. (tbk < tfreez) .or.                      &
     &       (stc(k) < tfreez) .or. (tbk1 < tfreez) ) then

          if (itavg) then

            call tmpavg                                                 &
!  ---  inputs:
     &     ( tbk, stc(k), tbk1, zsoil, nsoil, k,                        &
!  ---  outputs:
     &       tavg                                                       &
     &     )

          else
            tavg = stc(k)
          endif

          call snksrc                                                   &
!  ---  inputs:
     &     ( nsoil, k, tavg, smc(k), smcmax, psisat, bexp, dt,          &
     &       qtot, zsoil,                                               &
!  ---  input/outputs:
     &       sh2o(k),                                                   &
!  ---  outputs:
     &       tsnsr                                                      &
     &     )

          rhsts(k) = rhsts(k) - tsnsr/denom
        endif

!  --- ...  calc matrix coefs, ai, and bi for this layer.

        ai(k) = - df1 * ddz / ((zsoil(k-1) - zsoil(k)) * hcpct)
        bi(k) = -(ai(k) + ci(k))

!  --- ...  reset values of df1, dtsdz, ddz, and tbk for loop to next soil layer.

        tbk   = tbk1
        df1k  = df1n
        dtsdz = dtsdz2
        ddz   = ddz2

      enddo   ! end do_k_loop

!
      return
!...................................
      end subroutine hrt
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates the right hand side of the time tendency
!! term of the soil thermal diffusion equation for sea-ice (ice = 1) or
!! glacial-ice (ice).
      subroutine hrtice                                                 &
!...................................
!  ---  inputs:
     &     ( nsoil, stc, zsoil, yy, zz1, df1, ice,                      &
!  ---  input/outputs:
     &       tbot,                                                      &
!  ---  outputs:
     &       rhsts, ai, bi, ci                                          &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine hrtice calculates the right hand side of the time tendency!
!  term of the soil thermal diffusion equation for sea-ice (ice = 1) or !
!  glacial-ice (ice). compute (prepare) the matrix coefficients for the !
!  tri-diagonal matrix of the implicit time scheme.                     !
!  (note:  this subroutine only called for sea-ice or glacial ice, but  !
!  not for non-glacial land (ice = 0).                                  !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     stc      - real, soil temperature                          nsoil  !
!     zsoil    - real, soil depth (negative sign, as below grd)  nsoil  !
!     yy       - real, soil temperature at the top of column       1    !
!     zz1      - real,                                             1    !
!     df1      - real, thermal diffusivity and conductivity        1    !
!     ice      - integer, sea-ice flag (=1: sea-ice, =0: land)     1    !
!                                                                       !
!  input/outputs:                                                       !
!     tbot     - real, bottom soil temperature                     1    !
!                                                                       !
!  outputs:                                                             !
!     rhsts    - real, time tendency of soil thermal diffusion   nsoil  !
!     ai       - real, matrix coefficients                       nsold  !
!     bi       - real, matrix coefficients                       nsold  !
!     ci       - real, matrix coefficients                       nsold  !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
      integer, intent(in) :: nsoil, ice

      real (kind=kind_phys), intent(in) :: stc(nsoil), zsoil(nsoil),    &
     &       yy, zz1, df1

!  ---  input/outputs:
      real (kind=kind_phys), intent(inout) :: tbot

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: rhsts(nsoil), ai(nsold),    &
     &       bi(nsold), ci(nsold)

!  ---  locals:
      real (kind=kind_phys) :: ddz, ddz2, denom, dtsdz, dtsdz2,         &
     &       hcpct, ssoil, zbot

      integer :: k

!
!===> ...  begin here
!
!  --- ...  set a nominal universal value of the sea-ice specific heat capacity,
!           hcpct = 1880.0*917.0 = 1.72396e+6 (source:  fei chen, 1995)
!           set bottom of sea-ice pack temperature: tbot = 271.16
!           set a nominal universal value of glacial-ice specific heat capacity,
!           hcpct = 2100.0*900.0 = 1.89000e+6 (source:  bob grumbine, 2005)
!           tbot passed in as argument, value from global data set

      if (ice == 1) then
!  --- ...  sea-ice
        hcpct = 1.72396e+6
        tbot = 271.16
      else
!  --- ...  glacial-ice
        hcpct = 1.89000e+6
      endif

!  --- ...  the input argument df1 is a universally constant value of sea-ice
!           and glacial-ice thermal diffusivity, set in sflx as df1 = 2.2.

!  --- ...  set ice pack depth.  use tbot as ice pack lower boundary temperature
!           (that of unfrozen sea water at bottom of sea ice pack).  assume ice
!           pack is of n=nsoil layers spanning a uniform constant ice pack
!           thickness as defined by zsoil(nsoil) in routine sflx.
!           if glacial-ice, set zbot = -25 meters

      if (ice == 1) then
!  --- ...  sea-ice
        zbot = zsoil(nsoil)
      else
!  --- ...  glacial-ice
        zbot = -25.0
      endif

!  --- ...  calc the matrix coefficients ai, bi, and ci for the top layer

      ddz = 1.0 / (-0.5*zsoil(2))
      ai(1) = 0.0
      ci(1) = (df1*ddz) / (zsoil(1)*hcpct)
      bi(1) = -ci(1) + df1 / (0.5*zsoil(1)*zsoil(1)*hcpct*zz1)

!  --- ...  calc the vertical soil temp gradient btwn the top and 2nd soil
!           layers. recalc/adjust the soil heat flux.  use the gradient and
!           flux to calc rhsts for the top soil layer.

      dtsdz = (stc(1) - stc(2)) / (-0.5*zsoil(2))
      ssoil = df1 * (stc(1) - yy) / (0.5*zsoil(1)*zz1)
      rhsts(1) = (df1*dtsdz - ssoil) / (zsoil(1)*hcpct)

!  --- ...  initialize ddz2

      ddz2 = 0.0

!  --- ...  loop thru the remaining soil layers, repeating the above process

      do k = 2, nsoil

        if (k /= nsoil) then

!  --- ...  calc the vertical soil temp gradient thru this layer.

          denom = 0.5 * (zsoil(k-1) - zsoil(k+1))
          dtsdz2 = (stc(k) - stc(k+1)) / denom

!  --- ...  calc the matrix coef, ci, after calc'ng its partial product.

          ddz2 = 2.0 / (zsoil(k-1) - zsoil(k+1))
          ci(k) = -df1*ddz2 / ((zsoil(k-1) - zsoil(k))*hcpct)

        else

!  --- ...  calc the vertical soil temp gradient thru the lowest layer.

          dtsdz2 = (stc(k) - tbot)                                      &
     &           / (0.5*(zsoil(k-1) + zsoil(k)) - zbot)

!  --- ...  set matrix coef, ci to zero.

          ci(k) = 0.0

        endif   ! end if_k_block

!  --- ...  calc rhsts for this layer after calc'ng a partial product.

        denom = (zsoil(k) - zsoil(k-1)) * hcpct
        rhsts(k) = (df1*dtsdz2 - df1*dtsdz) / denom

!  --- ...  calc matrix coefs, ai, and bi for this layer.

        ai(k) = - df1*ddz / ((zsoil(k-1) - zsoil(k)) * hcpct)
        bi(k) = -(ai(k) + ci(k))

!  --- ...  reset values of dtsdz and ddz for loop to next soil lyr.

        dtsdz = dtsdz2
        ddz   = ddz2

      enddo   ! end do_k_loop
!
      return
!...................................
      end subroutine hrtice
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates/updates the soil temperature field.
      subroutine hstep                                                  &
!...................................
!  ---  inputs:
     &     ( nsoil, stcin, dt,                                          &
!  ---  input/outputs:
     &       rhsts, ai, bi, ci,                                         &
!  ---  outputs:
     &       stcout                                                     &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine hstep calculates/updates the soil temperature field.      !
!                                                                       !
!  subprogram called:  rosr12                                           !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     stcin    - real, soil temperature                          nsoil  !
!     dt       - real, time step                                   1    !
!                                                                       !
!  input/outputs:                                                       !
!     rhsts    - real, time tendency of soil thermal diffusion   nsoil  !
!     ai       - real, matrix coefficients                       nsold  !
!     bi       - real, matrix coefficients                       nsold  !
!     ci       - real, matrix coefficients                       nsold  !
!                                                                       !
!  outputs:                                                             !
!     stcout   - real, updated soil temperature                  nsoil  !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
      integer, intent(in) :: nsoil

      real (kind=kind_phys),  intent(in) :: stcin(nsoil), dt

!  ---  input/outputs:
      real (kind=kind_phys),  intent(inout) :: rhsts(nsoil),            &
     &      ai(nsold), bi(nsold), ci(nsold)

!  ---  outputs:
      real (kind=kind_phys),  intent(out) :: stcout(nsoil)

!  ---  locals:
      integer :: k

      real (kind=kind_phys) :: ciin(nsold), rhstsin(nsoil)

!
!===> ...  begin here
!
!  --- ...  create finite difference values for use in rosr12 routine

      do k = 1, nsoil
        rhsts(k) = rhsts(k) * dt
        ai(k) = ai(k) * dt
        bi(k) = 1.0 + bi(k)*dt
        ci(k) = ci(k) * dt
      enddo

!  --- ...  copy values for input variables before call to rosr12

      do k = 1, nsoil
         rhstsin(k) = rhsts(k)
      enddo

      do k = 1, nsold
        ciin(k) = ci(k)
      enddo

!  --- ...  solve the tri-diagonal matrix equation

      call rosr12                                                       &
!  ---  inputs:
     &     ( nsoil, ai, bi, rhstsin,                                    &
!  ---  input/outputs:
     &       ciin,                                                      &
!  ---  outputs:
     &       ci, rhsts                                                  &
     &     )

!  --- ...  calc/update the soil temps using matrix solution

      do k = 1, nsoil
        stcout(k) = stcin(k) + ci(k)
      enddo
!
      return
!...................................
      end subroutine hstep
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine inverts (solve) the tri-diagonal matrix problem.
      subroutine rosr12                                                 &
!...................................
!  ---  inputs:
     &     ( nsoil, a, b, d,                                            &
!  ---  input/outputs:
     &       c,                                                         &
!  ---  outputs:
     &       p, delta                                                   &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine rosr12 inverts (solve) the tri-diagonal matrix problem    !
!  shown below:                                                         !
!                                                                       !
! ###                                            ### ###  ###   ###  ###!
! #b(1), c(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #!
! #a(2), b(2), c(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #!
! # 0  , a(3), b(3), c(3),  0  ,   . . .  ,    0   # #      #   # d(3) #!
! # 0  ,  0  , a(4), b(4), c(4),   . . .  ,    0   # # p(4) #   # d(4) #!
! # 0  ,  0  ,  0  , a(5), b(5),   . . .  ,    0   # # p(5) #   # d(5) #!
! # .                                          .   # #  .   # = #   .  #!
! # .                                          .   # #  .   #   #   .  #!
! # .                                          .   # #  .   #   #   .  #!
! # 0  , . . . , 0 , a(m-2), b(m-2), c(m-2),   0   # #p(m-2)#   #d(m-2)#!
! # 0  , . . . , 0 ,   0   , a(m-1), b(m-1), c(m-1)# #p(m-1)#   #d(m-1)#!
! # 0  , . . . , 0 ,   0   ,   0   ,  a(m) ,  b(m) # # p(m) #   # d(m) #!
! ###                                            ### ###  ###   ###  ###!
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     a        - real, matrix coefficients                       nsoil  !
!     b        - real, matrix coefficients                       nsoil  !
!     d        - real, soil water time tendency                  nsoil  !
!                                                                       !
!  input/outputs:                                                       !
!     c        - real, matrix coefficients                       nsoil  !
!                                                                       !
!  outputs:                                                             !
!     p        - real,                                           nsoil  !
!     delta    - real,                                           nsoil  !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
      integer, intent(in) :: nsoil

      real (kind=kind_phys),  dimension(nsoil), intent(in) :: a, b, d

!  ---  input/outputs:
      real (kind=kind_phys),  dimension(nsoil), intent(inout) :: c

!  ---  outputs:
      real (kind=kind_phys),  dimension(nsoil), intent(out) :: p, delta

!  ---  locals:
      integer :: k, kk

!
!===> ...  begin here
!
!  --- ...  initialize eqn coef c for the lowest soil layer

      c(nsoil) = 0.0

!  --- ...  solve the coefs for the 1st soil layer

      p(1) = -c(1) / b(1)
      delta(1) = d(1) / b(1)

!  --- ...  solve the coefs for soil layers 2 thru nsoil

      do k = 2, nsoil
        p(k) = -c(k) * ( 1.0 / (b(k) + a (k)*p(k-1)) )
        delta(k) = (d(k) - a(k)*delta(k-1))                              &
     &           * ( 1.0 / (b(k) + a(k)*p(k-1)) )
      enddo

!  --- ...  set p to delta for lowest soil layer

      p(nsoil) = delta(nsoil)

!  --- ...  adjust p for soil layers 2 thru nsoil

      do k = 2, nsoil
         kk = nsoil - k + 1
         p(kk) = p(kk)*p(kk+1) + delta(kk)
      enddo
!
      return
!...................................
      end subroutine rosr12
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates sink/source term of the termal diffusion equation.
      subroutine snksrc                                                 &
!...................................
!  ---  inputs:
     &     ( nsoil, k, tavg, smc, smcmax, psisat, bexp, dt,             &
     &       qtot, zsoil,                                               &
!  ---  input/outputs:
     &       sh2o,                                                      &
!  ---  outputs:
     &       tsrc                                                       &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  subroutine snksrc calculates sink/source term of the termal          !
!  diffusion equation.                                                  !
!                                                                       !
!  subprograms called: frh2o                                            !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     k        - integer, index of soil layers                     1    !
!     tavg     - real, soil layer average temperature              1    !
!     smc      - real, total soil moisture                         1    !
!     smcmax   - real, porosity                                    1    !
!     psisat   - real, saturated soil potential                    1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     dt       - real, time step                                   1    !
!     qtot     - real, tot vertical diff of heat flux              1    !
!     zsoil    - real, soil layer depth below ground (negative)  nsoil  !
!                                                                       !
!  input/outputs:                                                       !
!     sh2o     - real, available liqued water                      1    !
!                                                                       !
!  outputs:                                                             !
!     tsrc     - real, heat source/sink                            1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  parameter constants:
      real (kind=kind_phys), parameter :: dh2o = 1.0000e3

!  ---  inputs:
      integer, intent(in) :: nsoil, k

      real (kind=kind_phys), intent(in) :: tavg, smc, smcmax, psisat,   &
     &       bexp, dt, qtot, zsoil(nsoil)

!  ---  input/outputs:
      real (kind=kind_phys), intent(inout) :: sh2o

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: tsrc

!  ---  locals:
      real (kind=kind_phys) :: dz, free, xh2o

!  ---  external functions:
!     real (kind=kind_phys) :: frh2o

!urban
       if (ivegsrc == 1)then
            if ( vegtyp == 13 ) df1=3.24
       endif
!
!===> ...  begin here
!
      if (k == 1) then
        dz = -zsoil(1)
      else
        dz = zsoil(k-1) - zsoil(k)
      endif

!  --- ...  via function frh2o, compute potential or 'equilibrium' unfrozen
!           supercooled free water for given soil type and soil layer temperature.
!           function frh20 invokes eqn (17) from v. koren et al (1999, jgr, vol.
!           104, pg 19573).  (aside:  latter eqn in journal in centigrade units.
!           routine frh2o use form of eqn in kelvin units.)

!     free = frh2o( tavg,smc,sh2o,smcmax,bexp,psisat )

      call frh2o                                                        &
!  ---  inputs:
     &     ( tavg, smc, sh2o, smcmax, bexp, psisat,                     &
!  ---  outputs:
     &       free                                                       &
     &     )


!  --- ...  in next block of code, invoke eqn 18 of v. koren et al (1999, jgr,
!           vol. 104, pg 19573.)  that is, first estimate the new amountof liquid
!           water, 'xh2o', implied by the sum of (1) the liquid water at the begin
!           of current time step, and (2) the freeze of thaw change in liquid
!           water implied by the heat flux 'qtot' passed in from routine hrt.
!           second, determine if xh2o needs to be bounded by 'free' (equil amt) or
!           if 'free' needs to be bounded by xh2o.

      xh2o = sh2o + qtot*dt / (dh2o*lsubf*dz)

!  --- ...  first, if freezing and remaining liquid less than lower bound, then
!           reduce extent of freezing, thereby letting some or all of heat flux
!           qtot cool the soil temp later in routine hrt.

      if ( xh2o < sh2o .and. xh2o < free) then
        if ( free > sh2o ) then
          xh2o = sh2o
        else
          xh2o = free
        endif
      endif

!  --- ...  second, if thawing and the increase in liquid water greater than
!           upper bound, then reduce extent of thaw, thereby letting some or
!           all of heat flux qtot warm the soil temp later in routine hrt.

      if ( xh2o > sh2o .and. xh2o > free )  then
        if ( free < sh2o ) then
          xh2o = sh2o
        else
          xh2o = free
        endif
      endif

      xh2o = max( min( xh2o, smc ), 0.0 )

!  --- ...  calculate phase-change heat source/sink term for use in routine hrt
!           and update liquid water to reflcet final freeze/thaw increment.

      tsrc = -dh2o * lsubf * dz * (xh2o - sh2o) / dt
      sh2o = xh2o
!
      return
!...................................
      end subroutine snksrc
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates the right hand side of the time tendency
!! term of the soil water diffusion equation. Also to compute
!! (prepare) the matrix coefficients for the tri-diagonal matrix of 
!! the implicit time scheme.
      subroutine srt                                                    &
!...................................
!  ---  inputs:
     &     ( nsoil, edir, et, sh2o, sh2oa, pcpdrp, zsoil, dwsat,        &
     &       dksat, smcmax, bexp, dt, smcwlt, slope, kdt, frzx, sice,   &
!  ---  outputs:
     &       rhstt, runoff1, runoff2, ai, bi, ci                        &
     &     )

! ===================================================================== !
!  description:                                                         !
!    subroutine srt calculates the right hand side of the time tendency !
!    term of the soil water diffusion equation.  also to compute        !
!    ( prepare ) the matrix coefficients for the tri-diagonal matrix    !
!    of the implicit time scheme.                                       !
!                                                                       !
!  subprogram called:  wdfcnd                                           !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     edir     - real, direct soil evaporation                     1    !
!     et       - real, plant transpiration                       nsoil  !
!     sh2o     - real, unfrozen soil moisture                    nsoil  !
!     sh2oa    - real,                                           nsoil  !
!     pcpdrp   - real, combined prcp and drip                      1    !
!     zsoil    - real, soil layer depth below ground             nsoil  !
!     dwsat    - real, saturated soil diffusivity                  1    !
!     dksat    - real, saturated soil hydraulic conductivity       1    !
!     smcmax   - real, porosity                                    1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     dt       - real, time step                                   1    !
!     smcwlt   - real, wilting point                               1    !
!     slope    - real, linear reservoir coefficient                1    !
!     kdt      - real,                                             1    !
!     frzx     - real, frozen ground parameter                     1    !
!     sice     - real, ice content at each soil layer            nsoil  !
!                                                                       !
!  outputs:                                                             !
!     rhstt    - real, soil water time tendency                  nsoil  !
!     runoff1  - real, surface runoff not infiltrating sfc         1    !
!     runoff2  - real, sub surface runoff (baseflow)               1    !
!     ai       - real, matrix coefficients                       nsold  !
!     bi       - real, matrix coefficients                       nsold  !
!     ci       - real, matrix coefficients                       nsold  !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  inputs:
      integer, intent(in) :: nsoil

      real (kind=kind_phys), dimension(nsoil), intent(in) :: et,        &
     &       sh2o, sh2oa, zsoil, sice

      real (kind=kind_phys), intent(in) :: edir, pcpdrp, dwsat, dksat,  &
     &       smcmax, smcwlt, bexp, dt, slope, kdt, frzx

!  --- outputs:
      real (kind=kind_phys), intent(out) :: runoff1, runoff2,           &
     &      rhstt(nsoil), ai(nsold), bi(nsold), ci(nsold)


!  ---  locals:
      real (kind=kind_phys) :: acrt, dd, ddt, ddz, ddz2, denom, denom2, &
     &       dice, dsmdz, dsmdz2, dt1, fcr, infmax, mxsmc, mxsmc2, px,  &
     &       numer, pddum, sicemax, slopx, smcav, sstt, sum, val, wcnd, &
     &       wcnd2, wdf, wdf2, dmax(nsold)

      integer :: cvfrz, ialp1, iohinf, j, jj, k, ks
!
!===> ...  begin here
!
!  --- ...  frozen ground version:
!           reference frozen ground parameter, cvfrz, is a shape parameter
!           of areal distribution function of soil ice content which equals
!           1/cv. cv is a coefficient of spatial variation of soil ice content.
!           based on field data cv depends on areal mean of frozen depth, and
!           it close to constant = 0.6 if areal mean frozen depth is above 20 cm.
!           that is why parameter cvfrz = 3 (int{1/0.6*0.6}). current logic
!           doesn't allow cvfrz be bigger than 3

        parameter (cvfrz = 3)

c ----------------------------------------------------------------------
!  --- ...  determine rainfall infiltration rate and runoff.  include
!           the infiltration formule from schaake and koren model.
!           modified by q duan

      iohinf = 1

!  --- ... let sicemax be the greatest, if any, frozen water content within
!          soil layers.

      sicemax = 0.0
      do ks = 1, nsoil
       if (sice(ks) > sicemax) sicemax = sice(ks)
      enddo

!  --- ...  determine rainfall infiltration rate and runoff

      pddum = pcpdrp
      runoff1 = 0.0

      if (pcpdrp /= 0.0) then

!  --- ...  modified by q. duan, 5/16/94

        dt1 = dt/86400.
        smcav = smcmax - smcwlt
        dmax(1) = -zsoil(1) * smcav

!  --- ...  frozen ground version:

        dice = -zsoil(1) * sice(1)

        dmax(1) = dmax(1)*(1.0 - (sh2oa(1)+sice(1)-smcwlt)/smcav)
        dd = dmax(1)

        do ks = 2, nsoil

!  --- ...  frozen ground version:

          dice = dice + ( zsoil(ks-1) - zsoil(ks) ) * sice(ks)

          dmax(ks) = (zsoil(ks-1)-zsoil(ks))*smcav
          dmax(ks) = dmax(ks)*(1.0 - (sh2oa(ks)+sice(ks)-smcwlt)/smcav)
          dd = dd + dmax(ks)
        enddo

!  --- ...  val = (1.-exp(-kdt*sqrt(dt1)))
!           in below, remove the sqrt in above

        val = 1.0 - exp(-kdt*dt1)
        ddt = dd * val

        px = pcpdrp * dt
        if (px < 0.0) px = 0.0

        infmax = (px*(ddt/(px+ddt)))/dt

!  --- ...  frozen ground version:
!           reduction of infiltration based on frozen ground parameters

        fcr = 1.

        if (dice > 1.e-2) then
          acrt = cvfrz * frzx / dice
          sum = 1.

          ialp1 = cvfrz - 1
          do j = 1, ialp1
            k = 1

            do jj = j+1,ialp1
              k = k * jj
            enddo

            sum = sum + (acrt**( cvfrz-j)) / float (k)
          enddo

          fcr = 1.0 - exp(-acrt) * sum
        endif

        infmax = infmax * fcr

!  --- ...  correction of infiltration limitation:
!           if infmax .le. hydrolic conductivity assign infmax the value
!           of hydrolic conductivity

!       mxsmc = max ( sh2oa(1), sh2oa(2) )
        mxsmc = sh2oa(1)

        call wdfcnd                                                     &
!  ---  inputs:
     &     ( mxsmc, smcmax, bexp, dksat, dwsat, sicemax,                &
!  ---  outputs:
     &       wdf, wcnd                                                  &
     &     )

        infmax = max( infmax, wcnd )
        infmax = min( infmax, px )

        if (pcpdrp > infmax) then
          runoff1 = pcpdrp - infmax
          pddum   = infmax
        endif

      endif   ! end if_pcpdrp_block

!  --- ... to avoid spurious drainage behavior, 'upstream differencing'
!          in line below replaced with new approach in 2nd line:
!          'mxsmc = max(sh2oa(1), sh2oa(2))'

      mxsmc = sh2oa(1)

      call wdfcnd                                                       &
!  ---  inputs:
     &     ( mxsmc, smcmax, bexp, dksat, dwsat, sicemax,                &
!  ---  outputs:
     &       wdf, wcnd                                                  &
     &     )

!  --- ...  calc the matrix coefficients ai, bi, and ci for the top layer

      ddz = 1.0 / ( -.5*zsoil(2) )
      ai(1) = 0.0
      bi(1) = wdf * ddz / ( -zsoil(1) )
      ci(1) = -bi(1)

!  --- ...  calc rhstt for the top layer after calc'ng the vertical soil
!           moisture gradient btwn the top and next to top layers.

      dsmdz = ( sh2o(1) - sh2o(2) ) / ( -.5*zsoil(2) )
      rhstt(1) = (wdf*dsmdz + wcnd - pddum + edir + et(1)) / zsoil(1)
      sstt = wdf * dsmdz + wcnd + edir + et(1)

!  --- ...  initialize ddz2

      ddz2 = 0.0

!  --- ...  loop thru the remaining soil layers, repeating the abv process

      do k = 2, nsoil
        denom2 = (zsoil(k-1) - zsoil(k))

        if (k /= nsoil) then
          slopx = 1.0

!  --- ...  again, to avoid spurious drainage behavior, 'upstream differencing'
!           in line below replaced with new approach in 2nd line:
!           'mxsmc2 = max (sh2oa(k), sh2oa(k+1))'

          mxsmc2 = sh2oa(k)

          call wdfcnd                                                   &
!  ---  inputs:
     &     ( mxsmc2, smcmax, bexp, dksat, dwsat, sicemax,               &
!  ---  outputs:
     &       wdf2, wcnd2                                                &
     &     )

!  --- ...  calc some partial products for later use in calc'ng rhstt

          denom = (zsoil(k-1) - zsoil(k+1))
          dsmdz2 = (sh2o(k) - sh2o(k+1)) / (denom * 0.5)

!  --- ...  calc the matrix coef, ci, after calc'ng its partial product

          ddz2 = 2.0 / denom
          ci(k) = -wdf2 * ddz2 / denom2

        else   ! if_k_block

!  --- ...  slope of bottom layer is introduced

          slopx = slope

!  --- ...  retrieve the soil water diffusivity and hydraulic conductivity
!           for this layer

          call wdfcnd                                                   &
!  ---  inputs:
     &     ( sh2oa(nsoil), smcmax, bexp, dksat, dwsat, sicemax,         &
!  ---  outputs:
     &       wdf2, wcnd2                                                &
     &     )

!  --- ...  calc a partial product for later use in calc'ng rhstt
          dsmdz2 = 0.0

!  --- ...  set matrix coef ci to zero

          ci(k) = 0.0

        endif   ! end if_k_block

!  --- ...  calc rhstt for this layer after calc'ng its numerator

        numer = wdf2*dsmdz2 + slopx*wcnd2 - wdf*dsmdz - wcnd + et(k)
        rhstt(k) = numer / (-denom2)

!  --- ...  calc matrix coefs, ai, and bi for this layer

        ai(k) = -wdf * ddz / denom2
        bi(k) = -( ai(k) + ci(k) )

!  --- ...  reset values of wdf, wcnd, dsmdz, and ddz for loop to next lyr
!      runoff2:  sub-surface or baseflow runoff

        if (k == nsoil) then
          runoff2 = slopx * wcnd2
        endif

        if (k /= nsoil) then
          wdf  = wdf2
          wcnd = wcnd2
          dsmdz= dsmdz2
          ddz  = ddz2
        endif
      enddo   ! end do_k_loop
!
      return
!...................................
      end subroutine srt
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates/updates soil moisture content values and
!! canopy moisture content values.
      subroutine sstep                                                  &
!...................................
!  ---  inputs:
     &     ( nsoil, sh2oin, rhsct, dt, smcmax, cmcmax, zsoil, sice,     &
!  ---  input/outputs:
     &       cmc, rhstt, ai, bi, ci,                                    &
!  ---  outputs:
     &       sh2oout, runoff3, smc                                      &
     &     )

! ===================================================================== !
!  description:                                                         !
!    subroutine sstep calculates/updates soil moisture content values   !
!    and canopy moisture content values.                                !
!                                                                       !
!  subprogram called:  rosr12                                           !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     sh2oin   - real, unfrozen soil moisture                    nsoil  !
!     rhsct    - real,                                             1    !
!     dt       - real, time step                                   1    !
!     smcmax   - real, porosity                                    1    !
!     cmcmax   - real, maximum canopy water parameters             1    !
!     zsoil    - real, soil layer depth below ground             nsoil  !
!     sice     - real, ice content at each soil layer            nsoil  !
!                                                                       !
!  input/outputs:                                                       !
!     cmc      - real, canopy moisture content                     1    !
!     rhstt    - real, soil water time tendency                  nsoil  !
!     ai       - real, matrix coefficients                       nsold  !
!     bi       - real, matrix coefficients                       nsold  !
!     ci       - real, matrix coefficients                       nsold  !
!                                                                       !
!  outputs:                                                             !
!     sh2oout  - real, updated soil moisture content             nsoil  !
!     runoff3  - real, excess of porosity                          1    !
!     smc      - real, total soil moisture                       nsoil  !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  input:
      integer, intent(in) :: nsoil

      real (kind=kind_phys), dimension(nsoil), intent(in) :: sh2oin,    &
     &       zsoil, sice

      real (kind=kind_phys), intent(in) :: rhsct, dt, smcmax, cmcmax

!  ---  inout/outputs:
      real (kind=kind_phys), intent(inout) :: cmc, rhstt(nsoil),        &
     &       ai(nsold), bi(nsold), ci(nsold)

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: sh2oout(nsoil), runoff3,    &
     &       smc(nsoil)

!  ---  locals:
      real (kind=kind_phys) :: ciin(nsold), rhsttin(nsoil), ddz, stot,  &
     &       wplus

      integer :: i, k, kk11
!
!===> ...  begin here
!
!  --- ...  create 'amount' values of variables to be input to the
!           tri-diagonal matrix routine.

      do k = 1, nsoil
        rhstt(k) = rhstt(k) * dt
        ai(k) = ai(k) * dt
        bi(k) = 1. + bi(k) * dt
        ci(k) = ci(k) * dt
      enddo

!  --- ...  copy values for input variables before call to rosr12

      do k = 1, nsoil
        rhsttin(k) = rhstt(k)
      enddo

      do k = 1, nsold
        ciin(k) = ci(k)
      enddo

!  --- ...  call rosr12 to solve the tri-diagonal matrix

      call rosr12                                                       &
!  ---  inputs:
     &     ( nsoil, ai, bi, rhsttin,                                    &
!  ---  input/outputs:
     &       ciin,                                                      &
!  ---  outputs:
     &       ci, rhstt                                                  &
     &     )

!  --- ...  sum the previous smc value and the matrix solution to get
!           a new value.  min allowable value of smc will be 0.02.
!      runoff3: runoff within soil layers

      wplus   = 0.0
      runoff3 = 0.0
      ddz     = -zsoil(1)

      do k = 1, nsoil
        if (k /= 1) ddz = zsoil(k - 1) - zsoil(k)

        sh2oout(k) = sh2oin(k) + ci(k) + wplus/ddz

        stot = sh2oout(k) + sice(k)
        if (stot > smcmax) then
          if (k == 1) then
            ddz = -zsoil(1)
          else
            kk11 = k - 1
            ddz = -zsoil(k) + zsoil(kk11)
          endif

          wplus = (stot - smcmax) * ddz
        else
          wplus = 0.0
        endif

        smc(k) = max( min( stot, smcmax ), 0.02 )
        sh2oout(k) = max( smc(k)-sice(k), 0.0 )
      enddo

      runoff3 = wplus

!  --- ...  update canopy water content/interception (cmc).  convert rhsct to
!           an 'amount' value and add to previous cmc value to get new cmc.

      cmc = cmc + dt*rhsct
      if (cmc < 1.e-20) cmc = 0.0
      cmc = min( cmc, cmcmax )
!
      return
!...................................
      end subroutine sstep
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates temperature on the boundary of the 
!! layer by interpolation of the middle layer temperatures.
      subroutine tbnd                                                   &
!...................................
!  ---  inputs:
     &     ( tu, tb, zsoil, zbot, k, nsoil,                             &
!  ---  outputs:
     &       tbnd1                                                      &
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!    subroutine tbnd calculates temperature on the boundary of the      !
!    layer by interpolation of the middle layer temperatures            !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     tu       - real, soil temperature                            1    !
!     tb       - real, bottom soil temp                            1    !
!     zsoil    - real, soil layer depth                          nsoil  !
!     zbot     - real, specify depth of lower bd soil              1    !
!     k        - integer, soil layer index                         1    !
!     nsoil    - integer, number of soil layers                    1    !
!                                                                       !
!  outputs:                                                             !
!     tbnd1    - real, temperature at bottom interface of the lyr  1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  input:
      integer, intent(in) :: k, nsoil

      real (kind=kind_phys), intent(in) :: tu, tb, zbot, zsoil(nsoil)

!  ---  output:
      real (kind=kind_phys), intent(out) :: tbnd1

!  ---  locals:
      real (kind=kind_phys) :: zb, zup

!  --- ...  use surface temperature on the top of the first layer

      if (k == 1) then
        zup = 0.0
      else
        zup = zsoil(k-1)
      endif

!  --- ...  use depth of the constant bottom temperature when interpolate
!           temperature into the last layer boundary

      if (k == nsoil) then
        zb = 2.0*zbot - zsoil(k)
      else
        zb = zsoil(k+1)
      endif

!  --- ...  linear interpolation between the average layer temperatures

      tbnd1 = tu + (tb-tu)*(zup-zsoil(k))/(zup-zb)
!
      return
!...................................
      end subroutine tbnd
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates soil layer average temperature (tavg)
!! in freezing/thawing layer using up, down, and middle layer
!! temperature (tup, tdn, tm), where tup is at top boundary of layer,
!! tdn is at bottom boundary of layer. tm is layer prognostic state
!! temperature.
      subroutine tmpavg                                                 &
!...................................
!  ---  inputs:
     &     ( tup, tm, tdn, zsoil, nsoil, k,                             &
!  ---  outputs:
     &       tavg                                                       &
     &     )

! ===================================================================== !
!  description:                                                         !
!    subroutine tmpavg calculates soil layer average temperature (tavg) !
!    in freezing/thawing layer using up, down, and middle layer         !
!    temperatures (tup, tdn, tm), where tup is at top boundary of       !
!    layer, tdn is at bottom boundary of layer.  tm is layer prognostic !
!    state temperature.                                                 !
!                                                                       !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     tup      - real, temperature ar top boundary of layer        1    !
!     tm       - real, layer prognostic state temperature          1    !
!     tdn      - real, temperature ar bottom boundary of layer     1    !
!     zsoil    - real, soil layer depth                          nsoil  !
!     nsoil    - integer, number of soil layers                    1    !
!     k        - integer, layer index                              1    !
!  outputs:                                                             !
!     tavg     - real, soil layer average temperature              1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  input:
      integer, intent(in) :: nsoil, k

      real (kind=kind_phys), intent(in) :: tup, tm, tdn, zsoil(nsoil)

!  ---  output:
      real (kind=kind_phys), intent(out) :: tavg

!  ---  locals:
      real (kind=kind_phys) :: dz, dzh, x0, xdn, xup
!
!===> ...  begin here
!
      if (k == 1) then
        dz = -zsoil(1)
      else
        dz = zsoil(k-1) - zsoil(k)
      endif

      dzh = dz * 0.5

      if (tup < tfreez) then

        if (tm < tfreez) then
          if (tdn < tfreez) then               ! tup, tm, tdn < t0
            tavg = (tup + 2.0*tm + tdn) / 4.0
          else                                 ! tup & tm < t0,  tdn >= t0
            x0 = (tfreez - tm) * dzh / (tdn - tm)
            tavg = 0.5*(tup*dzh + tm*(dzh+x0)+tfreez*(2.*dzh-x0)) / dz
          endif
        else
          if (tdn < tfreez) then               ! tup < t0, tm >= t0, tdn < t0
            xup  = (tfreez-tup) * dzh / (tm-tup)
            xdn  = dzh - (tfreez-tm) * dzh / (tdn-tm)
            tavg = 0.5*(tup*xup + tfreez*(2.*dz-xup-xdn)+tdn*xdn) / dz
          else                                 ! tup < t0, tm >= t0, tdn >= t0
            xup  = (tfreez-tup) * dzh / (tm-tup)
            tavg = 0.5*(tup*xup + tfreez*(2.*dz-xup)) / dz
          endif
        endif

      else    ! if_tup_block

        if (tm < tfreez) then
          if (tdn < tfreez) then               ! tup >= t0, tm < t0, tdn < t0
            xup  = dzh - (tfreez-tup) * dzh / (tm-tup)
            tavg = 0.5*(tfreez*(dz-xup) + tm*(dzh+xup)+tdn*dzh) / dz
          else                                 ! tup >= t0, tm < t0, tdn >= t0
            xup  = dzh - (tfreez-tup) * dzh / (tm-tup)
            xdn  = (tfreez-tm) * dzh / (tdn-tm)
            tavg = 0.5 * (tfreez*(2.*dz-xup-xdn) + tm*(xup+xdn)) / dz
          endif
        else
          if (tdn < tfreez) then               ! tup >= t0, tm >= t0, tdn < t0
            xdn  = dzh - (tfreez-tm) * dzh / (tdn-tm)
            tavg = (tfreez*(dz-xdn) + 0.5*(tfreez+tdn)*xdn) / dz
          else                                 ! tup >= t0, tm >= t0, tdn >= t0
            tavg = (tup + 2.0*tm + tdn) / 4.0
          endif
        endif

      endif   ! end if_tup_block
!
      return
!...................................
      end subroutine tmpavg
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates transpiration for the veg class.
      subroutine transp                                                    &
!...................................
!  ---  inputs:
     &     ( nsoil, nroot, etp1, smc, smcwlt, smcref,                   &
     &       cmc, cmcmax, zsoil, shdfac, pc, cfactr, rtdis,             &
!  ---  outputs:
     &       et1                                                        &
     &     )

! ===================================================================== !
!  description:                                                         !
!     subroutine transp calculates transpiration for the veg class.     !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     nsoil    - integer, number of soil layers                    1    !
!     nroot    - integer, number of root layers                    1    !
!     etp1     - real, potential evaporation                       1    !
!     smc      - real, unfrozen soil moisture                    nsoil  !
!     smcwlt   - real, wilting point                               1    !
!     smcref   - real, soil mois threshold                         1    !
!     cmc      - real, canopy moisture content                     1    !
!     cmcmax   - real, maximum canopy water parameters             1    !
!     zsoil    - real, soil layer depth below ground             nsoil  !
!     shdfac   - real, aeral coverage of green vegetation          1    !
!     pc       - real, plant coeff                                 1    !
!     cfactr   - real, canopy water parameters                     1    !
!     rtdis    - real, root distribution                         nsoil  !
!                                                                       !
!  outputs:                                                             !
!     et1      - real, plant transpiration                       nsoil  !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  input:
      integer, intent(in) :: nsoil, nroot

      real (kind=kind_phys), intent(in) :: etp1, smcwlt, smcref,        &
     &       cmc, cmcmax, shdfac, pc, cfactr

      real (kind=kind_phys), dimension(nsoil), intent(in) :: smc,       &
     &       zsoil, rtdis

!  ---  output:
      real (kind=kind_phys), dimension(nsoil), intent(out) :: et1

!  ---  locals:
      real (kind=kind_phys) :: denom, etp1a, rtx, sgx, gx(7)

      integer :: i, k
!
!===> ...  begin here
!
!  --- ...  initialize plant transp to zero for all soil layers.

      do k = 1, nsoil
        et1(k) = 0.0
      enddo

!  --- ...  calculate an 'adjusted' potential transpiration
!           if statement below to avoid tangent linear problems near zero
!           note: gx and other terms below redistribute transpiration by layer,
!           et(k), as a function of soil moisture availability, while preserving
!           total etp1a.

      if (cmc /= 0.0) then
        etp1a = shdfac * pc * etp1 * (1.0 - (cmc /cmcmax) ** cfactr)
      else
        etp1a = shdfac * pc * etp1
      endif

      sgx = 0.0
      do i = 1, nroot
        gx(i) = ( smc(i) - smcwlt ) / ( smcref - smcwlt )
        gx(i) = max ( min ( gx(i), 1.0 ), 0.0 )
        sgx = sgx + gx(i)
      enddo
      sgx = sgx / nroot

      denom = 0.0
      do i = 1, nroot
        rtx = rtdis(i) + gx(i) - sgx
        gx(i) = gx(i) * max ( rtx, 0.0 )
        denom = denom + gx(i)
      enddo
      if (denom <= 0.0) denom = 1.0

      do i = 1, nroot
        et1(i) = etp1a * gx(i) / denom
      enddo

!  --- ...  above code assumes a vertically uniform root distribution
!           code below tests a variable root distribution

!     et(1) = ( zsoil(1) / zsoil(nroot) ) * gx * etp1a
!     et(1) = ( zsoil(1) / zsoil(nroot) ) * etp1a

!  --- ...  using root distribution as weighting factor

!     et(1) = rtdis(1) * etp1a
!     et(1) = etp1a * part(1)

!  --- ...  loop down thru the soil layers repeating the operation above,
!           but using the thickness of the soil layer (rather than the
!           absolute depth of each layer) in the final calculation.

!     do k = 2, nroot
!       gx = ( smc(k) - smcwlt ) / ( smcref - smcwlt )
!       gx = max ( min ( gx, 1.0 ), 0.0 )
!  --- ...  test canopy resistance
!       gx = 1.0
!       et(k) = ((zsoil(k)-zsoil(k-1))/zsoil(nroot))*gx*etp1a
!       et(k) = ((zsoil(k)-zsoil(k-1))/zsoil(nroot))*etp1a

!  --- ...  using root distribution as weighting factor

!       et(k) = rtdis(k) * etp1a
!       et(k) = etp1a*part(k)
!     enddo

!
      return
!...................................
      end subroutine transp
!-----------------------------------


!-----------------------------------
!>\ingroup Noah_LSM
!> This subroutine calculates soil water diffusivity and soil
!! hydraulic conductivity.
      subroutine wdfcnd                                                 &
!...................................
!  ---  inputs:
     &     ( smc, smcmax, bexp, dksat, dwsat, sicemax,                  &
!  ---  outputs:
     &       wdf, wcnd                                                  &
     &     )

! ===================================================================== !
!  description:                                                         !
!     subroutine wdfcnd calculates soil water diffusivity and soil      !
!     hydraulic conductivity.                                           !
!                                                                       !
!  subprogram called:  none                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     smc      - real, layer total soil moisture                   1    !
!     smcmax   - real, porosity                                    1    !
!     bexp     - real, soil type "b" parameter                     1    !
!     dksat    - real, saturated soil hydraulic conductivity       1    !
!     dwsat    - real, saturated soil diffusivity                  1    !
!     sicemax  - real, max frozen water content in soil layer      1    !
!                                                                       !
!  outputs:                                                             !
!     wdf      - real, soil water diffusivity                      1    !
!     wcnd     - real, soil hydraulic conductivity                 1    !
!                                                                       !
!  ====================    end of description    =====================  !
!
!  ---  input:
      real (kind=kind_phys), intent(in)  :: smc, smcmax, bexp, dksat,   &
     &       dwsat, sicemax

!  ---  output:
      real (kind=kind_phys), intent(out) :: wdf, wcnd

!  ---  locals:
      real (kind=kind_phys) :: expon, factr1, factr2, vkwgt
!
!===> ...  begin here
!
!  --- ...  calc the ratio of the actual to the max psbl soil h2o content

      factr1 = 0.2 / smcmax
      factr2 = smc / smcmax

!  --- ...  prep an expntl coef and calc the soil water diffusivity

      expon = bexp + 2.0
      wdf = dwsat * factr2 ** expon

!  --- ...  frozen soil hydraulic diffusivity.  very sensitive to the vertical
!           gradient of unfrozen water. the latter gradient can become very
!           extreme in freezing/thawing situations, and given the relatively
!           few and thick soil layers, this gradient sufferes serious
!           trunction errors yielding erroneously high vertical transports of
!           unfrozen water in both directions from huge hydraulic diffusivity.
!           therefore, we found we had to arbitrarily constrain wdf
!
!           version d_10cm:  .......  factr1 = 0.2/smcmax
!           weighted approach.......  pablo grunmann, 28_sep_1999.

      if (sicemax > 0.0) then
        vkwgt = 1.0 / (1.0 + (500.0*sicemax)**3.0)
        wdf = vkwgt*wdf + (1.0- vkwgt)*dwsat*factr1**expon
      endif

!  --- ...  reset the expntl coef and calc the hydraulic conductivity

      expon = (2.0 * bexp) + 3.0
      wcnd = dksat * factr2 ** expon
!
      return
!...................................
      end subroutine wdfcnd
!-----------------------------------

! =========================== !
!     end contain programs    !
! =========================== !

!...................................
      end subroutine gfssflx
!-----------------------------------
!! @}
