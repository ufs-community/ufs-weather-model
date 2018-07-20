!>  \file radlw_main.f
!!  This file contains NCEP's modifications of the rrtmg-lw radiation
!!  code from AER.

!\defgroup RRTMG GFS RRTMG Shortwave/Longwave Radiation 
!
!!!!!  ==============================================================  !!!!!
!!!!!               lw-rrtm3 radiation package description             !!!!!
!!!!!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtm-lw radiation    !
!   code from aer inc.                                                     !
!                                                                          !
!    the lw-rrtm3 package includes these parts:                            !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    the 'radlw_rrtm3_param.f' contains:                                   !
!                                                                          !
!       'module_radlw_parameters'  -- band parameters set up               !
!                                                                          !
!    the 'radlw_rrtm3_datatb.f' contains:                                  !
!                                                                          !
!       'module_radlw_avplank'     -- plank flux data                      !
!       'module_radlw_ref'         -- reference temperature and pressure   !
!       'module_radlw_cldprlw'     -- cloud property coefficients          !
!       'module_radlw_kgbnn'       -- absorption coeffients for 16         !
!                                     bands, where nn = 01-16              !
!                                                                          !
!    the 'radlw_rrtm3_main.f' contains:                                    !
!                                                                          !
!       'rrtmg_lw'        -- main lw radiation transfer                    !
!                                                                          !
!    in the main module 'rrtmg_lw' there are only two                      !
!    externally callable subroutines:                                      !
!                                                                          !
!                                                                          !
!       'lwrad'     -- main lw radiation routine                           !
!          inputs:                                                         !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                         !
!            clouds,icseed,aerosols,sfemis,sfgtmp,                         !
!            npts, nlay, nlp1, lprnt,                                      !
!          outputs:                                                        !
!            hlwc,topflx,sfcflx,                                           !
!!         optional outputs:                                               !
!            HLW0,HLWB,FLXPRF)                                             !
!                                                                          !
!       'rlwinit'   -- initialization routine                              !
!          inputs:                                                         !
!           ( me )                                                         !
!          outputs:                                                        !
!           (none)                                                         !
!                                                                          !
!    all the lw radiation subprograms become contained subprograms         !
!    in module 'rrtmg_lw' and many of them are not directly                !
!    accessable from places outside the module.                            !
!                                                                          !
!    derived data type constructs used:                                    !
!                                                                          !
!     1. radiation flux at toa: (from module 'module_radlw_parameters')    !
!          topflw_type   -  derived data type for toa rad fluxes           !
!            upfxc              total sky upward flux at toa               !
!            upfx0              clear sky upward flux at toa               !
!                                                                          !
!     2. radiation flux at sfc: (from module 'module_radlw_parameters')    !
!          sfcflw_type   -  derived data type for sfc rad fluxes           !
!            upfxc              total sky upward flux at sfc               !
!            upfx0              clear sky upward flux at sfc               !
!            dnfxc              total sky downward flux at sfc             !
!            dnfx0              clear sky downward flux at sfc             !
!                                                                          !
!     3. radiation flux profiles(from module 'module_radlw_parameters')    !
!          proflw_type    -  derived data type for rad vertical prof       !
!            upfxc              level upward flux for total sky            !
!            dnfxc              level downward flux for total sky          !
!            upfx0              level upward flux for clear sky            !
!            dnfx0              level downward flux for clear sky          !
!                                                                          !
!    external modules referenced:                                          !
!                                                                          !
!       'module physparam'                                                 !
!       'module physcons'                                                  !
!       'mersenne_twister'                                                 !
!                                                                          !
!    compilation sequence is:                                              !
!                                                                          !
!       'radlw_rrtm3_param.f'                                              !
!       'radlw_rrtm3_datatb.f'                                             !
!       'radlw_rrtm3_main.f'                                               !
!                                                                          !
!    and all should be put in front of routines that use lw modules        !
!                                                                          !
!==========================================================================!
!                                                                          !
!    the original aer's program declarations:                              !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          |
!  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
!  This software may be used, copied, or redistributed as long as it is    |
!  not sold and this copyright notice is reproduced on each copy made.     |
!  This model is provided as is without any express or implied warranties. |
!                       (http://www.rtweb.aer.com/)                        |
!                                                                          |
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
! ************************************************************************ !
!                                                                          !
!                              rrtmg_lw                                    !
!                                                                          !
!                                                                          !
!                   a rapid radiative transfer model                       !
!                       for the longwave region                            !
!             for application to general circulation models                !
!                                                                          !
!                                                                          !
!            atmospheric and environmental research, inc.                  !
!                        131 hartwell avenue                               !
!                        lexington, ma 02421                               !
!                                                                          !
!                           eli j. mlawer                                  !
!                        jennifer s. delamere                              !
!                         michael j. iacono                                !
!                         shepard a. clough                                !
!                                                                          !
!                                                                          !
!                       email:  miacono@aer.com                            !
!                       email:  emlawer@aer.com                            !
!                       email:  jdelamer@aer.com                           !
!                                                                          !
!        the authors wish to acknowledge the contributions of the          !
!        following people:  steven j. taubman, karen cady-pereira,         !
!        patrick d. brown, ronald e. farren, luke chen, robert bergstrom.  !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    references:                                                           !
!    (rrtm_lw/rrtmg_lw):                                                   !
!      clough, s.A., m.w. shephard, e.j. mlawer, j.s. delamere,            !
!      m.j. iacono, k. cady-pereira, s. boukabara, and p.d. brown:         !
!      atmospheric radiative transfer modeling: a summary of the aer       !
!      codes, j. quant. spectrosc. radiat. transfer, 91, 233-244, 2005.    !
!                                                                          !
!      mlawer, e.j., s.j. taubman, p.d. brown, m.j. iacono, and s.a.       !
!      clough:  radiative transfer for inhomogeneous atmospheres: rrtm,    !
!      a validated correlated-k model for the longwave.  j. geophys. res., !
!      102, 16663-16682, 1997.                                             !
!                                                                          !
!    (mcica):                                                              !
!      pincus, r., h. w. barker, and j.-j. morcrette: a fast, flexible,    !
!      approximation technique for computing radiative transfer in         !
!      inhomogeneous cloud fields, j. geophys. res., 108(d13), 4376,       !
!      doi:10.1029/2002JD003322, 2003.                                     !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    aer's revision history:                                               !
!     this version of rrtmg_lw has been modified from rrtm_lw to use a     !
!     reduced set of g-points for application to gcms.                     !
!                                                                          !
! --  original version (derived from rrtm_lw), reduction of g-points,      !
!     other revisions for use with gcms.                                   !
!        1999: m. j. iacono, aer, inc.                                     !
! --  adapted for use with ncar/cam3.                                      !
!        may 2004: m. j. iacono, aer, inc.                                 !
! --  revised to add mcica capability.                                     !
!        nov 2005: m. j. iacono, aer, inc.                                 !
! --  conversion to f90 formatting for consistency with rrtmg_sw.          !
!        feb 2007: m. j. iacono, aer, inc.                                 !
! --  modifications to formatting to use assumed-shape arrays.             !
!        aug 2007: m. j. iacono, aer, inc.                                 !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    ncep modifications history log:                                       !
!                                                                          !
!       nov 1999,  ken campana       -- received the original code from    !
!                    aer (1998 ncar ccm version), updated to link up with  !
!                    ncep mrf model                                        !
!       jun 2000,  ken campana       -- added option to switch random and  !
!                    maximum/random cloud overlap                          !
!           2001,  shrinivas moorthi -- further updates for mrf model      !
!       may 2001,  yu-tai hou        -- updated on trace gases and cloud   !
!                    property based on rrtm_v3.0 codes.                    !
!       dec 2001,  yu-tai hou        -- rewritten code into fortran 90 std !
!                    set ncep radiation structure standard that contains   !
!                    three plug-in compatable fortran program files:       !
!                    'radlw_param.f', 'radlw_datatb.f', 'radlw_main.f'     !
!                    fixed bugs in subprograms taugb14, taugb2, etc. added !
!                    out-of-bounds protections. (a detailed note of        !
!                    up_to_date modifications/corrections by ncep was sent !
!                    to aer in 2002)                                       !
!       jun 2004,  yu-tai hou        -- added mike iacono's apr 2004       !
!                    modification of variable diffusivity angles.          !
!       apr 2005,  yu-tai hou        -- minor modifications on module      !
!                    structures include rain/snow effect (this version of  !
!                    code was given back to aer in jun 2006)               !
!       mar 2007,  yu-tai hou        -- added aerosol effect for ncep      !
!                    models using the generallized aerosol optical property!
!                    scheme for gfs model.                                 !
!       apr 2007,  yu-tai hou        -- added spectral band heating as an  !
!                    optional output to support the 500 km gfs model's     !
!                    upper stratospheric radiation calculations. and       !
!                    restructure optional outputs for easy access by       !
!                    different models.                                     !
!       oct 2008,  yu-tai hou        -- modified to include new features   !
!                    from aer's newer release v4.4-v4.7, including the     !
!                    mcica sub-grid cloud option. add rain/snow optical    !
!                    properties support to cloudy sky calculations.        !
!                    correct errors in mcica cloud optical properties for  !
!                    ebert & curry scheme (ilwcice=1) that needs band      !
!                    index conversion. simplified and unified sw and lw    !
!                    sub-column cloud subroutines into one module by using !
!                    optional parameters.                                  !
!       mar 2009,  yu-tai hou        -- replaced the original random number!
!                    generator coming from the original code with ncep w3  !
!                    library to simplify the program and moved sub-column  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       oct 2009,  yu-tai hou        -- modified subrtines "cldprop" and   !
!                    "rlwinit" according updats from aer's rrtmg_lw v4.8.  !
!       nov 2009,  yu-tai hou        -- modified subrtine "taumol" according
!                    updats from aer's rrtmg_lw version 4.82. notice the   !
!                    cloud ice/liquid are assumed as in-cloud quantities,  !
!                    not as grid averaged quantities.                      !
!       jun 2010,  yu-tai hou        -- optimized code to improve efficiency
!       apr 2012,  b. ferrier and y. hou -- added conversion factor to fu's!
!                    cloud-snow optical property scheme.                   !
!       nov 2012,  yu-tai hou        -- modified control parameters thru   !
!                     module 'physparam'.                                  !
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!

!========================================!
      module rrtmg_lw                    !
!........................................!
!
      use physparam,        only : ilwrate, ilwrgas, ilwcliq, ilwcice,  &
     &                             isubclw, icldflg, iovrlw,  ivflip,   &
     &                             kind_phys
      use physcons,         only : con_g, con_cp, con_avgd, con_amd,    &
     &                             con_amw, con_amo3
      use mersenne_twister, only : random_setseed, random_number,       &
     &                             random_stat

      use module_radlw_parameters
!
      use module_radlw_avplank, only : totplnk
      use module_radlw_ref,     only : preflog, tref, chi_mls
!
      implicit none
!
      private
!
!  ...  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGLW='NCEP LW v5.1  Nov 2012 -RRTMG-LW v4.82  '
!    &   VTAGLW='NCEP LW v5.0  Aug 2012 -RRTMG-LW v4.82  '
!    &   VTAGLW='RRTMG-LW v4.82  Nov 2009  '
!    &   VTAGLW='RRTMG-LW v4.8   Oct 2009  '
!    &   VTAGLW='RRTMG-LW v4.71  Mar 2009  '
!    &   VTAGLW='RRTMG-LW v4.4   Oct 2008  '
!    &   VTAGLW='RRTM-LW v2.3g   Mar 2007  '
!    &   VTAGLW='RRTM-LW v2.3g   Apr 2004  '

!  ---  constant values
      real (kind=kind_phys), parameter :: eps     = 1.0e-6
      real (kind=kind_phys), parameter :: oneminus= 1.0-eps
      real (kind=kind_phys), parameter :: cldmin  = tiny(cldmin)
      real (kind=kind_phys), parameter :: bpade   = 1.0/0.278  ! pade approx constant
      real (kind=kind_phys), parameter :: stpfac  = 296.0/1013.0
      real (kind=kind_phys), parameter :: wtdiff  = 0.5        ! weight for radiance to flux conversion
      real (kind=kind_phys), parameter :: tblint  = ntbl       ! lookup table conversion factor
      real (kind=kind_phys), parameter :: f_zero  = 0.0
      real (kind=kind_phys), parameter :: f_one   = 1.0

!  ...  atomic weights for conversion from mass to volume mixing ratios
      real (kind=kind_phys), parameter :: amdw    = con_amd/con_amw
      real (kind=kind_phys), parameter :: amdo3   = con_amd/con_amo3

!  ...  band indices
      integer, dimension(nbands) :: nspa, nspb

      data nspa / 1, 1, 9, 9, 9, 1, 9, 1, 9, 1, 1, 9, 9, 1, 9, 9 /
      data nspb / 1, 1, 5, 5, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0 /

!  ...  band wavenumber intervals
!     real (kind=kind_phys) :: wavenum1(nbands), wavenum2(nbands)
!     data wavenum1/                                                    &
!    &         10.,  350.,  500.,  630.,  700.,  820.,  980., 1080.,    &
!err &       1180., 1390., 1480., 1800., 2080., 2250., 2390., 2600. /
!    &       1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600. /
!     data wavenum2/                                                    &
!    &        350.,  500.,  630.,  700.,  820.,  980., 1080., 1180.,    &
!err &       1390., 1480., 1800., 2080., 2250., 2390., 2600., 3250. /
!    &       1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250. /
!     real (kind=kind_phys) :: delwave(nbands)
!     data delwave / 340., 150., 130.,  70., 120., 160., 100., 100.,    &
!    &               210.,  90., 320., 280., 170., 130., 220., 650. /

!  ---  reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
!       and 1.80) as a function of total column water vapor.  the function
!       has been defined to minimize flux and cooling rate errors in these bands
!       over a wide range of precipitable water values.
      real (kind=kind_phys), dimension(nbands) :: a0, a1, a2

      data a0 / 1.66,  1.55,  1.58,  1.66,  1.54, 1.454,  1.89,  1.33,  &
     &         1.668,  1.66,  1.66,  1.66,  1.66,  1.66,  1.66,  1.66 /
      data a1 / 0.00,  0.25,  0.22,  0.00,  0.13, 0.446, -0.10,  0.40,  &
     &        -0.006,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /
      data a2 / 0.00, -12.0, -11.7,  0.00, -0.72,-0.243,  0.19,-0.062,  &
     &         0.414,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00 /

!! ---  logical flags for optional output fields

      logical :: lhlwb  = .false.
      logical :: lhlw0  = .false.
      logical :: lflxprf= .false.

!  ---  those data will be set up only once by "rlwinit"

!  ...  fluxfac, heatfac are factors for fluxes (in w/m**2) and heating
!       rates (in k/day, or k/sec set by subroutine 'rlwinit')
!       semiss0 are default surface emissivity for each bands

      real (kind=kind_phys) :: fluxfac, heatfac, semiss0(nbands)
      data semiss0(:) / nbands*1.0 /

      real (kind=kind_phys) :: tau_tbl(0:ntbl)  !< clr-sky opt dep (for cldy transfer)
      real (kind=kind_phys) :: exp_tbl(0:ntbl)  !< transmittance lookup table
      real (kind=kind_phys) :: tfn_tbl(0:ntbl)  !< tau transition function; i.e. the
                                                !< transition of planck func from mean lyr
                                                !< temp to lyr boundary temp as a func of
                                                !< opt dep. "linear in tau" method is used.

!  ---  the following variables are used for sub-column cloud scheme

      integer, parameter :: ipsdlw0 = ngptlw     ! initial permutation seed

!  ---  public accessable subprograms

      public rrtmg_lw_init, rrtmg_lw_run, rrtmg_lw_finalize, rlwinit


! ================
      contains
! ================

         subroutine rrtmg_lw_init ()
         end subroutine rrtmg_lw_init

!> \defgroup module_radlw_main GFS radlw Main
!! \brief This module includes NCEP's modifications of the RRTMG-LW radiation
!! code from AER.
!!
!! The RRTM-LW package includes three files:
!! - radlw_param.f, which contains:
!!  - module_radlw_parameters: band parameters set up
!! - radlw_datatb.f, which contains modules:
!!  - module_radlw_avplank: plank flux data
!!  - module_radlw_ref: reference temperature and pressure
!!  - module_radlw_cldprlw: cloud property coefficients
!!  - module_radlw_kgbnn: absorption coeffients for 16 bands, where nn = 01-16
!! - radlw_main.f, which contains:
!!  - rrtmg_lw_run(): the main LW radiation routine
!!  - rlwinit(): the initialization routine
!!
!!\version NCEP LW v5.1  Nov 2012 -RRTMG-LW v4.82
!!
!!\copyright  2002-2007, Atmospheric & Environmental Research, Inc. (AER).
!!  This software may be used, copied, or redistributed as long as it is
!!  not sold and this copyright notice is reproduced on each copy made.
!!  This model is provided as is without any express or implied warranties.
!!  (http://www.rtweb.aer.com/)
!! \section arg_table_rrtmg_lw_run Argument Table
!! | local_name      | standard_name                                                                                 | long_name                                                 | units   | rank | type        |    kind   | intent | optional |
!! |-----------------|-----------------------------------------------------------------------------------------------|-----------------------------------------------------------|---------|------|-------------|-----------|--------|----------|
!! | plyr            | air_pressure_at_layer_for_radiation_in_hPa                                                    | air pressure layer                                        | hPa     |    2 | real        | kind_phys | in     | F        |
!! | plvl            | air_pressure_at_interface_for_radiation_in_hPa                                                | air pressure level                                        | hPa     |    2 | real        | kind_phys | in     | F        |
!! | tlyr            | air_temperature_at_layer_for_radiation                                                        | air temperature layer                                     | K       |    2 | real        | kind_phys | in     | F        |
!! | tlvl            | air_temperature_at_interface_for_radiation                                                    | air temperature level                                     | K       |    2 | real        | kind_phys | in     | F        |
!! | qlyr            | water_vapor_specific_humidity_at_layer_for_radiation                                          | specific humidity layer                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | olyr            | ozone_concentration_at_layer_for_radiation                                                    | ozone concentration layer                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_co2      | volume_mixing_ratio_co2                                                                       | volume mixing ratio co2                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_n2o      | volume_mixing_ratio_n2o                                                                       | volume mixing ratio no2                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_ch4      | volume_mixing_ratio_ch4                                                                       | volume mixing ratio ch4                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_o2       | volume_mixing_ratio_o2                                                                        | volume mixing ratio o2                                    | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_co       | volume_mixing_ratio_co                                                                        | volume mixing ratio co                                    | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_cfc11    | volume_mixing_ratio_cfc11                                                                     | volume mixing ratio cfc11                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_cfc12    | volume_mixing_ratio_cfc12                                                                     | volume mixing ratio cfc12                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_cfc22    | volume_mixing_ratio_cfc22                                                                     | volume mixing ratio cfc22                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_ccl4     | volume_mixing_ratio_ccl4                                                                      | volume mixing ratio ccl4                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | icseed          | seed_random_numbers_lw                                                                        | seed for random number generation for longwave radiation  | none    |    1 | integer     |           | in     | F        |
!! | aeraod          | aerosol_optical_depth_for_longwave_bands_01-16                                                | aerosol optical depth for longwave bands 01-16            | none    |    3 | real        | kind_phys | in     | F        |
!! | aerssa          | aerosol_single_scattering_albedo_for_longwave_bands_01-16                                     | aerosol single scattering albedo for longwave bands 01-16 | frac    |    3 | real        | kind_phys | in     | F        |
!! | sfemis          | surface_longwave_emissivity                                                                   | surface emissivity                                        | frac    |    1 | real        | kind_phys | in     | F        |
!! | sfgtmp          | surface_ground_temperature_for_radiation                                                      | surface ground temperature for radiation                  | K       |    1 | real        | kind_phys | in     | F        |
!! | npts            | horizontal_loop_extent                                                                        | horizontal dimension                                      | count   |    0 | integer     |           | in     | F        |
!! | nlay            | adjusted_vertical_layer_dimension_for_radiation                                               | number of vertical layers for radiation                   | count   |    0 | integer     |           | in     | F        |
!! | nlp1            | adjusted_vertical_level_dimension_for_radiation                                               | number of vertical levels for radiation                   | count   |    0 | integer     |           | in     | F        |
!! | lprnt           | flag_print                                                                                    | flag to print                                             | flag    |    0 | logical     |           | in     | F        |
!! | cld_cf          | total_cloud_fraction                                                                          | total cloud fraction                                      | frac    |    2 | real        | kind_phys | in     | F        |
!! | lslwr           | flag_to_calc_lw                                                                               | flag to calculate LW irradiances                          | flag    |    0 | logical     |           | in     | F        |
!! | hlwc            | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                    | longwave total sky heating rate                           | K s-1   |    2 | real        | kind_phys | inout  | F        |
!! | topflx          | lw_fluxes_top_atmosphere                                                                      | longwave total sky fluxes at the top of the atm           | W m-2   |    1 | topflw_type |           | inout  | F        |
!! | sfcflx          | lw_fluxes_sfc                                                                                 | longwave total sky fluxes at the Earth surface            | W m-2   |    1 | sfcflw_type |           | inout  | F        |
!! | hlw0            | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step | longwave clear sky heating rate                           | K s-1   |    2 | real        | kind_phys | inout  | T        |
!! | hlwb            | lw_heating_rate_spectral                                                                      | longwave total sky heating rate (spectral)                | K s-1   |    3 | real        | kind_phys | inout  | T        |
!! | flxprf          | lw_fluxes                                                                                     | lw fluxes total sky / csk and up / down at levels         | W m-2   |    2 | proflw_type |           | inout  | T        |
!! | cld_lwp         | cloud_liquid_water_path                                                                       | cloud liquid water path                                   | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_liq     | mean_effective_radius_for_liquid_cloud                                                        | mean effective radius for liquid cloud                    | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_iwp         | cloud_ice_water_path                                                                          | cloud ice water path                                      | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_ice     | mean_effective_radius_for_ice_cloud                                                           | mean effective radius for ice cloud                       | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_rwp         | cloud_rain_water_path                                                                         | cloud ice water path                                      | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_rain    | mean_effective_radius_for_rain_drop                                                           | mean effective radius for rain drop                       | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_swp         | cloud_snow_water_path                                                                         | cloud snow water path                                     | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_snow    | mean_effective_radius_for_snow_flake                                                          | mean effective radius for snow flake                      | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_od          | cloud_optical_depth                                                                           | cloud optical depth                                       | none    |    2 | real        | kind_phys | in     | T        |
!! | errmsg          | error_message                                                                                 | error message for error handling in CCPP                  | none    |    0 | character   | len=*     | out    | F        |
!! | errflg          | error_flag                                                                                    | error flag for error handling in CCPP                     | flag    |    0 | integer     |           | out    | F        |
!!
!> \section gen_lwrad RRTMG Longwave Radiation Scheme General Algorithm
!> @{
      subroutine rrtmg_lw_run                                           &
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr_co2, gasvmr_n2o,      &   !  ---  inputs
     &       gasvmr_ch4, gasvmr_o2, gasvmr_co, gasvmr_cfc11,            &
     &       gasvmr_cfc12, gasvmr_cfc22, gasvmr_ccl4,                   &
     &       icseed,aeraod,aerssa,sfemis,sfgtmp,                        &
     &       npts, nlay, nlp1, lprnt, cld_cf, lslwr,                    &
     &       hlwc,topflx,sfcflx,                                        &   !  ---  outputs
     &       HLW0,HLWB,FLXPRF,                                          &   !  ---  optional
     &       cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice,                &
     &       cld_rwp,cld_ref_rain, cld_swp, cld_ref_snow,               &
     &       cld_od, errmsg, errflg                                       &
     &     )


!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!     plyr (npts,nlay) : layer mean pressures (mb)                      !
!     plvl (npts,nlp1) : interface pressures (mb)                       !
!     tlyr (npts,nlay) : layer mean temperature (k)                     !
!     tlvl (npts,nlp1) : interface temperatures (k)                     !
!     qlyr (npts,nlay) : layer specific humidity (gm/gm)   *see inside  !
!     olyr (npts,nlay) : layer ozone concentration (gm/gm) *see inside  !
!     gasvmr(npts,nlay,:): atmospheric gases amount:                    !
!                       (check module_radiation_gases for definition)   !
!       gasvmr(:,:,1)  -   co2 volume mixing ratio                      !
!       gasvmr(:,:,2)  -   n2o volume mixing ratio                      !
!       gasvmr(:,:,3)  -   ch4 volume mixing ratio                      !
!       gasvmr(:,:,4)  -   o2  volume mixing ratio                      !
!       gasvmr(:,:,5)  -   co  volume mixing ratio                      !
!       gasvmr(:,:,6)  -   cfc11 volume mixing ratio                    !
!       gasvmr(:,:,7)  -   cfc12 volume mixing ratio                    !
!       gasvmr(:,:,8)  -   cfc22 volume mixing ratio                    !
!       gasvmr(:,:,9)  -   ccl4  volume mixing ratio                    !
!     clouds(npts,nlay,:): layer cloud profiles:                        !
!                       (check module_radiation_clouds for definition)  !
!                ---  for  ilwcliq > 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer in-cloud liq water path   (g/m**2)     !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)     !
!       clouds(:,:,4)  -   layer in-cloud ice water path   (g/m**2)     !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)     !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)     !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)     !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)     !
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)     !
!                ---  for  ilwcliq = 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer cloud optical depth                    !
!       clouds(:,:,3)  -   layer cloud single scattering albedo         !
!       clouds(:,:,4)  -   layer cloud asymmetry factor                 !
!     icseed(npts)   : auxiliary special cloud related array            !
!                      when module variable isubclw=2, it provides      !
!                      permutation seed for each column profile that    !
!                      are used for generating random numbers.          !
!                      when isubclw /=2, it will not be used.           !
!     aerosols(npts,nlay,nbands,:) : aerosol optical properties         !
!                       (check module_radiation_aerosols for definition)!
!        (:,:,:,1)     - optical depth                                  !
!        (:,:,:,2)     - single scattering albedo                       !
!        (:,:,:,3)     - asymmetry parameter                            !
!     sfemis (npts)  : surface emissivity                               !
!     sfgtmp (npts)  : surface ground temperature (k)                   !
!     npts           : total number of horizontal points                !
!     nlay, nlp1     : total number of vertical layers, levels          !
!     lprnt          : cntl flag for diagnostic print out               !
!                                                                       !
!  output variables:                                                    !
!     hlwc  (npts,nlay): total sky heating rate (k/day or k/sec)        !
!     topflx(npts)     : radiation fluxes at top, component:            !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux at top (w/m2)          !
!        upfx0           - clear sky upward flux at top (w/m2)          !
!     sfcflx(npts)     : radiation fluxes at sfc, component:            !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux at sfc (w/m2)          !
!        upfx0           - clear sky upward flux at sfc (w/m2)          !
!        dnfxc           - total sky downward flux at sfc (w/m2)        !
!        dnfx0           - clear sky downward flux at sfc (w/m2)        !
!                                                                       !
!! optional output variables:                                           !
!     hlwb(npts,nlay,nbands): spectral band total sky heating rates     !
!     hlw0  (npts,nlay): clear sky heating rate (k/day or k/sec)        !
!     flxprf(npts,nlp1): level radiative fluxes (w/m2), components:     !
!                        (check module_radlw_paramters for definition)  !
!        upfxc           - total sky upward flux                        !
!        dnfxc           - total sky dnward flux                        !
!        upfx0           - clear sky upward flux                        !
!        dnfx0           - clear sky dnward flux                        !
!                                                                       !
!  external module variables:  (in physparam)                            !
!   ilwrgas - control flag for rare gases (ch4,n2o,o2,cfcs, etc.)       !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   ilwcliq - control flag for liq-cloud optical properties             !
!           =0: input cloud optical depth, ignor ilwcice                !
!           =1: input cld liqp & reliq, hu & stamnes (1993)             !
!           =2: not used                                                !
!   ilwcice - control flag for ice-cloud optical properties             !
!           *** if ilwcliq==0, ilwcice is ignored                       !
!           =1: input cld icep & reice, ebert & curry (1997)            !
!           =2: input cld icep & reice, streamer (1996)                 !
!           =3: input cld icep & reice, fu (1998)                       !
!   isubclw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   iovrlw  - cloud overlapping control flag                            !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud (used for isubclw>0 only)         !
!   ivflip  - control flag for vertical index direction                 !
!           =0: vertical index from toa to surface                      !
!           =1: vertical index from surface to toa                      !
!                                                                       !
!  module parameters, control variables:                                !
!     nbands           - number of longwave spectral bands              !
!     maxgas           - maximum number of absorbing gaseous            !
!     maxxsec          - maximum number of cross-sections               !
!     ngptlw           - total number of g-point subintervals           !
!     ng##             - number of g-points in band (##=1-16)           !
!     ngb(ngptlw)      - band indices for each g-point                  !
!     bpade            - pade approximation constant (1/0.278)          !
!     nspa,nspb(nbands)- number of lower/upper ref atm's per band       !
!     delwave(nbands)  - longwave band width (wavenumbers)              !
!     ipsdlw0          - permutation seed for mcica sub-col clds        !
!                                                                       !
!  major local variables:                                               !
!     pavel  (nlay)         - layer pressures (mb)                      !
!     delp   (nlay)         - layer pressure thickness (mb)             !
!     tavel  (nlay)         - layer temperatures (k)                    !
!     tz     (0:nlay)       - level (interface) temperatures (k)        !
!     semiss (nbands)       - surface emissivity for each band          !
!     wx     (nlay,maxxsec) - cross-section molecules concentration     !
!     coldry (nlay)         - dry air column amount                     !
!                                   (1.e-20*molecules/cm**2)            !
!     cldfrc (0:nlp1)       - layer cloud fraction                      !
!     taucld (nbands,nlay)  - layer cloud optical depth for each band   !
!     cldfmc (ngptlw,nlay)  - layer cloud fraction for each g-point     !
!     tauaer (nbands,nlay)  - aerosol optical depths                    !
!     fracs  (ngptlw,nlay)  - planck fractions                          !
!     tautot (ngptlw,nlay)  - total optical depths (gaseous+aerosols)   !
!     colamt (nlay,maxgas)  - column amounts of absorbing gases         !
!                             1-maxgas are for watervapor, carbon       !
!                             dioxide, ozone, nitrous oxide, methane,   !
!                             oxigen, carbon monoxide, respectively     !
!                             (molecules/cm**2)                         !
!     pwvcm                 - column precipitable water vapor (cm)      !
!     secdiff(nbands)       - variable diffusivity angle defined as     !
!                             an exponential function of the column     !
!                             water amount in bands 2-3 and 5-9.        !
!                             this reduces the bias of several w/m2 in  !
!                             downward surface flux in high water       !
!                             profiles caused by using the constant     !
!                             diffusivity angle of 1.66.         (mji)  !
!     facij  (nlay)         - indicator of interpolation factors        !
!                             =0/1: indicate lower/higher temp & height !
!     selffac(nlay)         - scale factor for self-continuum, equals   !
!                          (w.v. density)/(atm density at 296K,1013 mb) !
!     selffrac(nlay)        - factor for temp interpolation of ref      !
!                             self-continuum data                       !
!     indself(nlay)         - index of the lower two appropriate ref    !
!                             temp for the self-continuum interpolation !
!     forfac (nlay)         - scale factor for w.v. foreign-continuum   !
!     forfrac(nlay)         - factor for temp interpolation of ref      !
!                             w.v. foreign-continuum data               !
!     indfor (nlay)         - index of the lower two appropriate ref    !
!                             temp for the foreign-continuum interp     !
!     laytrop               - tropopause layer index at which switch is !
!                             made from one conbination kew species to  !
!                             another.                                  !
!     jp(nlay),jt(nlay),jt1(nlay)                                       !
!                           - lookup table indexes                      !
!     totuflux(0:nlay)      - total-sky upward longwave flux (w/m2)     !
!     totdflux(0:nlay)      - total-sky downward longwave flux (w/m2)   !
!     htr(nlay)             - total-sky heating rate (k/day or k/sec)   !
!     totuclfl(0:nlay)      - clear-sky upward longwave flux (w/m2)     !
!     totdclfl(0:nlay)      - clear-sky downward longwave flux (w/m2)   !
!     htrcl(nlay)           - clear-sky heating rate (k/day or k/sec)   !
!     fnet    (0:nlay)      - net longwave flux (w/m2)                  !
!     fnetc   (0:nlay)      - clear-sky net longwave flux (w/m2)        !
!                                                                       !
!                                                                       !
!  ======================    end of definitions    ===================  !

!  ---  inputs:
      integer, intent(in) :: npts, nlay, nlp1
      integer, intent(in) :: icseed(npts)

      logical,  intent(in) :: lprnt

      real (kind=kind_phys), dimension(npts,nlp1), intent(in) :: plvl,  &
     &       tlvl
      real (kind=kind_phys), dimension(npts,nlay), intent(in) :: plyr,  &
     &       tlyr, qlyr, olyr

      real (kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_co2,&
     &     gasvmr_n2o, gasvmr_ch4, gasvmr_o2, gasvmr_co, gasvmr_cfc11,  &
     &     gasvmr_cfc12, gasvmr_cfc22, gasvmr_ccl4

      real (kind=kind_phys), dimension(npts,nlay),intent(in):: cld_cf
      real (kind=kind_phys), dimension(npts,nlay),intent(in),optional:: &
     &       cld_lwp, cld_ref_liq,  cld_iwp, cld_ref_ice,               &
     &       cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow,              &
     &       cld_od

      real (kind=kind_phys), dimension(npts), intent(in) :: sfemis,     &
     &       sfgtmp

      real (kind=kind_phys), dimension(npts,nlay,nbands),intent(in)::   &
     &       aeraod, aerssa

!  ---  outputs:
      real (kind=kind_phys), dimension(npts,nlay), intent(inout) :: hlwc

      type (topflw_type),    dimension(npts), intent(inout) :: topflx
      type (sfcflw_type),    dimension(npts), intent(inout) :: sfcflx

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!! ---  optional outputs:
      real (kind=kind_phys), dimension(npts,nlay,nbands),optional,      &
     &       intent(inout) :: hlwb
      real (kind=kind_phys), dimension(npts,nlay),       optional,      &
     &       intent(inout) :: hlw0
      type (proflw_type),    dimension(npts,nlp1),       optional,      &
     &       intent(inout) :: flxprf
      logical, intent(in) :: lslwr

!  ---  locals:
      real (kind=kind_phys), dimension(0:nlp1) :: cldfrc

      real (kind=kind_phys), dimension(0:nlay) :: totuflux, totdflux,   &
     &       totuclfl, totdclfl, tz

      real (kind=kind_phys), dimension(nlay)   :: htr, htrcl

      real (kind=kind_phys), dimension(nlay)   :: pavel, tavel, delp,   &
     &       clwp, ciwp, relw, reiw, cda1, cda2, cda3, cda4,            &
     &       coldry, colbrd, h2ovmr, o3vmr, fac00, fac01, fac10, fac11, &
     &       selffac, selffrac, forfac, forfrac, minorfrac, scaleminor, &
     &       scaleminorn2, temcol

      real (kind=kind_phys), dimension(nbands,0:nlay) :: pklev, pklay

      real (kind=kind_phys), dimension(nlay,nbands) :: htrb
      real (kind=kind_phys), dimension(nbands,nlay) :: taucld, tauaer
      real (kind=kind_phys), dimension(ngptlw,nlay) :: fracs, tautot,   &
     &       cldfmc

      real (kind=kind_phys), dimension(nbands) :: semiss, secdiff

!  ---  column amount of absorbing gases:
!       (:,m) m = 1-h2o, 2-co2, 3-o3, 4-n2o, 5-ch4, 6-o2, 7-co
      real (kind=kind_phys) :: colamt(nlay,maxgas)

!  ---  column cfc cross-section amounts:
!       (:,m) m = 1-ccl4, 2-cfc11, 3-cfc12, 4-cfc22
      real (kind=kind_phys) :: wx(nlay,maxxsec)

!  ---  reference ratios of binary species parameter in lower atmosphere:
!       (:,m,:) m = 1-h2o/co2, 2-h2o/o3, 3-h2o/n2o, 4-h2o/ch4, 5-n2o/co2, 6-o3/co2
      real (kind=kind_phys) :: rfrate(nlay,nrates,2)

      real (kind=kind_phys) :: tem0, tem1, tem2, pwvcm, summol, stemp

      integer, dimension(npts) :: ipseed
      integer, dimension(nlay) :: jp, jt, jt1, indself, indfor, indminor
      integer                  :: laytrop, iplon, i, j, k, k1
      logical :: lcf1

!
!===> ... begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
      if (.not. lslwr) return

!  --- ...  initialization

      lhlwb  = present ( hlwb )
      lhlw0  = present ( hlw0 )
      lflxprf= present ( flxprf )

      colamt(:,:) = f_zero

!! --- check for optional input arguments, depending on cloud method
      if (ilwcliq > 0) then    ! use prognostic cloud method
        if ( .not.present(cld_lwp) .or. .not.present(cld_ref_liq) .or.  &
     &       .not.present(cld_iwp) .or. .not.present(cld_ref_ice) .or.  &
     &       .not.present(cld_rwp) .or. .not.present(cld_ref_rain) .or. &
     &       .not.present(cld_swp) .or. .not.present(cld_ref_snow)) then
          write(errmsg,'(*(a))')                                        &
     &               'Logic error: ilwcliq>0 requires the following',   &
     &               ' optional arguments to be present:',              &
     &               ' cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice,',    &
     &               ' cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow'
          errflg = 1
          return
        end if
      else                     ! use diagnostic cloud method
        if ( .not.present(cld_od) ) then
          write(errmsg,'(*(a))')                                        &
     &               'Logic error: ilwcliq<=0 requires the following',  &
     &               ' optional argument to be present: cld_od'
          errflg = 1
          return
        end if
      endif                    ! end if_ilwcliq

!> -# Change random number seed value for each radiation invocation
!!    (isubclw =1 or 2).

      if     ( isubclw == 1 ) then     ! advance prescribed permutation seed
        do i = 1, npts
          ipseed(i) = ipsdlw0 + i
        enddo
      elseif ( isubclw == 2 ) then     ! use input array of permutaion seeds
        do i = 1, npts
          ipseed(i) = icseed(i)
        enddo
      endif

!     if ( lprnt ) then
!       print *,'  In rrtmg_lw, isubclw, ipsdlw0,ipseed =',             &
!    &          isubclw, ipsdlw0, ipseed
!     endif

!  --- ...  loop over horizontal npts profiles

      lab_do_iplon : do iplon = 1, npts

!> -# Read surface emissivity.
        if (sfemis(iplon) > eps .and. sfemis(iplon) <= 1.0) then  ! input surface emissivity
          do j = 1, nbands
            semiss(j) = sfemis(iplon)
          enddo
        else                                                      ! use default values
          do j = 1, nbands
            semiss(j) = semiss0(j)
          enddo
        endif

        stemp = sfgtmp(iplon)          ! surface ground temp

!> -# Prepare atmospheric profile for use in rrtm.
!           the vertical index of internal array is from surface to top

!  --- ...  molecular amounts are input or converted to volume mixing ratio
!           and later then converted to molecular amount (molec/cm2) by the
!           dry air column coldry (in molec/cm2) which is calculated from the
!           layer pressure thickness (in mb), based on the hydrostatic equation
!  --- ...  and includes a correction to account for h2o in the layer.

        if (ivflip == 0) then       ! input from toa to sfc

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd
          tz(0) = tlvl(iplon,nlp1)

          do k = 1, nlay
            k1 = nlp1 - k
            pavel(k)= plyr(iplon,k1)
            delp(k) = plvl(iplon,k1+1) - plvl(iplon,k1)
            tavel(k)= tlyr(iplon,k1)
            tz(k)   = tlvl(iplon,k1)

!> -# Set absorber amount for h2o, co2, and o3.

!test use
!           h2ovmr(k)= max(f_zero,qlyr(iplon,k1)*amdw)                  ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qlyr(iplon,k1))                       ! input vol mixing ratio
!           o3vmr (k)= max(f_zero,olyr(iplon,k1))                       ! input vol mixing ratio
!ncep model use
            h2ovmr(k)= max(f_zero,qlyr(iplon,k1)                        &
     &                           *amdw/(f_one-qlyr(iplon,k1)))          ! input specific humidity
            o3vmr (k)= max(f_zero,olyr(iplon,k1)*amdo3)                 ! input mass mixing ratio

!  --- ...  tem0 is the molecular weight of moist air
            tem0 = (f_one - h2ovmr(k))*con_amd + h2ovmr(k)*con_amw
            coldry(k) = tem2*delp(k) / (tem1*tem0*(f_one+h2ovmr(k)))
            temcol(k) = 1.0e-12 * coldry(k)

            colamt(k,1) = max(f_zero,    coldry(k)*h2ovmr(k))          ! h2o
            colamt(k,2) = max(temcol(k), coldry(k)*gasvmr_co2(iplon,k1)) ! co2
            colamt(k,3) = max(temcol(k), coldry(k)*o3vmr(k))           ! o3
          enddo

!> -# Set up column amount for rare gases n2o,ch4,o2,co,ccl4,cf11,cf12,
!!    cf22, convert from volume mixing ratio to molec/cm2 based on
!!    coldry (scaled to 1.0e-20).

          if (ilwrgas > 0) then
            do k = 1, nlay
              k1 = nlp1 - k
              colamt(k,4)=max(temcol(k), coldry(k)*gasvmr_n2o(iplon,k1))  ! n2o
              colamt(k,5)=max(temcol(k), coldry(k)*gasvmr_ch4(iplon,k1))  ! ch4
              colamt(k,6)=max(f_zero,    coldry(k)*gasvmr_o2(iplon,k1))  ! o2
              colamt(k,7)=max(f_zero,    coldry(k)*gasvmr_co(iplon,k1))  ! co

              wx(k,1) = max( f_zero, coldry(k)*gasvmr_ccl4(iplon,k1) )   ! ccl4
              wx(k,2) = max( f_zero, coldry(k)*gasvmr_cfc11(iplon,k1) )   ! cf11
              wx(k,3) = max( f_zero, coldry(k)*gasvmr_cfc12(iplon,k1) )   ! cf12
              wx(k,4) = max( f_zero, coldry(k)*gasvmr_cfc22(iplon,k1) )   ! cf22
            enddo
          else
            do k = 1, nlay
              colamt(k,4) = f_zero     ! n2o
              colamt(k,5) = f_zero     ! ch4
              colamt(k,6) = f_zero     ! o2
              colamt(k,7) = f_zero     ! co

              wx(k,1) = f_zero
              wx(k,2) = f_zero
              wx(k,3) = f_zero
              wx(k,4) = f_zero
            enddo
          endif

!> -# Set aerosol optical properties.

          do k = 1, nlay
            k1 = nlp1 - k
            do j = 1, nbands
              tauaer(j,k) = aeraod(iplon,k1,j)                          &
     &                    * (f_one - aerssa(iplon,k1,j))
            enddo
          enddo

!> -# Read cloud optical properties.
          if (ilwcliq > 0) then    ! use prognostic cloud method
            do k = 1, nlay
              k1 = nlp1 - k
              cldfrc(k)= cld_cf(iplon,k1)
              clwp(k)  = cld_lwp(iplon,k1)
              relw(k)  = cld_ref_liq(iplon,k1)
              ciwp(k)  = cld_iwp(iplon,k1)
              reiw(k)  = cld_ref_ice(iplon,k1)
              cda1(k)  = cld_rwp(iplon,k1)
              cda2(k)  = cld_ref_rain(iplon,k1)
              cda3(k)  = cld_swp(iplon,k1)
              cda4(k)  = cld_ref_snow(iplon,k1)
            enddo
          else                       ! use diagnostic cloud method
            do k = 1, nlay
              k1 = nlp1 - k
              cldfrc(k)= cld_cf(iplon,k1)
              cda1(k)  = cld_od(iplon,k1)
            enddo
          endif                      ! end if_ilwcliq

          cldfrc(0)    = f_one       ! padding value only
          cldfrc(nlp1) = f_zero      ! padding value only

!> -# Compute precipitable water vapor for diffusivity angle adjustments.

          tem1 = f_zero
          tem2 = f_zero
          do k = 1, nlay
            tem1 = tem1 + coldry(k) + colamt(k,1)
            tem2 = tem2 + colamt(k,1)
          enddo

          tem0 = 10.0 * tem2 / (amdw * tem1 * con_g)
          pwvcm = tem0 * plvl(iplon,nlp1)

        else                        ! input from sfc to toa

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd
          tz(0) = tlvl(iplon,1)

          do k = 1, nlay
            pavel(k)= plyr(iplon,k)
            delp(k) = plvl(iplon,k) - plvl(iplon,k+1)
            tavel(k)= tlyr(iplon,k)
            tz(k)   = tlvl(iplon,k+1)

!  --- ...  set absorber amount
!test use
!           h2ovmr(k)= max(f_zero,qlyr(iplon,k)*amdw)                   ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qlyr(iplon,k))                        ! input vol mixing ratio
!           o3vmr (k)= max(f_zero,olyr(iplon,k))                        ! input vol mixing ratio
!ncep model use
            h2ovmr(k)= max(f_zero,qlyr(iplon,k)                         &
     &                           *amdw/(f_one-qlyr(iplon,k)))           ! input specific humidity
            o3vmr (k)= max(f_zero,olyr(iplon,k)*amdo3)                  ! input mass mixing ratio

!  --- ...  tem0 is the molecular weight of moist air
            tem0 = (f_one - h2ovmr(k))*con_amd + h2ovmr(k)*con_amw
            coldry(k) = tem2*delp(k) / (tem1*tem0*(f_one+h2ovmr(k)))
            temcol(k) = 1.0e-12 * coldry(k)

            colamt(k,1) = max(f_zero,    coldry(k)*h2ovmr(k))          ! h2o
            colamt(k,2) = max(temcol(k), coldry(k)*gasvmr_co2(iplon,k))  ! co2
            colamt(k,3) = max(temcol(k), coldry(k)*o3vmr(k))           ! o3
          enddo

!  --- ...  set up col amount for rare gases, convert from volume mixing ratio
!           to molec/cm2 based on coldry (scaled to 1.0e-20)

          if (ilwrgas > 0) then
            do k = 1, nlay
              colamt(k,4)=max(temcol(k), coldry(k)*gasvmr_n2o(iplon,k))  ! n2o
              colamt(k,5)=max(temcol(k), coldry(k)*gasvmr_ch4(iplon,k))  ! ch4
              colamt(k,6)=max(f_zero,    coldry(k)*gasvmr_o2(iplon,k))  ! o2
              colamt(k,7)=max(f_zero,    coldry(k)*gasvmr_co(iplon,k))  ! co

              wx(k,1) = max( f_zero, coldry(k)*gasvmr_ccl4(iplon,k) )   ! ccl4
              wx(k,2) = max( f_zero, coldry(k)*gasvmr_cfc11(iplon,k) )   ! cf11
              wx(k,3) = max( f_zero, coldry(k)*gasvmr_cfc12(iplon,k) )   ! cf12
              wx(k,4) = max( f_zero, coldry(k)*gasvmr_cfc22(iplon,k) )   ! cf22
            enddo
          else
            do k = 1, nlay
              colamt(k,4) = f_zero     ! n2o
              colamt(k,5) = f_zero     ! ch4
              colamt(k,6) = f_zero     ! o2
              colamt(k,7) = f_zero     ! co

              wx(k,1) = f_zero
              wx(k,2) = f_zero
              wx(k,3) = f_zero
              wx(k,4) = f_zero
            enddo
          endif

!  --- ...  set aerosol optical properties

          do j = 1, nbands
            do k = 1, nlay
              tauaer(j,k) = aeraod(iplon,k,j)                           &
     &                    * (f_one - aerssa(iplon,k,j))
            enddo
          enddo

          if (ilwcliq > 0) then    ! use prognostic cloud method
            do k = 1, nlay
              cldfrc(k)= cld_cf(iplon,k)
              clwp(k)  = cld_lwp(iplon,k)
              relw(k)  = cld_ref_liq(iplon,k)
              ciwp(k)  = cld_iwp(iplon,k)
              reiw(k)  = cld_ref_ice(iplon,k)
              cda1(k)  = cld_rwp(iplon,k)
              cda2(k)  = cld_ref_rain(iplon,k)
              cda3(k)  = cld_swp(iplon,k)
              cda4(k)  = cld_ref_snow(iplon,k)
            enddo
          else                       ! use diagnostic cloud method
            do k = 1, nlay
              cldfrc(k)= cld_cf(iplon,k)
              cda1(k)  = cld_od(iplon,k)
            enddo
          endif                      ! end if_ilwcliq

          cldfrc(0)    = f_one       ! padding value only
          cldfrc(nlp1) = f_zero      ! padding value only

!  --- ...  compute precipitable water vapor for diffusivity angle adjustments

          tem1 = f_zero
          tem2 = f_zero
          do k = 1, nlay
            tem1 = tem1 + coldry(k) + colamt(k,1)
            tem2 = tem2 + colamt(k,1)
          enddo

          tem0 = 10.0 * tem2 / (amdw * tem1 * con_g)
          pwvcm = tem0 * plvl(iplon,1)

        endif                       ! if_ivflip

!> -# Compute column amount for broadening gases.

        do k = 1, nlay
          summol = f_zero
          do i = 2, maxgas
            summol = summol + colamt(k,i)
          enddo
          colbrd(k) = coldry(k) - summol
        enddo

!> -# Compute diffusivity angle adjustments.

        tem1 = 1.80
        tem2 = 1.50
        do j = 1, nbands
          if (j==1 .or. j==4 .or. j==10) then
            secdiff(j) = 1.66
          else
            secdiff(j) = min( tem1, max( tem2,                          &
     &                   a0(j)+a1(j)*exp(a2(j)*pwvcm) ))
          endif
        enddo

!     if (lprnt) then
!      print *,'  coldry',coldry
!      print *,' wx(*,1) ',(wx(k,1),k=1,NLAY)
!      print *,' wx(*,2) ',(wx(k,2),k=1,NLAY)
!      print *,' wx(*,3) ',(wx(k,3),k=1,NLAY)
!      print *,' wx(*,4) ',(wx(k,4),k=1,NLAY)
!      print *,' iplon ',iplon
!      print *,'  pavel ',pavel
!      print *,'  delp ',delp
!      print *,'  tavel ',tavel
!      print *,'  tz ',tz
!      print *,' h2ovmr ',h2ovmr
!      print *,' o3vmr ',o3vmr
!     endif

!> -# For cloudy atmosphere, call cldprop() to set cloud optical
!!    properties.

        lcf1 = .false.
        lab_do_k0 : do k = 1, nlay
          if ( cldfrc(k) > eps ) then
            lcf1 = .true.
            exit lab_do_k0
          endif
        enddo  lab_do_k0

        if ( lcf1 ) then

          call cldprop                                                  &
!  ---  inputs:
     &     ( cldfrc,clwp,relw,ciwp,reiw,cda1,cda2,cda3,cda4,            &
     &       nlay, nlp1, ipseed(iplon),                                 &
!  ---  outputs:
     &       cldfmc, taucld                                             &
     &     )

        else
          cldfmc = f_zero
          taucld = f_zero
        endif

!     if (lprnt) then
!      print *,' after cldprop'
!      print *,' clwp',clwp
!      print *,' ciwp',ciwp
!      print *,' relw',relw
!      print *,' reiw',reiw
!      print *,' taucl',cda1
!      print *,' cldfrac',cldfrc
!     endif

!> -# Calling setcoef() to compute various coefficients needed in
!!    radiative transfer calculations.
        call setcoef                                                    &
!  ---  inputs:
     &     ( pavel,tavel,tz,stemp,h2ovmr,colamt,coldry,colbrd,          &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       laytrop,pklay,pklev,jp,jt,jt1,                             &
     &       rfrate,fac00,fac01,fac10,fac11,                            &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor                 &
     &     )

!     if (lprnt) then
!      print *,'laytrop',laytrop
!      print *,'colh2o',(colamt(k,1),k=1,NLAY)
!      print *,'colco2',(colamt(k,2),k=1,NLAY)
!      print *,'colo3', (colamt(k,3),k=1,NLAY)
!      print *,'coln2o',(colamt(k,4),k=1,NLAY)
!      print *,'colch4',(colamt(k,5),k=1,NLAY)
!      print *,'fac00',fac00
!      print *,'fac01',fac01
!      print *,'fac10',fac10
!      print *,'fac11',fac11
!      print *,'jp',jp
!      print *,'jt',jt
!      print *,'jt1',jt1
!      print *,'selffac',selffac
!      print *,'selffrac',selffrac
!      print *,'indself',indself
!      print *,'forfac',forfac
!      print *,'forfrac',forfrac
!      print *,'indfor',indfor
!     endif

!> -# Call taumol() to calculte the gaseous optical depths and Plank
!! fractions for each longwave spectral band.

        call taumol                                                     &
!  ---  inputs:
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              &
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay,                                                      &
!  ---  outputs:
     &       fracs, tautot                                              &
     &     )

!     if (lprnt) then
!     print *,' after taumol'
!     do k = 1, nlay
!       write(6,121) k
!121    format(' k =',i3,5x,'FRACS')
!       write(6,122) (fracs(j,k),j=1,ngptlw)
!122    format(10e14.7)
!       write(6,123) k
!123    format(' k =',i3,5x,'TAUTOT')
!       write(6,122) (tautot(j,k),j=1,ngptlw)
!     enddo
!     endif

!> -# Call the radiative transfer routine based on cloud scheme
!!    selection. Compute the upward/downward radiative fluxes, and
!!    heating rates for both clear or cloudy atmosphere.
!!\n  - call rtrn(): clouds are assumed as randomly overlaping in a
!!                   vertical column
!!\n  - call rtrnmr(): clouds are assumed as in maximum-randomly
!!                     overlaping in a vertical column;
!!\n  - call rtrnmc(): clouds are treated with the mcica stochastic
!!                     approach.

        if (isubclw <= 0) then

          if (iovrlw <= 0) then

            call rtrn                                                   &
!  ---  inputs:
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdiff,nlay,nlp1,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
     &     )

          else

            call rtrnmr                                                 &
!  ---  inputs:
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdiff,nlay,nlp1,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
     &     )

          endif   ! end if_iovrlw_block

        else

          call rtrnmc                                                   &
!  ---  inputs:
     &     ( semiss,delp,cldfmc,taucld,tautot,pklay,pklev,              &
     &       fracs,secdiff,nlay,nlp1,                                   &
!  ---  outputs:
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &
     &     )

        endif   ! end if_isubclw_block

!> -# Save outputs.

        topflx(iplon)%upfxc = totuflux(nlay)
        topflx(iplon)%upfx0 = totuclfl(nlay)

        sfcflx(iplon)%upfxc = totuflux(0)
        sfcflx(iplon)%upfx0 = totuclfl(0)
        sfcflx(iplon)%dnfxc = totdflux(0)
        sfcflx(iplon)%dnfx0 = totdclfl(0)

        if (ivflip == 0) then       ! output from toa to sfc

!! --- ...  optional fluxes
          if ( lflxprf ) then
            do k = 0, nlay
              k1 = nlp1 - k
              flxprf(iplon,k1)%upfxc = totuflux(k)
              flxprf(iplon,k1)%dnfxc = totdflux(k)
              flxprf(iplon,k1)%upfx0 = totuclfl(k)
              flxprf(iplon,k1)%dnfx0 = totdclfl(k)
            enddo
          endif

          do k = 1, nlay
            k1 = nlp1 - k
            hlwc(iplon,k1) = htr(k)
          enddo

!! --- ...  optional clear sky heating rate
          if ( lhlw0 ) then
            do k = 1, nlay
              k1 = nlp1 - k
              hlw0(iplon,k1) = htrcl(k)
            enddo
          endif

!! --- ...  optional spectral band heating rate
          if ( lhlwb ) then
            do j = 1, nbands
            do k = 1, nlay
              k1 = nlp1 - k
              hlwb(iplon,k1,j) = htrb(k,j)
            enddo
            enddo
          endif

        else                        ! output from sfc to toa

!! --- ...  optional fluxes
          if ( lflxprf ) then
            do k = 0, nlay
              flxprf(iplon,k+1)%upfxc = totuflux(k)
              flxprf(iplon,k+1)%dnfxc = totdflux(k)
              flxprf(iplon,k+1)%upfx0 = totuclfl(k)
              flxprf(iplon,k+1)%dnfx0 = totdclfl(k)
            enddo
          endif

          do k = 1, nlay
            hlwc(iplon,k) = htr(k)
          enddo

!! --- ...  optional clear sky heating rate
          if ( lhlw0 ) then
            do k = 1, nlay
              hlw0(iplon,k) = htrcl(k)
            enddo
          endif

!! --- ...  optional spectral band heating rate
          if ( lhlwb ) then
            do j = 1, nbands
            do k = 1, nlay
              hlwb(iplon,k,j) = htrb(k,j)
            enddo
            enddo
          endif

        endif                       ! if_ivflip

      enddo  lab_do_iplon

!...................................
      end subroutine rrtmg_lw_run
!-----------------------------------
!> @}
      subroutine rrtmg_lw_finalize ()
      end subroutine rrtmg_lw_finalize 



!> \ingroup module_radlw_main
!> \brief This subroutine performs calculations necessary for the initialization
!! of the longwave model, which includes non-varying model variables, conversion
!! factors, and look-up tables  
!!
!! Lookup tables are computed for use in the lw
!! radiative transfer, and input absorption coefficient data for each
!! spectral band are reduced from 256 g-point intervals to 140.
!!\param me        print control for parallel process
!!\section rlwinit_gen rlwinit General Algorithm
!! @{
      subroutine rlwinit                                                &
     &     ( me ) !  ---  inputs
!  ---  outputs: (none)

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  initialize non-varying module variables, conversion factors,!
! and look-up tables.                                                   !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!    me       - print control for parallel process                      !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  external module variables:  (in physparam)                            !
!   ilwrate - heating rate unit selections                              !
!           =1: output in k/day                                         !
!           =2: output in k/second                                      !
!   ilwrgas - control flag for rare gases (ch4,n2o,o2,cfcs, etc.)       !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   ilwcliq - liquid cloud optical properties contrl flag               !
!           =0: input cloud opt depth from diagnostic scheme            !
!           >0: input cwp,rew, and other cloud content parameters       !
!   isubclw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   icldflg - cloud scheme control flag                                 !
!           =0: diagnostic scheme gives cloud tau, omiga, and g.        !
!           =1: prognostic scheme gives cloud liq/ice path, etc.        !
!   iovrlw  - clouds vertical overlapping control flag                  !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud (isubcol>0 only)                  !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:       michael j. iacono; july, 1998                !
!  first revision for ncar ccm:               september, 1998           !
!  second revision for rrtm_v3.0:             september, 2002           !
!                                                                       !
!  this subroutine performs calculations necessary for the initialization
!  of the longwave model.  lookup tables are computed for use in the lw !
!  radiative transfer, and input absorption coefficient data for each   !
!  spectral band are reduced from 256 g-point intervals to 140.         !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
! definitions:                                                          !
!   arrays for 10000-point look-up tables:                              !
!   tau_tbl - clear-sky optical depth (used in cloudy radiative transfer!
!   exp_tbl - exponential lookup table for tansmittance                 !
!   tfn_tbl - tau transition function; i.e. the transition of the Planck!
!             function from that for the mean layer temperature to that !
!             for the layer boundary temperature as a function of optical
!             depth. the "linear in tau" method is used to make the table
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: me

!  ---  outputs: none

!  ---  locals:
      real (kind=kind_phys), parameter :: expeps = 1.e-20

      real (kind=kind_phys) :: tfn, pival, explimit

      integer               :: i

!
!===> ... begin here
!
      if ( iovrlw<0 .or. iovrlw>2 ) then
        print *,'  *** Error in specification of cloud overlap flag',   &
     &          ' IOVRLW=',iovrlw,' in RLWINIT !!'
        stop
      elseif ( iovrlw==2 .and. isubclw==0 ) then
        if (me == 0) then
          print *,'  *** IOVRLW=2 - maximum cloud overlap, is not yet', &
     &          ' available for ISUBCLW=0 setting!!'
          print *,'      The program uses maximum/random overlap',      &
     &          ' instead.'
        endif

        iovrlw = 1
      endif

      if (me == 0) then
        print *,' - Using AER Longwave Radiation, Version: ', VTAGLW

        if (ilwrgas > 0) then
          print *,'   --- Include rare gases N2O, CH4, O2, CFCs ',      &
     &            'absorptions in LW'
        else
          print *,'   --- Rare gases effect is NOT included in LW'
        endif

        if ( isubclw == 0 ) then
          print *,'   --- Using standard grid average clouds, no ',     &
     &            'sub-column clouds approximation applied'
        elseif ( isubclw == 1 ) then
          print *,'   --- Using MCICA sub-colum clouds approximation ', &
     &            'with a prescribed sequence of permutaion seeds'
        elseif ( isubclw == 2 ) then
          print *,'   --- Using MCICA sub-colum clouds approximation ', &
     &            'with provided input array of permutation seeds'
        else
          print *,'  *** Error in specification of sub-column cloud ',  &
     &            ' control flag isubclw =',isubclw,' !!'
          stop
        endif
      endif

!> -# Check cloud flags for consistency.

      if ((icldflg == 0 .and. ilwcliq /= 0) .or.                        &
     &    (icldflg == 1 .and. ilwcliq == 0)) then
        print *,'  *** Model cloud scheme inconsistent with LW',        &
     &          ' radiation cloud radiative property setup !!'
        stop
      endif

!> -# Setup default surface emissivity for each band.

      semiss0(:) = f_one

!> -# Setup constant factors for flux and heating rate
!! the 1.0e-2 is to convert pressure from mb to \f$N/m^2\f$.

      pival = 2.0 * asin(f_one)
      fluxfac = pival * 2.0d4
!     fluxfac = 62831.85307179586                   ! = 2 * pi * 1.0e4

      if (ilwrate == 1) then
!       heatfac = 8.4391
!       heatfac = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
        heatfac = con_g * 864.0 / con_cp            !   (in k/day)
      else
        heatfac = con_g * 1.0e-2 / con_cp           !   (in k/second)
      endif

!> -# Compute lookup tables for transmittance, tau transition
!! function, and clear sky tau (for the cloudy sky radiative
!! transfer).  tau is computed as a function of the tau
!! transition function, transmittance is calculated as a
!! function of tau, and the tau transition function is
!! calculated using the linear in tau formulation at values of
!! tau above 0.01.  tf is approximated as tau/6 for tau < 0.01.
!! all tables are computed at intervals of 0.001.  the inverse
!! of the constant used in the pade approximation to the tau
!! transition function is set to b.

      tau_tbl(0) = f_zero
      exp_tbl(0) = f_one
      tfn_tbl(0) = f_zero

      tau_tbl(ntbl) = 1.e10
      exp_tbl(ntbl) = expeps
      tfn_tbl(ntbl) = f_one

      explimit = aint( -log(tiny(exp_tbl(0))) )

      do i = 1, ntbl-1
!org    tfn = float(i) / float(ntbl)
!org    tau_tbl(i) = bpade * tfn / (f_one - tfn)
        tfn = real(i, kind_phys) / real(ntbl-i, kind_phys)
        tau_tbl(i) = bpade * tfn
        if (tau_tbl(i) >= explimit) then
          exp_tbl(i) = expeps
        else
          exp_tbl(i) = exp( -tau_tbl(i) )
        endif

        if (tau_tbl(i) < 0.06) then
          tfn_tbl(i) = tau_tbl(i) / 6.0
        else
          tfn_tbl(i) = f_one - 2.0*( (f_one / tau_tbl(i))               &
     &               - ( exp_tbl(i) / (f_one - exp_tbl(i)) ) )
        endif
      enddo

!...................................
      end subroutine rlwinit
!! @}
!-----------------------------------


!>\ingroup module_radlw_main
!> \brief This subroutine computes the cloud optical depth(s) for each cloudy
!! layer and g-point interval.
!!\param cfrac           layer cloud fraction
!!\n     ---  for  ilwcliq > 0 (prognostic cloud scheme)  - - -
!!\param cliqp           layer in-cloud liq water path (\f$g/m^2\f$)
!!\param reliq           mean eff radius for liq cloud (micron)
!!\param cicep           layer in-cloud ice water path (\f$g/m^2\f$)
!!\param reice           mean eff radius for ice cloud (micron)
!!\param cdat1           layer rain drop water path (\f$g/m^2\f$)
!!\param cdat2           effective radius for rain drop (micron)
!!\param cdat3           layer snow flake water path(\f$g/m^2\f$)
!!\param cdat4           mean effective radius for snow flake(micron)
!!\param cliqp           not used
!!\param cicep           not used
!!\param reliq           not used
!!\param reice           not used
!!\param cdat1           layer cloud optical depth
!!\param cdat2           layer cloud single scattering albedo
!!\param cdat3           layer cloud asymmetry factor
!!\param cdat4           optional use
!!\param nlay            number of layer number
!!\param nlp1            number of veritcal levels
!!\param ipseed          permutation seed for generating random numbers (isubclw>0)
!!\param cldfmc          cloud fraction for each sub-column
!!\param taucld          cloud optical depth for bands (non-mcica)
!!\section gen_cldprop cldprop General Algorithm
!> @{
      subroutine cldprop                                                &
     &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     & !  ---  inputs
     &       nlay, nlp1, ipseed,                                        &
     &       cldfmc, taucld                                             & !  ---  outputs
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the cloud optical depth(s) for each cloudy layer    !
! and g-point interval.                                                 !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       -size- !
!    cfrac - real, layer cloud fraction                          0:nlp1 !
!        .....  for ilwcliq > 0  (prognostic cloud sckeme)  - - -       !
!    cliqp - real, layer in-cloud liq water path (g/m**2)          nlay !
!    reliq - real, mean eff radius for liq cloud (micron)          nlay !
!    cicep - real, layer in-cloud ice water path (g/m**2)          nlay !
!    reice - real, mean eff radius for ice cloud (micron)          nlay !
!    cdat1 - real, layer rain drop water path  (g/m**2)            nlay !
!    cdat2 - real, effective radius for rain drop (microm)         nlay !
!    cdat3 - real, layer snow flake water path (g/m**2)            nlay !
!    cdat4 - real, effective radius for snow flakes (micron)       nlay !
!        .....  for ilwcliq = 0  (diagnostic cloud sckeme)  - - -       !
!    cdat1 - real, input cloud optical depth                       nlay !
!    cdat2 - real, layer cloud single scattering albedo            nlay !
!    cdat3 - real, layer cloud asymmetry factor                    nlay !
!    cdat4 - real, optional use                                    nlay !
!    cliqp - not used                                              nlay !
!    reliq - not used                                              nlay !
!    cicep - not used                                              nlay !
!    reice - not used                                              nlay !
!                                                                       !
!    nlay  - integer, number of vertical layers                      1  !
!    nlp1  - integer, number of vertical levels                      1  !
!    ipseed- permutation seed for generating random numbers (isubclw>0) !
!                                                                       !
!  outputs:                                                             !
!    cldfmc - real, cloud fraction for each sub-column       ngptlw*nlay!
!    taucld - real, cld opt depth for bands (non-mcica)      nbands*nlay!
!                                                                       !
!  explanation of the method for each value of ilwcliq, and ilwcice.    !
!    set up in module "module_radlw_cntr_para"                          !
!                                                                       !
!     ilwcliq=0  : input cloud optical property (tau, ssa, asy).        !
!                  (used for diagnostic cloud method)                   !
!     ilwcliq>0  : input cloud liq/ice path and effective radius, also  !
!                  require the user of 'ilwcice' to specify the method  !
!                  used to compute aborption due to water/ice parts.    !
!  ...................................................................  !
!                                                                       !
!     ilwcliq=1:   the water droplet effective radius (microns) is input!
!                  and the opt depths due to water clouds are computed  !
!                  as in hu and stamnes, j., clim., 6, 728-742, (1993). !
!                  the values for absorption coefficients appropriate for
!                  the spectral bands in rrtm have been obtained for a  !
!                  range of effective radii by an averaging procedure   !
!                  based on the work of j. pinto (private communication).
!                  linear interpolation is used to get the absorption   !
!                  coefficients for the input effective radius.         !
!                                                                       !
!     ilwcice=1:   the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in ebert and curry, jgr, 97,  !
!                  3831-3836 (1992).  the spectral regions in this work !
!                  have been matched with the spectral bands in rrtm to !
!                  as great an extent as possible:                      !
!                     e&c 1      ib = 5      rrtm bands 9-16            !
!                     e&c 2      ib = 4      rrtm bands 6-8             !
!                     e&c 3      ib = 3      rrtm bands 3-5             !
!                     e&c 4      ib = 2      rrtm band 2                !
!                     e&c 5      ib = 1      rrtm band 1                !
!     ilwcice=2:   the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are input and the optical depths due to ice!
!                  clouds are computed as in rt code, streamer v3.0     !
!                  (ref: key j., streamer user's guide, cooperative     !
!                  institute for meteorological satellite studies, 2001,!
!                  96 pp.) valid range of values for re are between 5.0 !
!                  and 131.0 micron.                                    !
!     ilwcice=3:   the ice generalized effective size (dge) is input and!
!                  the optical properties, are calculated as in q. fu,  !
!                  j. climate, (1998). q. fu provided high resolution   !
!                  tales which were appropriately averaged for the bands!
!                  in rrtm_lw. linear interpolation is used to get the  !
!                  coeff from the stored tables. valid range of values  !
!                  for deg are between 5.0 and 140.0 micron.            !
!                                                                       !
!  other cloud control module variables:                                !
!     isubclw =0: standard cloud scheme, no sub-col cloud approximation !
!             >0: mcica sub-col cloud scheme using ipseed as permutation!
!                 seed for generating rundom numbers                    !
!                                                                       !
!  ======================  end of description block  =================  !
!
      use module_radlw_cldprlw

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1, ipseed

      real (kind=kind_phys), dimension(0:nlp1), intent(in) :: cfrac
      real (kind=kind_phys), dimension(nlay),   intent(in) :: cliqp,    &
     &       reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4

!  ---  outputs:
      real (kind=kind_phys), dimension(ngptlw,nlay),intent(out):: cldfmc
      real (kind=kind_phys), dimension(nbands,nlay),intent(out):: taucld

!  ---  locals:
      real (kind=kind_phys), dimension(nbands) :: tauliq, tauice
      real (kind=kind_phys), dimension(nlay)   :: cldf

      real (kind=kind_phys) :: dgeice, factor, fint, tauran, tausnw,    &
     &       cldliq, refliq, cldice, refice

      logical :: lcloudy(ngptlw,nlay)
      integer :: ia, ib, ig, k, index

!
!===> ...  begin here
!
      do k = 1, nlay
        do ib = 1, nbands
          taucld(ib,k) = f_zero
        enddo
      enddo

      do k = 1, nlay
        do ig = 1, ngptlw
          cldfmc(ig,k) = f_zero
        enddo
      enddo

!> -# Compute cloud radiative properties for a cloudy column:
!!\n  - Compute cloud radiative properties for rain and snow (tauran,tausnw)
!!\n  - Calculation of absorption coefficients due to water clouds(tauliq)
!!\n  - Calculation of absorption coefficients due to ice clouds (tauice).
!!\n  - For prognostic cloud scheme: sum up the cloud optical property:
!!\n    \f$ taucld=tauice+tauliq+tauran+tausnw \f$

!  --- ...  compute cloud radiative properties for a cloudy column

      lab_if_ilwcliq : if (ilwcliq > 0) then

        lab_do_k : do k = 1, nlay
          lab_if_cld : if (cfrac(k) > cldmin) then

            tauran = absrain * cdat1(k)                      ! ncar formula
!!          tausnw = abssnow1 * cdat3(k)                     ! ncar formula
!  ---  if use fu's formula it needs to be normalized by snow density
!       !not use snow density = 0.1 g/cm**3 = 0.1 g/(mu * m**2)
!       use ice density = 0.9167 g/cm**3 = 0.9167 g/(mu * m**2)
!       factor 1.5396=8/(3*sqrt(3)) converts reff to generalized ice particle size
!       use newer factor value 1.0315
!       1/(0.9167*1.0315) = 1.05756
            if (cdat3(k)>f_zero .and. cdat4(k)>10.0_kind_phys) then
              tausnw = abssnow0*1.05756*cdat3(k)/cdat4(k)      ! fu's formula
            else
              tausnw = f_zero
            endif

            cldliq = cliqp(k)
            cldice = cicep(k)
!           refliq = max(2.5e0, min(60.0e0, reliq(k) ))
!           refice = max(5.0e0, reice(k) )
            refliq = reliq(k)
            refice = reice(k)

!  --- ...  calculation of absorption coefficients due to water clouds.

            if ( cldliq <= f_zero ) then
              do ib = 1, nbands
                tauliq(ib) = f_zero
              enddo
            else
              if ( ilwcliq == 1 ) then

                factor = refliq - 1.5
                index  = max( 1, min( 57, int( factor ) ))
                fint   = factor - float(index)

                do ib = 1, nbands
                  tauliq(ib) = max(f_zero, cldliq*(absliq1(index,ib)    &
     &              + fint*(absliq1(index+1,ib)-absliq1(index,ib)) ))
                enddo
              endif   ! end if_ilwcliq_block
            endif   ! end if_cldliq_block

!  --- ...  calculation of absorption coefficients due to ice clouds.

            if ( cldice <= f_zero ) then
              do ib = 1, nbands
                tauice(ib) = f_zero
              enddo
            else

!  --- ...  ebert and curry approach for all particle sizes though somewhat
!           unjustified for large ice particles

              if ( ilwcice == 1 ) then
                refice = min(130.0, max(13.0, real(refice) ))

                do ib = 1, nbands
                  ia = ipat(ib)             ! eb_&_c band index for ice cloud coeff
                  tauice(ib) = max(f_zero, cldice*(absice1(1,ia)        &
     &                         + absice1(2,ia)/refice) )
                enddo

!  --- ...  streamer approach for ice effective radius between 5.0 and 131.0 microns
!           and ebert and curry approach for ice eff radius greater than 131.0 microns.
!           no smoothing between the transition of the two methods.

              elseif ( ilwcice == 2 ) then

                factor = (refice - 2.0) / 3.0
                index  = max( 1, min( 42, int( factor ) ))
                fint   = factor - float(index)

                do ib = 1, nbands
                  tauice(ib) = max(f_zero, cldice*(absice2(index,ib)    &
     &              + fint*(absice2(index+1,ib) - absice2(index,ib)) ))
                enddo

!  --- ...  fu's approach for ice effective radius between 4.8 and 135 microns
!           (generalized effective size from 5 to 140 microns)

              elseif ( ilwcice == 3 ) then

!               dgeice = max(5.0, 1.5396*refice)              ! v4.4 value
                dgeice = max(5.0, 1.0315*refice)              ! v4.71 value
                factor = (dgeice - 2.0) / 3.0
                index  = max( 1, min( 45, int( factor ) ))
                fint   = factor - float(index)

                do ib = 1, nbands
                  tauice(ib) = max(f_zero, cldice*(absice3(index,ib)    &
     &              + fint*(absice3(index+1,ib) - absice3(index,ib)) ))
                enddo

              endif   ! end if_ilwcice_block
            endif   ! end if_cldice_block

            do ib = 1, nbands
              taucld(ib,k) = tauice(ib) + tauliq(ib) + tauran + tausnw
            enddo

          endif  lab_if_cld
        enddo  lab_do_k

      else  lab_if_ilwcliq

        do k = 1, nlay
          if (cfrac(k) > cldmin) then
            do ib = 1, nbands
              taucld(ib,k) = cdat1(k)
            enddo
          endif
        enddo

      endif  lab_if_ilwcliq

!> -# if physparam::isubclw > 0, call mcica_subcol() to distribute
!!    cloud properties to each g-point.

      if ( isubclw > 0 ) then      ! mcica sub-col clouds approx
        do k = 1, nlay
          if ( cfrac(k) < cldmin ) then
            cldf(k) = f_zero
          else
            cldf(k) = cfrac(k)
          endif
        enddo

!  --- ...  call sub-column cloud generator

        call mcica_subcol                                               &
!  ---  inputs:
     &     ( cldf, nlay, ipseed,                                        &
!  ---  output:
     &       lcloudy                                                    &
     &     )

        do k = 1, nlay
          do ig = 1, ngptlw
            if ( lcloudy(ig,k) ) then
              cldfmc(ig,k) = f_one
            else
              cldfmc(ig,k) = f_zero
            endif
          enddo
        enddo

      endif   ! end if_isubclw_block

      return
! ..................................
      end subroutine cldprop
! ----------------------------------
!> @}

!>\ingroup module_radlw_main
!>\brief This suroutine computes sub-colum cloud profile flag array.
!!\param cldf        layer cloud fraction
!!\param nlay        number of model vertical layers
!!\param ipseed      permute seed for random num generator
!!\param lcloudy     sub-colum cloud profile flag array
!!\section mcica_subcol_gen mcica_subcol General Algorithm
!! @{
      subroutine mcica_subcol                                           &
     &    ( cldf, nlay, ipseed,                                         &!  ---  inputs
     &      lcloudy                                                     & !  ---  outputs
     &    )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                size !
!   cldf    - real, layer cloud fraction                           nlay !
!   nlay    - integer, number of model vertical layers               1  !
!   ipseed  - integer, permute seed for random num generator         1  !
!    ** note : if the cloud generator is called multiple times, need    !
!              to permute the seed between each call; if between calls  !
!              for lw and sw, use values differ by the number of g-pts. !
!                                                                       !
!  output variables:                                                    !
!   lcloudy - logical, sub-colum cloud profile flag array    ngptlw*nlay!
!                                                                       !
!  other control flags from module variables:                           !
!     iovrlw    : control flag for cloud overlapping method             !
!                 =0:random; =1:maximum/random: =2:maximum              !
!                                                                       !
!  =====================    end of definitions    ====================  !

      implicit none

!  ---  inputs:
      integer, intent(in) :: nlay, ipseed

      real (kind=kind_phys), dimension(nlay), intent(in) :: cldf

!  ---  outputs:
      logical, dimension(ngptlw,nlay), intent(out) :: lcloudy

!  ---  locals:
      real (kind=kind_phys) :: cdfunc(ngptlw,nlay), rand1d(ngptlw),     &
     &       rand2d(nlay*ngptlw), tem1

      type (random_stat) :: stat          ! for thread safe random generator

      integer :: k, n, k1
!
!===> ...  begin here
!
!> -# Call random_setseed() to advance randum number generator by ipseed values.

      call random_setseed                                               &
!  ---  inputs:
     &    ( ipseed,                                                     &
!  ---  outputs:
     &      stat                                                        &
     &    )

!> -# Sub-column set up according to overlapping assumption:
!!  - For random overlap, pick a random value at every level 
!!  - For max-random overlap, pick a random value at every level
!!  - For maximum overlap, pick same random numebr at every level

      select case ( iovrlw )

        case( 0 )        ! random overlap, pick a random value at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptlw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(n,k) = rand2d(k1)
            enddo
          enddo

        case( 1 )        ! max-ran overlap

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptlw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(n,k) = rand2d(k1)
            enddo
          enddo

!  ---  first pick a random number for bottom (or top) layer.
!       then walk up the column: (aer's code)
!       if layer below is cloudy, use the same rand num in the layer below
!       if layer below is clear,  use a new random number

!  ---  from bottom up
          do k = 2, nlay
            k1 = k - 1
            tem1 = f_one - cldf(k1)

            do n = 1, ngptlw
              if ( cdfunc(n,k1) > tem1 ) then
                cdfunc(n,k) = cdfunc(n,k1)
              else
                cdfunc(n,k) = cdfunc(n,k) * tem1
              endif
            enddo
          enddo

!  ---  or walk down the column: (if use original author's method)
!       if layer above is cloudy, use the same rand num in the layer above
!       if layer above is clear,  use a new random number

!  ---  from top down
!         do k = nlay-1, 1, -1
!           k1 = k + 1
!           tem1 = f_one - cldf(k1)

!           do n = 1, ngptlw
!             if ( cdfunc(n,k1) > tem1 ) then
!               cdfunc(n,k) = cdfunc(n,k1)
!             else
!               cdfunc(n,k) = cdfunc(n,k) * tem1
!             endif
!           enddo
!         enddo

        case( 2 )        !<  - For maximum overlap, pick same random numebr at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand1d, stat )

          do n = 1, ngptlw
            tem1 = rand1d(n)

            do k = 1, nlay
              cdfunc(n,k) = tem1
            enddo
          enddo

      end select

!> -# Generate subcolumns for homogeneous clouds.

      do k = 1, nlay
        tem1 = f_one - cldf(k)

        do n = 1, ngptlw
          lcloudy(n,k) = cdfunc(n,k) >= tem1
        enddo
      enddo

      return
! ..................................
      end subroutine mcica_subcol
!! @}
! ----------------------------------

!>\ingroup module_radlw_main
!> This subroutine computes various coefficients needed in radiative
!! transfer calculations.
!!\param pavel           layer pressure (mb)
!!\param tavel           layer temperature (K)
!!\param tz              level(interface) temperatures (K)
!!\param stemp           surface ground temperature (K)
!!\param h2ovmr          layer w.v. volumn mixing ratio (kg/kg)
!!\param colamt           column amounts of absorbing gases.
!! 2nd indices range: 1-maxgas, for watervapor,carbon dioxide, ozone,
!! nitrous oxide, methane,oxigen, carbon monoxide,etc. \f$(mol/cm^2)\f$
!!\param coldry          dry air column amount
!!\param colbrd          column amount of broadening gases
!!\param nlay            total number of vertical layers
!!\param nlp1            total number of vertical levels
!!\param laytrop         tropopause layer index (unitless)
!!\param pklay           integrated planck func at lay temp
!!\param pklev           integrated planck func at lev temp
!!\param jp              indices of lower reference pressure
!!\param jt, jt1         indices of lower reference temperatures
!!\param rfrate          ref ratios of binary species param
!!\n                     (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,
!!                                4-h2o/ch4,5-n2o/co2,6-o3/co2
!!\n                     (:,:,n)n=1,2: the rates of ref press at
!!                                the 2 sides of the layer
!!\param facij           factors multiply the reference ks, i,j=0/1 for
!!                       lower/higher of the 2 appropriate temperatures
!!                       and altitudes.
!!\param selffac         scale factor for w. v. self-continuum equals
!!                       (w. v. density)/(atmospheric density at 296k and 1013 mb)
!!\param selffrac        factor for temperature interpolation of
!!                       reference w. v. self-continuum data
!!\param indself         index of lower ref temp for selffac
!!\param forfac          scale factor for w. v. foreign-continuum
!!\param forfrac         factor for temperature interpolation of
!!                       reference w.v. foreign-continuum data
!!\param indfor          index of lower ref temp for forfac
!!\param minorfrac       factor for minor gases
!!\param scaleminor,scaleminorn2         scale factors for minor gases
!!\param indminor        index of lower ref temp for minor gases
!>\section setcoef_gen setcoef General Algorithm
!! 
! ----------------------------------
      subroutine setcoef                                                &
     &     ( pavel,tavel,tz,stemp,h2ovmr,colamt,coldry,colbrd,          & !  ---  inputs:
     &       nlay, nlp1,                                                &
     &       laytrop,pklay,pklev,jp,jt,jt1,                             & !  ---  outputs:
     &       rfrate,fac00,fac01,fac10,fac11,                            &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor                 &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute various coefficients needed in radiative transfer   !
!    calculations.                                                      !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       -size- !
!   pavel     - real, layer pressures (mb)                         nlay !
!   tavel     - real, layer temperatures (k)                       nlay !
!   tz        - real, level (interface) temperatures (k)         0:nlay !
!   stemp     - real, surface ground temperature (k)                1   !
!   h2ovmr    - real, layer w.v. volum mixing ratio (kg/kg)        nlay !
!   colamt    - real, column amounts of absorbing gases      nlay*maxgas!
!                 2nd indices range: 1-maxgas, for watervapor,          !
!                 carbon dioxide, ozone, nitrous oxide, methane,        !
!                 oxigen, carbon monoxide,etc. (molecules/cm**2)        !
!   coldry    - real, dry air column amount                        nlay !
!   colbrd    - real, column amount of broadening gases            nlay !
!   nlay/nlp1 - integer, total number of vertical layers, levels    1   !
!                                                                       !
!  outputs:                                                             !
!   laytrop   - integer, tropopause layer index (unitless)          1   !
!   pklay     - real, integrated planck func at lay temp   nbands*0:nlay!
!   pklev     - real, integrated planck func at lev temp   nbands*0:nlay!
!   jp        - real, indices of lower reference pressure          nlay !
!   jt, jt1   - real, indices of lower reference temperatures      nlay !
!   rfrate    - real, ref ratios of binary species param   nlay*nrates*2!
!     (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,4-h2o/ch4,5-n2o/co2,6-o3/co2!
!     (:,:,n)n=1,2: the rates of ref press at the 2 sides of the layer  !
!   facij     - real, factors multiply the reference ks,           nlay !
!                 i,j=0/1 for lower/higher of the 2 appropriate         !
!                 temperatures and altitudes.                           !
!   selffac   - real, scale factor for w. v. self-continuum        nlay !
!                 equals (w. v. density)/(atmospheric density           !
!                 at 296k and 1013 mb)                                  !
!   selffrac  - real, factor for temperature interpolation of      nlay !
!                 reference w. v. self-continuum data                   !
!   indself   - integer, index of lower ref temp for selffac       nlay !
!   forfac    - real, scale factor for w. v. foreign-continuum     nlay !
!   forfrac   - real, factor for temperature interpolation of      nlay !
!                 reference w.v. foreign-continuum data                 !
!   indfor    - integer, index of lower ref temp for forfac        nlay !
!   minorfrac - real, factor for minor gases                       nlay !
!   scaleminor,scaleminorn2                                             !
!             - real, scale factors for minor gases                nlay !
!   indminor  - integer, index of lower ref temp for minor gases   nlay !
!                                                                       !
!  ======================    end of definitions    ===================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nlay,maxgas),intent(in):: colamt
      real (kind=kind_phys), dimension(0:nlay),     intent(in):: tz

      real (kind=kind_phys), dimension(nlay), intent(in) :: pavel,      &
     &       tavel, h2ovmr, coldry, colbrd

      real (kind=kind_phys), intent(in) :: stemp

!  ---  outputs:
      integer, dimension(nlay), intent(out) :: jp, jt, jt1, indself,    &
     &       indfor, indminor

      integer, intent(out) :: laytrop

      real (kind=kind_phys), dimension(nlay,nrates,2), intent(out) ::   &
     &       rfrate
      real (kind=kind_phys), dimension(nbands,0:nlay), intent(out) ::   &
     &       pklev, pklay

      real (kind=kind_phys), dimension(nlay),          intent(out) ::   &
     &       fac00, fac01, fac10, fac11, selffac, selffrac, forfac,     &
     &       forfrac, minorfrac, scaleminor, scaleminorn2

!  ---  locals:
      real (kind=kind_phys) :: tlvlfr, tlyrfr, plog, fp, ft, ft1,       &
     &       tem1, tem2

      integer :: i, k, jp1, indlev, indlay
!
!===> ... begin here
!
!> -# Calculate information needed by the radiative transfer routine
!! that is specific to this atmosphere, especially some of the
!! coefficients and indices needed to compute the optical depths
!! by interpolating data from stored reference atmospheres.

      indlay = min(180, max(1, int(stemp-159.0) ))
      indlev = min(180, max(1, int(tz(0)-159.0) ))
      tlyrfr = stemp - int(stemp)
      tlvlfr = tz(0) - int(tz(0))
      do i = 1, nbands
        tem1 = totplnk(indlay+1,i) - totplnk(indlay,i)
        tem2 = totplnk(indlev+1,i) - totplnk(indlev,i)
        pklay(i,0) = delwave(i) * (totplnk(indlay,i) + tlyrfr*tem1)
        pklev(i,0) = delwave(i) * (totplnk(indlev,i) + tlvlfr*tem2)
      enddo

!  --- ...  begin layer loop
!> -# Calculate the integrated Planck functions for each band at the
!! surface, level, and layer temperatures.

      laytrop = 0

      do k = 1, nlay

        indlay = min(180, max(1, int(tavel(k)-159.0) ))
        tlyrfr = tavel(k) - int(tavel(k))

        indlev = min(180, max(1, int(tz(k)-159.0) ))
        tlvlfr = tz(k) - int(tz(k))

!  --- ...  begin spectral band loop

        do i = 1, nbands
          pklay(i,k) = delwave(i) * (totplnk(indlay,i) + tlyrfr         &
     &               * (totplnk(indlay+1,i) - totplnk(indlay,i)) )
          pklev(i,k) = delwave(i) * (totplnk(indlev,i) + tlvlfr         &
     &               * (totplnk(indlev+1,i) - totplnk(indlev,i)) )
        enddo

!> -# Find the two reference pressures on either side of the
!! layer pressure. store them in jp and jp1. store in fp the
!! fraction of the difference (in ln(pressure)) between these
!! two values that the layer pressure lies.

        plog = log(pavel(k))
        jp(k)= max(1, min(58, int(36.0 - 5.0*(plog+0.04)) ))
        jp1  = jp(k) + 1
!  --- ...  limit pressure extrapolation at the top
        fp   = max(f_zero, min(f_one, 5.0*(preflog(jp(k))-plog) ))
!org    fp   = 5.0 * (preflog(jp(k)) - plog)

!> -# Determine, for each reference pressure (jp and jp1), which
!! reference temperature (these are different for each
!! reference pressure) is nearest the layer temperature but does
!! not exceed it. store these indices in jt and jt1, resp.
!! store in ft (resp. ft1) the fraction of the way between jt
!! (jt1) and the next highest reference temperature that the
!! layer temperature falls.

        tem1 = (tavel(k)-tref(jp(k))) / 15.0
        tem2 = (tavel(k)-tref(jp1  )) / 15.0
        jt (k) = max(1, min(4, int(3.0 + tem1) ))
        jt1(k) = max(1, min(4, int(3.0 + tem2) ))
!  --- ...  restrict extrapolation ranges by limiting abs(det t) < 37.5 deg
        ft  = max(-0.5, min(1.5, tem1 - float(jt (k) - 3) ))
        ft1 = max(-0.5, min(1.5, tem2 - float(jt1(k) - 3) ))
!org    ft  = tem1 - float(jt (k) - 3)
!org    ft1 = tem2 - float(jt1(k) - 3)

!> -# We have now isolated the layer ln pressure and temperature,
!! between two reference pressures and two reference temperatures
!!(for each reference pressure).  we multiply the pressure
!! fraction fp with the appropriate temperature fractions to get
!! the factors that will be needed for the interpolation that yields
!! the optical depths (performed in routines taugbn for band n).

        tem1 = f_one - fp
        fac10(k) = tem1 * ft
        fac00(k) = tem1 * (f_one - ft)
        fac11(k) = fp * ft1
        fac01(k) = fp * (f_one - ft1)

        forfac(k) = pavel(k)*stpfac / (tavel(k)*(1.0 + h2ovmr(k)))
        selffac(k) = h2ovmr(k) * forfac(k)

!> -# Set up factors needed to separately include the minor gases
!! in the calculation of absorption coefficient.

        scaleminor(k) = pavel(k) / tavel(k)
        scaleminorn2(k) = (pavel(k) / tavel(k))                         &
     &                  * (colbrd(k)/(coldry(k) + colamt(k,1)))
        tem1 = (tavel(k) - 180.8) / 7.2
        indminor(k) = min(18, max(1, int(tem1)))
        minorfrac(k) = tem1 - float(indminor(k))

!> -# If the pressure is less than ~100mb, perform a different
!! set of species interpolations.

        if (plog > 4.56) then

          laytrop =  laytrop + 1

          tem1 = (332.0 - tavel(k)) / 36.0
          indfor(k) = min(2, max(1, int(tem1)))
          forfrac(k) = tem1 - float(indfor(k))

!> -# Set up factors needed to separately include the water vapor
!! self-continuum in the calculation of absorption coefficient.

          tem1 = (tavel(k) - 188.0) / 7.2
          indself(k) = min(9, max(1, int(tem1)-7))
          selffrac(k) = tem1 - float(indself(k) + 7)

!> -# Setup reference ratio to be used in calculation of binary
!! species parameter in lower atmosphere.

          rfrate(k,1,1) = chi_mls(1,jp(k)) / chi_mls(2,jp(k))
          rfrate(k,1,2) = chi_mls(1,jp(k)+1) / chi_mls(2,jp(k)+1)

          rfrate(k,2,1) = chi_mls(1,jp(k)) / chi_mls(3,jp(k))
          rfrate(k,2,2) = chi_mls(1,jp(k)+1) / chi_mls(3,jp(k)+1)

          rfrate(k,3,1) = chi_mls(1,jp(k)) / chi_mls(4,jp(k))
          rfrate(k,3,2) = chi_mls(1,jp(k)+1) / chi_mls(4,jp(k)+1)

          rfrate(k,4,1) = chi_mls(1,jp(k)) / chi_mls(6,jp(k))
          rfrate(k,4,2) = chi_mls(1,jp(k)+1) / chi_mls(6,jp(k)+1)

          rfrate(k,5,1) = chi_mls(4,jp(k)) / chi_mls(2,jp(k))
          rfrate(k,5,2) = chi_mls(4,jp(k)+1) / chi_mls(2,jp(k)+1)

        else

          tem1 = (tavel(k) - 188.0) / 36.0
          indfor(k) = 3
          forfrac(k) = tem1 - f_one

          indself(k) = 0
          selffrac(k) = f_zero

!> -# Setup reference ratio to be used in calculation of binary
!! species parameter in upper atmosphere.

          rfrate(k,1,1) = chi_mls(1,jp(k)) / chi_mls(2,jp(k))
          rfrate(k,1,2) = chi_mls(1,jp(k)+1) / chi_mls(2,jp(k)+1)

          rfrate(k,6,1) = chi_mls(3,jp(k)) / chi_mls(2,jp(k))
          rfrate(k,6,2) = chi_mls(3,jp(k)+1) / chi_mls(2,jp(k)+1)

        endif

!> -# Rescale \a selffac and \a forfac for use in taumol.

        selffac(k) = colamt(k,1) * selffac(k)
        forfac(k)  = colamt(k,1) * forfac(k)

      enddo   ! end do_k layer loop

      return
! ..................................
      end subroutine setcoef
! ----------------------------------

!>\ingroup module_radlw_main
!> This subroutine computes the upward/downward radiative fluxes, and
!! heating rates for both clear or cloudy atmosphere. Clouds assumed as
!! randomly overlaping in a vertical column.
!!\brief Original Code Description: this program calculates the upward
!! fluxes, downward fluxes, and heating rates for an arbitrary clear or
!! cloudy atmosphere. The input to this program is the atmospheric
!! profile, all Planck function information, and the cloud fraction by
!! layer. A variable diffusivity angle (secdif) is used for the angle
!! integration. Bands 2-3 and 5-9 use a value for secdif that varies
!! from 1.50 to 1.80 as a function of the column water vapor, and other
!! bands use a value of 1.66. The gaussian weight appropriate to this
!! angle (wtdiff =0.5) is applied here. Note that use of the emissivity
!! angle for the flux integration can cause errors of 1 to 4 \f$W/m^2\f$
!! within cloudy layers. Clouds are treated with a random cloud overlap
!! method.
!!\param semiss      lw surface emissivity
!!\param delp        layer pressure thickness (mb)
!!\param cldfrc      layer cloud fraction
!!\param taucld      layer cloud opt depth
!!\param tautot      total optical depth (gas+aerosols)
!!\param pklay       integrated planck function at lay temp
!!\param pklev       integrated planck func at lev temp
!!\param fracs       planck fractions
!!\param secdif      secant of diffusivity angle
!!\param nlay        number of vertical layers
!!\param nlp1        number of vertical levels (interfaces)
!!\param totuflux    total sky upward flux \f$(w/m^2)\f$
!!\param totdflux    total sky downward flux \f$(w/m^2)\f$
!!\param htr         total sky heating rate (k/sec or k/day)
!!\param totuclfl    clear sky upward flux \f$(w/m^2)\f$
!!\param totdclfl    clear sky downward flux \f$(w/m^2)\f$
!!\param htrcl       clear sky heating rate (k/sec or k/day)
!!\param htrb        spectral band lw heating rate (k/day)
!>\section gen_rtrn rtrn General Algorithm
!! @{
! ----------------------------------
      subroutine rtrn                                                   &
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              & !  ---  inputs
     &       fracs,secdif, nlay,nlp1,                                   &
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       & !  ---  outputs
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are assumed as     !
! randomly overlaping in a vertical colum.                              !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfrc  - real, layer cloud fraction                         0:nlp1 !
!   taucld  - real, layer cloud opt depth                    nbands,nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw,nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw,nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    efclrfr- real, effective clear sky fraction (1-efcldfr)       nlay !
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck fn         1  !
!    totfac - real, gas+cld pade factor, used for planck fn          1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas only            nlay!
!    totsrcu- real, upwd source radiance due to gas+cld             nlay!
!    gassrcd- real, dnwd source radiance due to gas only             1  !
!    totsrcd- real, dnwd source radiance due to gas+cld              1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all Planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 W/m2 within cloudy layers.                !
!  clouds are treated with a random cloud overlap method.               !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(0:nlp1), intent(in) :: cldfrc
      real (kind=kind_phys), dimension(nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(nbands,nlay),intent(in):: taucld
      real (kind=kind_phys), dimension(ngptlw,nlay),intent(in):: fracs, &
     &       tautot

      real (kind=kind_phys), dimension(nbands,0:nlay), intent(in) ::    &
     &       pklev, pklay

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(nlay,nbands),intent(out) :: htrb

      real (kind=kind_phys), dimension(0:nlay), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      real (kind=kind_phys), dimension(0:nlay,nbands) :: clrurad,       &
     &       clrdrad, toturad, totdrad

      real (kind=kind_phys), dimension(nlay)   :: gassrcu, totsrcu,     &
     &       trngas, efclrfr, rfdelp
      real (kind=kind_phys), dimension(0:nlay) :: fnet, fnetc

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, radtotd, radclrd, rad0,          &
     &       clfr, trng, gasu

      integer :: ittot, itgas, ib, ig, k
!
!===> ...  begin here
!
      do ib = 1, NBANDS
        do k = 0, NLAY
          toturad(k,ib) = f_zero
          totdrad(k,ib) = f_zero
          clrurad(k,ib) = f_zero
          clrdrad(k,ib) = f_zero
        enddo
      enddo

      do k = 0, nlay
        totuflux(k) = f_zero
        totdflux(k) = f_zero
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

!  --- ...  loop over all g-points

      do ig = 1, ngptlw
        ib = ngb(ig)

        radtotd = f_zero
        radclrd = f_zero

!> -# Downward radiative transfer loop.

        do k = nlay, 1, -1

!!\n  - clear sky, gases contribution

          odepth = max( f_zero, secdif(ib)*tautot(ig,k) )
          if (odepth <= 0.06) then
            atrgas = odepth - 0.5*odepth*odepth
            trng   = f_one - atrgas
            gasfac = rec_6 * odepth
          else
            tblind = odepth / (bpade + odepth)
            itgas = tblint*tblind + 0.5
            trng  = exp_tbl(itgas)
            atrgas = f_one - trng
            gasfac = tfn_tbl(itgas)
            odepth = tau_tbl(itgas)
          endif

          plfrac = fracs(ig,k)
          blay = pklay(ib,k)

          dplnku = pklev(ib,k  ) - blay
          dplnkd = pklev(ib,k-1) - blay
          bbdgas = plfrac * (blay + dplnkd*gasfac)
          bbugas = plfrac * (blay + dplnku*gasfac)
          gassrcd= bbdgas * atrgas
          gassrcu(k)= bbugas * atrgas
          trngas(k) = trng

!!\n  - total sky, gases+clouds contribution

          clfr = cldfrc(k)
          if (clfr >= eps) then
!!\n  - cloudy layer

            odcld = secdif(ib) * taucld(ib,k)
            efclrfr(k) = f_one-(f_one - exp(-odcld))*clfr
            odtot = odepth + odcld
            if (odtot < 0.06) then
              totfac = rec_6 * odtot
              atrtot = odtot - 0.5*odtot*odtot
            else
              tblind = odtot / (bpade + odtot)
              ittot  = tblint*tblind + 0.5
              totfac = tfn_tbl(ittot)
              atrtot = f_one - exp_tbl(ittot)
            endif

            bbdtot = plfrac * (blay + dplnkd*totfac)
            bbutot = plfrac * (blay + dplnku*totfac)
            totsrcd= bbdtot * atrtot
            totsrcu(k)= bbutot * atrtot

!  --- ...  total sky radiance
            radtotd = radtotd*trng*efclrfr(k) + gassrcd                 &
     &              + clfr*(totsrcd - gassrcd)
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotd = radtotd*trng + gassrcd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

!> -# Compute spectral emissivity & reflectance, include the
!!    contribution of spectrally varying longwave emissivity and
!!     reflection from the surface to the upward radiative transfer.

!     note: spectral and Lambertian reflection are identical for the
!           diffusivity angle flux integration used here.

        reflct = f_one - semiss(ib)
        rad0 = semiss(ib) * fracs(ig,1) * pklay(ib,0)

!> -# Compute total sky radiance.
        radtotu = rad0 + reflct*radtotd
        toturad(0,ib) = toturad(0,ib) + radtotu

!> -# Compute clear sky radiance
        radclru = rad0 + reflct*radclrd
        clrurad(0,ib) = clrurad(0,ib) + radclru

!> -# Upward radiative transfer loop.

        do k = 1, nlay
          clfr = cldfrc(k)
          trng = trngas(k)
          gasu = gassrcu(k)

          if (clfr >= eps) then
!  --- ...  cloudy layer

!  --- ... total sky radiance
            radtotu = radtotu*trng*efclrfr(k) + gasu                    &
     &            + clfr*(totsrcu(k) - gasu)
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          else
!  --- ...  clear layer

!  --- ... total sky radiance
            radtotu = radtotu*trng + gasu
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

      enddo   ! end do_ig_loop

!> -# Process longwave output from band for total and clear streams.
!!    Calculate upward, downward, and net flux.

      flxfac = wtdiff * fluxfac

      do k = 0, nlay
        do ib = 1, nbands
          totuflux(k) = totuflux(k) + toturad(k,ib)
          totdflux(k) = totdflux(k) + totdrad(k,ib)
          totuclfl(k) = totuclfl(k) + clrurad(k,ib)
          totdclfl(k) = totdclfl(k) + clrdrad(k,ib)
        enddo

        totuflux(k) = totuflux(k) * flxfac
        totdflux(k) = totdflux(k) * flxfac
        totuclfl(k) = totuclfl(k) * flxfac
        totdclfl(k) = totdclfl(k) * flxfac
      enddo

!  --- ...  calculate net fluxes and heating rates
      fnet(0) = totuflux(0) - totdflux(0)

      do k = 1, nlay
        rfdelp(k) = heatfac / delp(k)
        fnet(k) = totuflux(k) - totdflux(k)
        htr (k) = (fnet(k-1) - fnet(k)) * rfdelp(k)
      enddo

!! --- ...  optional clear sky heating rates
      if ( lhlw0 ) then
        fnetc(0) = totuclfl(0) - totdclfl(0)

        do k = 1, nlay
          fnetc(k) = totuclfl(k) - totdclfl(k)
          htrcl(k) = (fnetc(k-1) - fnetc(k)) * rfdelp(k)
        enddo
      endif

!! --- ...  optional spectral band heating rates
      if ( lhlwb ) then
        do ib = 1, nbands
          fnet(0) = (toturad(0,ib) - totdrad(0,ib)) * flxfac

          do k = 1, nlay
            fnet(k) = (toturad(k,ib) - totdrad(k,ib)) * flxfac
            htrb(k,ib) = (fnet(k-1) - fnet(k)) * rfdelp(k)
          enddo
        enddo
      endif

! ..................................
      end subroutine rtrn
!! @}
! ----------------------------------


!>\ingroup module_radlw_main
!> This subroutine computes the upward/downward radiative fluxes, and
!! heating rates for both clear or cloudy atmosphere. Clouds are
!! assumed as in maximum-randomly overlaping in a vertical column.
!!\param semiss        lw surface emissivity
!!\param delp          layer pressure thickness (mb)
!!\param cldfrc        layer cloud fraction
!!\param taucld        layer cloud opt depth
!!\param tautot        total optical depth (gas+aerosols)
!!\param pklay         integrated planck func at lay temp
!!\param pklev         integrated planck func at lev temp
!!\param fracs         planck fractions
!!\param secdif        secant of diffusivity angle
!!\param nlay          number of vertical layers
!!\param nlp1          number of vertical levels (interfaces)
!!\param totuflux      total sky upward flux (\f$w/m^2\f$)
!!\param totdflux      total sky downward flux (\f$w/m^2\f$)
!!\param htr           total sky heating rate (k/sec or k/day)
!!\param totuclfl      clear sky upward flux (\f$w/m^2\f$)
!!\param totdclfl      clear sky downward flux (\f$w/m^2\f$)
!!\param htrcl         clear sky heating rate (k/sec or k/day)
!!\param htrb          spectral band lw heating rate (k/day)
!!\section gen_rtrnmr rtrnmr General Algorithm
!> @{
! ----------------------------------
      subroutine rtrnmr                                                 &
     &     ( semiss,delp,cldfrc,taucld,tautot,pklay,pklev,              &!  ---  inputs
     &       fracs,secdif, nlay,nlp1,                                   &
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       &  !  ---  outputs:
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are assumed as in  !
! maximum-randomly overlaping in a vertical colum.                      !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfrc  - real, layer cloud fraction                         0:nlp1 !
!   taucld  - real, layer cloud opt depth                    nbands,nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw,nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw,nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck fn         1  !
!    totfac - real, gas+cld pade factor, used for planck fn          1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas only            nlay!
!    totsrcu- real, upwd source radiance due to gas + cld           nlay!
!    gassrcd- real, dnwd source radiance due to gas only             1  !
!    totsrcd- real, dnwd source radiance due to gas + cld            1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all Planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 W/m2 within cloudy layers.                !
!  clouds are treated with a maximum-random cloud overlap method.       !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(0:nlp1), intent(in) :: cldfrc
      real (kind=kind_phys), dimension(nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(nbands,nlay),intent(in):: taucld
      real (kind=kind_phys), dimension(ngptlw,nlay),intent(in):: fracs, &
     &       tautot

      real (kind=kind_phys), dimension(nbands,0:nlay), intent(in) ::    &
     &       pklev, pklay

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(nlay,nbands),intent(out) :: htrb

      real (kind=kind_phys), dimension(0:nlay), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      real (kind=kind_phys), dimension(0:nlay,nbands) :: clrurad,       &
     &       clrdrad, toturad, totdrad

      real (kind=kind_phys), dimension(nlay)   :: gassrcu, totsrcu,     &
     &       trngas, trntot, rfdelp
      real (kind=kind_phys), dimension(0:nlay) :: fnet, fnetc

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, radtotd, radclrd, rad0, rad,     &
     &       totradd, clrradd, totradu, clrradu, fmax, fmin, rat1, rat2,&
     &       radmod, clfr, trng, trnt, gasu, totu

      integer :: ittot, itgas, ib, ig, k

!  dimensions for cloud overlap adjustment
      real (kind=kind_phys), dimension(nlp1) :: faccld1u, faccld2u,     &
     &        facclr1u, facclr2u, faccmb1u, faccmb2u
      real (kind=kind_phys), dimension(0:nlay) :: faccld1d, faccld2d,   &
     &        facclr1d, facclr2d, faccmb1d, faccmb2d

      logical :: lstcldu(nlay), lstcldd(nlay)
!
!===> ...  begin here
!
      do k = 1, nlp1
        faccld1u(k) = f_zero
        faccld2u(k) = f_zero
        facclr1u(k) = f_zero
        facclr2u(k) = f_zero
        faccmb1u(k) = f_zero
        faccmb2u(k) = f_zero
      enddo

      lstcldu(1) = cldfrc(1) > eps
      rat1 = f_zero
      rat2 = f_zero

      do k = 1, nlay-1

        lstcldu(k+1) = cldfrc(k+1)>eps .and. cldfrc(k)<=eps

        if (cldfrc(k) > eps) then

!> -# Setup maximum/random cloud overlap.

          if (cldfrc(k+1) >= cldfrc(k)) then
            if (lstcldu(k)) then
              if (cldfrc(k) < f_one) then
                facclr2u(k+1) = (cldfrc(k+1) - cldfrc(k))               &
     &                        / (f_one - cldfrc(k))
              endif
              facclr2u(k) = f_zero
              faccld2u(k) = f_zero
            else
              fmax = max(cldfrc(k), cldfrc(k-1))
              if (cldfrc(k+1) > fmax) then
                facclr1u(k+1) = rat2
                facclr2u(k+1) = (cldfrc(k+1) - fmax)/(f_one - fmax)
              elseif (cldfrc(k+1) < fmax) then
                facclr1u(k+1) = (cldfrc(k+1) - cldfrc(k))               &
     &                        / (cldfrc(k-1) - cldfrc(k))
              else
                facclr1u(k+1) = rat2
              endif
            endif

            if (facclr1u(k+1)>f_zero .or. facclr2u(k+1)>f_zero) then
              rat1 = f_one
              rat2 = f_zero
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          else
            if (lstcldu(k)) then
              faccld2u(k+1) = (cldfrc(k) - cldfrc(k+1)) / cldfrc(k)
              facclr2u(k) = f_zero
              faccld2u(k) = f_zero
            else
              fmin = min(cldfrc(k), cldfrc(k-1))
              if (cldfrc(k+1) <= fmin) then
                faccld1u(k+1) = rat1
                faccld2u(k+1) = (fmin - cldfrc(k+1)) / fmin
              else
                faccld1u(k+1) = (cldfrc(k) - cldfrc(k+1))               &
     &                        / (cldfrc(k) - fmin)
              endif
            endif

            if (faccld1u(k+1)>f_zero .or. faccld2u(k+1)>f_zero) then
              rat1 = f_zero
              rat2 = f_one
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          endif

          faccmb1u(k+1) = facclr1u(k+1) * faccld2u(k) * cldfrc(k-1)
          faccmb2u(k+1) = faccld1u(k+1) * facclr2u(k)                   &
     &                  * (f_one - cldfrc(k-1))
        endif

      enddo

      do k = 0, nlay
        faccld1d(k) = f_zero
        faccld2d(k) = f_zero
        facclr1d(k) = f_zero
        facclr2d(k) = f_zero
        faccmb1d(k) = f_zero
        faccmb2d(k) = f_zero
      enddo

      lstcldd(nlay) = cldfrc(nlay) > eps
      rat1 = f_zero
      rat2 = f_zero

      do k = nlay, 2, -1

        lstcldd(k-1) = cldfrc(k-1) > eps .and. cldfrc(k)<=eps

        if (cldfrc(k) > eps) then

          if (cldfrc(k-1) >= cldfrc(k)) then
            if (lstcldd(k)) then
              if (cldfrc(k) < f_one) then
                facclr2d(k-1) = (cldfrc(k-1) - cldfrc(k))               &
     &                        / (f_one - cldfrc(k))
              endif

              facclr2d(k) = f_zero
              faccld2d(k) = f_zero
            else
              fmax = max(cldfrc(k), cldfrc(k+1))

              if (cldfrc(k-1) > fmax) then
                facclr1d(k-1) = rat2
                facclr2d(k-1) = (cldfrc(k-1) - fmax) / (f_one - fmax)
              elseif (cldfrc(k-1) < fmax) then
                facclr1d(k-1) = (cldfrc(k-1) - cldfrc(k))               &
     &                        / (cldfrc(k+1) - cldfrc(k))
              else
                facclr1d(k-1) = rat2
              endif
            endif

            if (facclr1d(k-1)>f_zero .or. facclr2d(k-1)>f_zero) then
              rat1 = f_one
              rat2 = f_zero
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          else
            if (lstcldd(k)) then
              faccld2d(k-1) = (cldfrc(k) - cldfrc(k-1)) / cldfrc(k)
              facclr2d(k) = f_zero
              faccld2d(k) = f_zero
            else
              fmin = min(cldfrc(k), cldfrc(k+1))

              if (cldfrc(k-1) <= fmin) then
                faccld1d(k-1) = rat1
                faccld2d(k-1) = (fmin - cldfrc(k-1)) / fmin
              else
                faccld1d(k-1) = (cldfrc(k) - cldfrc(k-1))               &
     &                        / (cldfrc(k) - fmin)
              endif
            endif

            if (faccld1d(k-1)>f_zero .or. faccld2d(k-1)>f_zero) then
              rat1 = f_zero
              rat2 = f_one
            else
              rat1 = f_zero
              rat2 = f_zero
            endif
          endif

          faccmb1d(k-1) = facclr1d(k-1) * faccld2d(k) * cldfrc(k+1)
          faccmb2d(k-1) = faccld1d(k-1) * facclr2d(k)                   &
     &                  * (f_one - cldfrc(k+1))
        endif

      enddo

!> -# Initialize for radiative transfer

      do ib = 1, NBANDS
        do k = 0, NLAY
          toturad(k,ib) = f_zero
          totdrad(k,ib) = f_zero
          clrurad(k,ib) = f_zero
          clrdrad(k,ib) = f_zero
        enddo
      enddo

      do k = 0, nlay
        totuflux(k) = f_zero
        totdflux(k) = f_zero
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

!  --- ...  loop over all g-points

      do ig = 1, ngptlw
        ib = ngb(ig)

        radtotd = f_zero
        radclrd = f_zero

!> -# Downward radiative transfer loop:

        do k = nlay, 1, -1

!  --- ...  clear sky, gases contribution

          odepth = max( f_zero, secdif(ib)*tautot(ig,k) )
          if (odepth <= 0.06) then
            atrgas = odepth - 0.5*odepth*odepth
            trng   = f_one - atrgas
            gasfac = rec_6 * odepth
          else
            tblind = odepth / (bpade + odepth)
            itgas = tblint*tblind + 0.5
            trng  = exp_tbl(itgas)
            atrgas = f_one - trng
            gasfac = tfn_tbl(itgas)
            odepth = tau_tbl(itgas)
          endif

          plfrac = fracs(ig,k)
          blay = pklay(ib,k)

          dplnku = pklev(ib,k  ) - blay
          dplnkd = pklev(ib,k-1) - blay
          bbdgas = plfrac * (blay + dplnkd*gasfac)
          bbugas = plfrac * (blay + dplnku*gasfac)
          gassrcd   = bbdgas * atrgas
          gassrcu(k)= bbugas * atrgas
          trngas(k) = trng

!  --- ...  total sky, gases+clouds contribution

          clfr = cldfrc(k)
          if (lstcldd(k)) then
            totradd = clfr * radtotd
            clrradd = radtotd - totradd
            rad = f_zero
          endif

          if (clfr >= eps) then
!>  - cloudy layer

            odcld = secdif(ib) * taucld(ib,k)
            odtot = odepth + odcld
            if (odtot < 0.06) then
              totfac = rec_6 * odtot
              atrtot = odtot - 0.5*odtot*odtot
              trnt   = f_one - atrtot
            else
              tblind = odtot / (bpade + odtot)
              ittot  = tblint*tblind + 0.5
              totfac = tfn_tbl(ittot)
              trnt   = exp_tbl(ittot)
              atrtot = f_one - trnt
            endif

            bbdtot = plfrac * (blay + dplnkd*totfac)
            bbutot = plfrac * (blay + dplnku*totfac)
            totsrcd   = bbdtot * atrtot
            totsrcu(k)= bbutot * atrtot
            trntot(k) = trnt

            totradd = totradd*trnt + clfr*totsrcd
            clrradd = clrradd*trng + (f_one - clfr)*gassrcd

!>  - total sky radiance
            radtotd = totradd + clrradd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!>  - clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

            radmod = rad*(facclr1d(k-1)*trng + faccld1d(k-1)*trnt)      &
     &             - faccmb1d(k-1)*gassrcd + faccmb2d(k-1)*totsrcd

            rad = -radmod + facclr2d(k-1)*(clrradd + radmod)            &
     &                    - faccld2d(k-1)*(totradd - radmod)
            totradd = totradd + rad
            clrradd = clrradd - rad

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotd = radtotd*trng + gassrcd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

!> -# Compute spectral emissivity & reflectance, include the
!!    contribution of spectrally varying longwave emissivity and
!!    reflection from the surface to the upward radiative transfer.

!     note: spectral and Lambertian reflection are identical for the
!           diffusivity angle flux integration used here.

        reflct = f_one - semiss(ib)
        rad0 = semiss(ib) * fracs(ig,1) * pklay(ib,0)

!> -# Compute total sky radiance.
        radtotu = rad0 + reflct*radtotd
        toturad(0,ib) = toturad(0,ib) + radtotu

!> -# Compute clear sky radiance.
        radclru = rad0 + reflct*radclrd
        clrurad(0,ib) = clrurad(0,ib) + radclru

!> -# Upward radiative transfer loop:

        do k = 1, nlay

          clfr = cldfrc(k)
          trng = trngas(k)
          gasu = gassrcu(k)

          if (lstcldu(k)) then
            totradu = clfr * radtotu
            clrradu = radtotu - totradu
            rad = f_zero
          endif

          if (clfr >= eps) then
!>  - cloudy layer radiance

            trnt = trntot(k)
            totu = totsrcu(k)
            totradu = totradu*trnt + clfr*totu
            clrradu = clrradu*trng + (f_one - clfr)*gasu

!>  - total sky radiance
            radtotu = totradu + clrradu
            toturad(k,ib) = toturad(k,ib) + radtotu

!>  - clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

            radmod = rad*(facclr1u(k+1)*trng + faccld1u(k+1)*trnt)      &
     &             - faccmb1u(k+1)*gasu + faccmb2u(k+1)*totu
            rad = -radmod + facclr2u(k+1)*(clrradu + radmod)            &
     &                    - faccld2u(k+1)*(totradu - radmod)
            totradu = totradu + rad
            clrradu = clrradu - rad

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotu = radtotu*trng + gasu
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ...  clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          endif   ! end if_clfr_block

        enddo   ! end do_k_loop

      enddo   ! end do_ig_loop

!> -# Process longwave output from band for total and clear streams.
!! calculate upward, downward, and net flux.

      flxfac = wtdiff * fluxfac

      do k = 0, nlay
        do ib = 1, nbands
          totuflux(k) = totuflux(k) + toturad(k,ib)
          totdflux(k) = totdflux(k) + totdrad(k,ib)
          totuclfl(k) = totuclfl(k) + clrurad(k,ib)
          totdclfl(k) = totdclfl(k) + clrdrad(k,ib)
        enddo

        totuflux(k) = totuflux(k) * flxfac
        totdflux(k) = totdflux(k) * flxfac
        totuclfl(k) = totuclfl(k) * flxfac
        totdclfl(k) = totdclfl(k) * flxfac
      enddo

!  --- ...  calculate net fluxes and heating rates
      fnet(0) = totuflux(0) - totdflux(0)

      do k = 1, nlay
        rfdelp(k) = heatfac / delp(k)
        fnet(k) = totuflux(k) - totdflux(k)
        htr (k) = (fnet(k-1) - fnet(k)) * rfdelp(k)
      enddo

!! --- ...  optional clear sky heating rates
      if ( lhlw0 ) then
        fnetc(0) = totuclfl(0) - totdclfl(0)

        do k = 1, nlay
          fnetc(k) = totuclfl(k) - totdclfl(k)
          htrcl(k) = (fnetc(k-1) - fnetc(k)) * rfdelp(k)
        enddo
      endif

!! --- ...  optional spectral band heating rates
      if ( lhlwb ) then
        do ib = 1, nbands
          fnet(0) = (toturad(0,ib) - totdrad(0,ib)) * flxfac

          do k = 1, nlay
            fnet(k) = (toturad(k,ib) - totdrad(k,ib)) * flxfac
            htrb(k,ib) = (fnet(k-1) - fnet(k)) * rfdelp(k)
          enddo
        enddo
      endif

! .................................
      end subroutine rtrnmr
! ---------------------------------
!> @}

!>\ingroup module_radlw_main
!> \brief This subroutine computes the upward/downward radiative fluxes, and
!! heating rates for both clear or cloudy atmosphere.Clouds are treated
!! with the mcica stochastic approach.
!!
!!\param semiss       lw surface emissivity
!!\param delp         layer pressure thickness (mb)
!!\param cldfmc       layer cloud fraction (sub-column)
!!\param taucld       layer cloud opt depth
!!\param tautot       total optical depth (gas+aerosols)
!!\param pklay        integrated planck func at lay temp
!!\param pklev        integrated planck func at lev temp
!!\param fracs        planck fractions
!!\param secdif       secant of diffusivity angle
!!\param nlay         number of vertical layers
!!\param nlp1         number of vertical levels (interfaces)
!!\param totuflux     total sky upward flux \f$(w/m^2)\f$
!!\param totdflux     total sky downward flux \f$(w/m^2)\f$
!!\param htr          total sky heating rate (k/sec or k/day)
!!\param totuclfl     clear sky upward flux \f$(w/m^2)\f$
!!\param totdclfl     clear sky downward flux \f$(w/m^2)\f$
!!\param htrcl        clear sky heating rate (k/sec or k/day)
!!\param htrb         spectral band lw heating rate (k/day)
!!\section gen_rtrnmc rtrnmc General Algorithm
!> @{
! ---------------------------------
      subroutine rtrnmc                                                 &
     &     ( semiss,delp,cldfmc,taucld,tautot,pklay,pklev,              & !  ---  inputs:
     &       fracs,secdif, nlay,nlp1,                                   &
     &       totuflux,totdflux,htr, totuclfl,totdclfl,htrcl, htrb       & !  ---  outputs:
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! purpose:  compute the upward/downward radiative fluxes, and heating   !
! rates for both clear or cloudy atmosphere.  clouds are treated with   !
! the mcica stochastic approach.                                        !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                     -size-   !
!   semiss  - real, lw surface emissivity                         nbands!
!   delp    - real, layer pressure thickness (mb)                  nlay !
!   cldfmc  - real, layer cloud fraction (sub-column)        ngptlw*nlay!
!   taucld  - real, layer cloud opt depth                    nbands*nlay!
!   tautot  - real, total optical depth (gas+aerosols)       ngptlw*nlay!
!   pklay   - real, integrated planck func at lay temp     nbands*0:nlay!
!   pklev   - real, integrated planck func at lev temp     nbands*0:nlay!
!   fracs   - real, planck fractions                         ngptlw*nlay!
!   secdif  - real, secant of diffusivity angle                   nbands!
!   nlay    - integer, number of vertical layers                    1   !
!   nlp1    - integer, number of vertical levels (interfaces)       1   !
!                                                                       !
!  outputs:                                                             !
!   totuflux- real, total sky upward flux (w/m2)                 0:nlay !
!   totdflux- real, total sky downward flux (w/m2)               0:nlay !
!   htr     - real, total sky heating rate (k/sec or k/day)        nlay !
!   totuclfl- real, clear sky upward flux (w/m2)                 0:nlay !
!   totdclfl- real, clear sky downward flux (w/m2)               0:nlay !
!   htrcl   - real, clear sky heating rate (k/sec or k/day)        nlay !
!   htrb    - real, spectral band lw heating rate (k/day)    nlay*nbands!
!                                                                       !
!  module veriables:                                                    !
!   ngb     - integer, band index for each g-value                ngptlw!
!   fluxfac - real, conversion factor for fluxes (pi*2.e4)           1  !
!   heatfac - real, conversion factor for heating rates (g/cp*1e-2)  1  !
!   tblint  - real, conversion factor for look-up tbl (float(ntbl)   1  !
!   bpade   - real, pade approx constant (1/0.278)                   1  !
!   wtdiff  - real, weight for radiance to flux conversion           1  !
!   ntbl    - integer, dimension of look-up tables                   1  !
!   tau_tbl - real, clr-sky opt dep lookup table                 0:ntbl !
!   exp_tbl - real, transmittance lookup table                   0:ntbl !
!   tfn_tbl - real, tau transition function                      0:ntbl !
!                                                                       !
!  local variables:                                                     !
!    itgas  - integer, index for gases contribution look-up table    1  !
!    ittot  - integer, index for gases plus clouds  look-up table    1  !
!    reflct - real, surface reflectance                              1  !
!    atrgas - real, gaseous absorptivity                             1  !
!    atrtot - real, gaseous and cloud absorptivity                   1  !
!    odcld  - real, cloud optical depth                              1  !
!    efclrfr- real, effective clear sky fraction (1-efcldfr)        nlay!
!    odepth - real, optical depth of gaseous only                    1  !
!    odtot  - real, optical depth of gas and cloud                   1  !
!    gasfac - real, gas-only pade factor, used for planck function   1  !
!    totfac - real, gas and cloud pade factor, used for planck fn    1  !
!    bbdgas - real, gas-only planck function for downward rt         1  !
!    bbugas - real, gas-only planck function for upward rt           1  !
!    bbdtot - real, gas and cloud planck function for downward rt    1  !
!    bbutot - real, gas and cloud planck function for upward rt      1  !
!    gassrcu- real, upwd source radiance due to gas                 nlay!
!    totsrcu- real, upwd source radiance due to gas+cld             nlay!
!    gassrcd- real, dnwd source radiance due to gas                  1  !
!    totsrcd- real, dnwd source radiance due to gas+cld              1  !
!    radtotu- real, spectrally summed total sky upwd radiance        1  !
!    radclru- real, spectrally summed clear sky upwd radiance        1  !
!    radtotd- real, spectrally summed total sky dnwd radiance        1  !
!    radclrd- real, spectrally summed clear sky dnwd radiance        1  !
!    toturad- real, total sky upward radiance by layer     0:nlay*nbands!
!    clrurad- real, clear sky upward radiance by layer     0:nlay*nbands!
!    totdrad- real, total sky downward radiance by layer   0:nlay*nbands!
!    clrdrad- real, clear sky downward radiance by layer   0:nlay*nbands!
!    fnet   - real, net longwave flux (w/m2)                     0:nlay !
!    fnetc  - real, clear sky net longwave flux (w/m2)           0:nlay !
!                                                                       !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  original version:   e. j. mlawer, et al. rrtm_v3.0                   !
!  revision for gcms:  michael j. iacono; october, 2002                 !
!  revision for f90:   michael j. iacono; june, 2006                    !
!                                                                       !
!  this program calculates the upward fluxes, downward fluxes, and      !
!  heating rates for an arbitrary clear or cloudy atmosphere. the input !
!  to this program is the atmospheric profile, all Planck function      !
!  information, and the cloud fraction by layer.  a variable diffusivity!
!  angle (secdif) is used for the angle integration. bands 2-3 and 5-9  !
!  use a value for secdif that varies from 1.50 to 1.80 as a function   !
!  of the column water vapor, and other bands use a value of 1.66.  the !
!  gaussian weight appropriate to this angle (wtdiff=0.5) is applied    !
!  here.  note that use of the emissivity angle for the flux integration!
!  can cause errors of 1 to 4 W/m2 within cloudy layers.                !
!  clouds are treated with the mcica stochastic approach and            !
!  maximum-random cloud overlap.                                        !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nbands), intent(in) :: semiss,   &
     &       secdif
      real (kind=kind_phys), dimension(nlay),   intent(in) :: delp

      real (kind=kind_phys), dimension(nbands,nlay),intent(in):: taucld
      real (kind=kind_phys), dimension(ngptlw,nlay),intent(in):: fracs, &
     &       tautot, cldfmc

      real (kind=kind_phys), dimension(nbands,0:nlay), intent(in) ::    &
     &       pklev, pklay

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay), intent(out) :: htr, htrcl

      real (kind=kind_phys), dimension(nlay,nbands),intent(out) :: htrb

      real (kind=kind_phys), dimension(0:nlay), intent(out) ::          &
     &       totuflux, totdflux, totuclfl, totdclfl

!  ---  locals:
      real (kind=kind_phys), parameter :: rec_6 = 0.166667

      real (kind=kind_phys), dimension(0:nlay,nbands) :: clrurad,       &
     &       clrdrad, toturad, totdrad

      real (kind=kind_phys), dimension(nlay)   :: gassrcu, totsrcu,     &
     &       trngas, efclrfr, rfdelp
      real (kind=kind_phys), dimension(0:nlay) :: fnet, fnetc

      real (kind=kind_phys) :: totsrcd, gassrcd, tblind, odepth, odtot, &
     &       odcld, atrtot, atrgas, reflct, totfac, gasfac, flxfac,     &
     &       plfrac, blay, bbdgas, bbdtot, bbugas, bbutot, dplnku,      &
     &       dplnkd, radtotu, radclru, radtotd, radclrd, rad0,          &
     &       clfm, trng, gasu

      integer :: ittot, itgas, ib, ig, k
!
!===> ...  begin here
!
      do ib = 1, NBANDS
        do k = 0, NLAY
          toturad(k,ib) = f_zero
          totdrad(k,ib) = f_zero
          clrurad(k,ib) = f_zero
          clrdrad(k,ib) = f_zero
        enddo
      enddo

      do k = 0, nlay
        totuflux(k) = f_zero
        totdflux(k) = f_zero
        totuclfl(k) = f_zero
        totdclfl(k) = f_zero
      enddo

!  --- ...  loop over all g-points

      do ig = 1, ngptlw
        ib = ngb(ig)

        radtotd = f_zero
        radclrd = f_zero

!> -# Downward radiative transfer loop.
!!\n  - Clear sky, gases contribution
!!\n  - Total sky, gases+clouds contribution
!!\n  - Cloudy layer
!!\n  - Total sky radiance
!!\n  - Clear sky radiance

        do k = nlay, 1, -1

!  --- ...  clear sky, gases contribution

          odepth = max( f_zero, secdif(ib)*tautot(ig,k) )
          if (odepth <= 0.06) then
            atrgas = odepth - 0.5*odepth*odepth
            trng   = f_one - atrgas
            gasfac = rec_6 * odepth
          else
            tblind = odepth / (bpade + odepth)
            itgas = tblint*tblind + 0.5
            trng  = exp_tbl(itgas)
            atrgas = f_one - trng
            gasfac = tfn_tbl(itgas)
            odepth = tau_tbl(itgas)
          endif

          plfrac = fracs(ig,k)
          blay = pklay(ib,k)

          dplnku = pklev(ib,k  ) - blay
          dplnkd = pklev(ib,k-1) - blay
          bbdgas = plfrac * (blay + dplnkd*gasfac)
          bbugas = plfrac * (blay + dplnku*gasfac)
          gassrcd= bbdgas * atrgas
          gassrcu(k)= bbugas * atrgas
          trngas(k) = trng

!  --- ...  total sky, gases+clouds contribution

          clfm = cldfmc(ig,k)
          if (clfm >= eps) then
!  --- ...  cloudy layer

            odcld = secdif(ib) * taucld(ib,k)
            efclrfr(k) = f_one - (f_one - exp(-odcld))*clfm
            odtot = odepth + odcld
            if (odtot < 0.06) then
              totfac = rec_6 * odtot
              atrtot = odtot - 0.5*odtot*odtot
            else
              tblind = odtot / (bpade + odtot)
              ittot  = tblint*tblind + 0.5
              totfac = tfn_tbl(ittot)
              atrtot = f_one - exp_tbl(ittot)
            endif

            bbdtot = plfrac * (blay + dplnkd*totfac)
            bbutot = plfrac * (blay + dplnku*totfac)
            totsrcd= bbdtot * atrtot
            totsrcu(k)= bbutot * atrtot

!  --- ...  total sky radiance
            radtotd = radtotd*trng*efclrfr(k) + gassrcd                 &
     &              + clfm*(totsrcd - gassrcd)
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          else
!  --- ...  clear layer

!  --- ...  total sky radiance
            radtotd = radtotd*trng + gassrcd
            totdrad(k-1,ib) = totdrad(k-1,ib) + radtotd

!  --- ...  clear sky radiance
            radclrd = radclrd*trng + gassrcd
            clrdrad(k-1,ib) = clrdrad(k-1,ib) + radclrd

          endif   ! end if_clfm_block

        enddo   ! end do_k_loop

!> -# Compute spectral emissivity & reflectance, include the
!!    contribution of spectrally varying longwave emissivity and
!!    reflection from the surface to the upward radiative transfer.

!     note: spectral and Lambertian reflection are identical for the
!           diffusivity angle flux integration used here.

        reflct = f_one - semiss(ib)
        rad0 = semiss(ib) * fracs(ig,1) * pklay(ib,0)

!> -# Compute total sky radiance.
        radtotu = rad0 + reflct*radtotd
        toturad(0,ib) = toturad(0,ib) + radtotu

!> -# Compute clear sky radiance.
        radclru = rad0 + reflct*radclrd
        clrurad(0,ib) = clrurad(0,ib) + radclru

!> -# Upward radiative transfer loop.
!!\n  - Compute total sky radiance
!!\n  - Compute clear sky radiance

!          toturad holds summed radiance for total sky stream
!          clrurad holds summed radiance for clear sky stream

        do k = 1, nlay
          clfm = cldfmc(ig,k)
          trng = trngas(k)
          gasu = gassrcu(k)

          if (clfm > eps) then
!  --- ...  cloudy layer

!  --- ... total sky radiance
            radtotu = radtotu*trng*efclrfr(k) + gasu                    &
     &              + clfm*(totsrcu(k) - gasu)
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          else
!  --- ...  clear layer

!  --- ... total sky radiance
            radtotu = radtotu*trng + gasu
            toturad(k,ib) = toturad(k,ib) + radtotu

!  --- ... clear sky radiance
            radclru = radclru*trng + gasu
            clrurad(k,ib) = clrurad(k,ib) + radclru

          endif   ! end if_clfm_block

        enddo   ! end do_k_loop

      enddo   ! end do_ig_loop

!> -# Process longwave output from band for total and clear streams.
!!    Calculate upward, downward, and net flux.

      flxfac = wtdiff * fluxfac

      do k = 0, nlay
        do ib = 1, nbands
          totuflux(k) = totuflux(k) + toturad(k,ib)
          totdflux(k) = totdflux(k) + totdrad(k,ib)
          totuclfl(k) = totuclfl(k) + clrurad(k,ib)
          totdclfl(k) = totdclfl(k) + clrdrad(k,ib)
        enddo

        totuflux(k) = totuflux(k) * flxfac
        totdflux(k) = totdflux(k) * flxfac
        totuclfl(k) = totuclfl(k) * flxfac
        totdclfl(k) = totdclfl(k) * flxfac
      enddo

!> -# Calculate net fluxes and heating rates.
      fnet(0) = totuflux(0) - totdflux(0)

      do k = 1, nlay
        rfdelp(k) = heatfac / delp(k)
        fnet(k) = totuflux(k) - totdflux(k)
        htr (k) = (fnet(k-1) - fnet(k)) * rfdelp(k)
      enddo

!> -# Optional clear sky heating rates.
      if ( lhlw0 ) then
        fnetc(0) = totuclfl(0) - totdclfl(0)

        do k = 1, nlay
          fnetc(k) = totuclfl(k) - totdclfl(k)
          htrcl(k) = (fnetc(k-1) - fnetc(k)) * rfdelp(k)
        enddo
      endif

!> -# Optional spectral band heating rates.
      if ( lhlwb ) then
        do ib = 1, nbands
          fnet(0) = (toturad(0,ib) - totdrad(0,ib)) * flxfac

          do k = 1, nlay
            fnet(k) = (toturad(k,ib) - totdrad(k,ib)) * flxfac
            htrb(k,ib) = (fnet(k-1) - fnet(k)) * rfdelp(k)
          enddo
        enddo
      endif

! ..................................
      end subroutine rtrnmc
! ----------------------------------
!> @}

!>\ingroup module_radlw_main
!>\brief This subroutine contains optical depths developed for the rapid
!! radiative transfer model.
!!
!! It contains the subroutines \a taugbn (where n goes from
!! 1 to 16). \a taugbn calculates the optical depths and planck fractions
!! per g-value and layer for band n.
!!\param laytrop          tropopause layer index (unitless) layer at
!!                        which switch is made for key species
!!\param pavel            layer pressures (mb)
!!\param coldry           column amount for dry air \f$(mol/cm^2)\f$
!!\param colamt           column amounts of h2o, co2, o3, n2o, ch4,o2,
!!                        co \f$(mol/cm^2)\f$
!!\param colbrd           column amount of broadening gases
!!\param wx               cross-section amounts \f$(mol/cm^2)\f$
!!\param tauaer           aerosol optical depth
!!\param rfrate           reference ratios of binary species parameter
!!\n                      (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,4-h2o/ch4,
!!                                 5-n2o/co2,6-o3/co2
!!\n                      (:,:,n)n=1,2: the rates of ref press at the 2
!!                                 sides of the layer
!!\param facij            factors multiply the reference ks, i,j of 0/1
!!                        for lower/higher of the 2 appropriate
!!                        temperatures and altitudes
!!\param jp               index of lower reference pressure
!!\param jt, jt1          indices of lower reference temperatures for
!!                        pressure levels jp and jp+1, respectively
!!\param selffac          scale factor for water vapor self-continuum
!!                        equals (water vapor density)/(atmospheric
!!                        density at 296k and 1013 mb)
!!\param selffrac         factor for temperature interpolation of
!!                        reference water vapor self-continuum data
!!\param indself          index of lower reference temperature for the
!!                        self-continuum interpolation
!!\param forfac           scale factor for w. v. foreign-continuum
!!\param forfrac          factor for temperature interpolation of
!!                        reference w.v. foreign-continuum data
!!\param indfor           index of lower reference temperature for the
!!                        foreign-continuum interpolation
!!\param minorfrac        factor for minor gases
!!\param scaleminor,scaleminorn2       scale factors for minor gases
!!\param indminor         index of lower reference temperature for
!!                        minor gases
!!\param nlay             total number of layers
!!\param fracs            planck fractions
!!\param tautot           total optical depth (gas+aerosols)
!>\section taumol_gen taumol General Algorithm
!! @{
!! subprograms called:  taugb## (## = 01 -16) 
      subroutine taumol                                                 &
     &     ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              & !  ---  inputs
     &       rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  &
     &       selffac,selffrac,indself,forfac,forfrac,indfor,            &
     &       minorfrac,scaleminor,scaleminorn2,indminor,                &
     &       nlay,                                                      &
     &       fracs, tautot                                              & !  ---  outputs
     &     )

!  ************    original subprogram description    ***************   !
!                                                                       !
!                  optical depths developed for the                     !
!                                                                       !
!                rapid radiative transfer model (rrtm)                  !
!                                                                       !
!            atmospheric and environmental research, inc.               !
!                        131 hartwell avenue                            !
!                        lexington, ma 02421                            !
!                                                                       !
!                           eli j. mlawer                               !
!                         jennifer delamere                             !
!                         steven j. taubman                             !
!                         shepard a. clough                             !
!                                                                       !
!                       email:  mlawer@aer.com                          !
!                       email:  jdelamer@aer.com                        !
!                                                                       !
!        the authors wish to acknowledge the contributions of the       !
!        following people:  karen cady-pereira, patrick d. brown,       !
!        michael j. iacono, ronald e. farren, luke chen,                !
!        robert bergstrom.                                              !
!                                                                       !
!  revision for g-point reduction: michael j. iacono; aer, inc.         !
!                                                                       !
!     taumol                                                            !
!                                                                       !
!     this file contains the subroutines taugbn (where n goes from      !
!     1 to 16).  taugbn calculates the optical depths and planck        !
!     fractions per g-value and layer for band n.                       !
!                                                                       !
!  *******************************************************************  !
!  ==================   program usage description   ==================  !
!                                                                       !
!    call  taumol                                                       !
!       inputs:                                                         !
!          ( laytrop,pavel,coldry,colamt,colbrd,wx,tauaer,              !
!            rfrate,fac00,fac01,fac10,fac11,jp,jt,jt1,                  !
!            selffac,selffrac,indself,forfac,forfrac,indfor,            !
!            minorfrac,scaleminor,scaleminorn2,indminor,                !
!            nlay,                                                      !
!       outputs:                                                        !
!            fracs, tautot )                                            !
!                                                                       !
!  subprograms called:  taugb## (## = 01 -16)                           !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!     laytrop   - integer, tropopause layer index (unitless)        1   !
!                   layer at which switch is made for key species       !
!     pavel     - real, layer pressures (mb)                       nlay !
!     coldry    - real, column amount for dry air (mol/cm2)        nlay !
!     colamt    - real, column amounts of h2o, co2, o3, n2o, ch4,       !
!                   o2, co (mol/cm**2)                       nlay*maxgas!
!     colbrd    - real, column amount of broadening gases          nlay !
!     wx        - real, cross-section amounts(mol/cm2)      nlay*maxxsec!
!     tauaer    - real, aerosol optical depth               nbands*nlay !
!     rfrate    - real, reference ratios of binary species parameter    !
!     (:,m,:)m=1-h2o/co2,2-h2o/o3,3-h2o/n2o,4-h2o/ch4,5-n2o/co2,6-o3/co2!
!     (:,:,n)n=1,2: the rates of ref press at the 2 sides of the layer  !
!                                                          nlay*nrates*2!
!     facij     - real, factors multiply the reference ks, i,j of 0/1   !
!                   for lower/higher of the 2 appropriate temperatures  !
!                   and altitudes                                  nlay !
!     jp        - real, index of lower reference pressure          nlay !
!     jt, jt1   - real, indices of lower reference temperatures    nlay !
!                   for pressure levels jp and jp+1, respectively       !
!     selffac   - real, scale factor for water vapor self-continuum     !
!                   equals (water vapor density)/(atmospheric density   !
!                   at 296k and 1013 mb)                           nlay !
!     selffrac  - real, factor for temperature interpolation of         !
!                   reference water vapor self-continuum data      nlay !
!     indself   - integer, index of lower reference temperature for     !
!                   the self-continuum interpolation               nlay !
!     forfac    - real, scale factor for w. v. foreign-continuum   nlay !
!     forfrac   - real, factor for temperature interpolation of         !
!                   reference w.v. foreign-continuum data          nlay !
!     indfor    - integer, index of lower reference temperature for     !
!                   the foreign-continuum interpolation            nlay !
!     minorfrac - real, factor for minor gases                     nlay !
!     scaleminor,scaleminorn2                                           !
!               - real, scale factors for minor gases              nlay !
!     indminor  - integer, index of lower reference temperature for     !
!                   minor gases                                    nlay !
!     nlay      - integer, total number of layers                   1   !
!                                                                       !
!  outputs:                                                             !
!     fracs     - real, planck fractions                     ngptlw,nlay!
!     tautot    - real, total optical depth (gas+aerosols)   ngptlw,nlay!
!                                                                       !
!  internal variables:                                                  !
!     ng##      - integer, number of g-values in band ## (##=01-16) 1   !
!     nspa      - integer, for lower atmosphere, the number of ref      !
!                   atmos, each has different relative amounts of the   !
!                   key species for the band                      nbands!
!     nspb      - integer, same but for upper atmosphere          nbands!
!     absa      - real, k-values for lower ref atmospheres (no w.v.     !
!                   self-continuum) (cm**2/molecule)  nspa(##)*5*13*ng##!
!     absb      - real, k-values for high ref atmospheres (all sources) !
!                   (cm**2/molecule)               nspb(##)*5*13:59*ng##!
!     ka_m'mgas'- real, k-values for low ref atmospheres minor species  !
!                   (cm**2/molecule)                          mmn##*ng##!
!     kb_m'mgas'- real, k-values for high ref atmospheres minor species !
!                   (cm**2/molecule)                          mmn##*ng##!
!     selfref   - real, k-values for w.v. self-continuum for ref atmos  !
!                   used below laytrop (cm**2/mol)               10*ng##!
!     forref    - real, k-values for w.v. foreign-continuum for ref atmos
!                   used below/above laytrop (cm**2/mol)          4*ng##!
!                                                                       !
!  ******************************************************************   !

!  ---  inputs:
      integer, intent(in) :: nlay, laytrop

      integer, dimension(nlay), intent(in) :: jp, jt, jt1, indself,     &
     &       indfor, indminor

      real (kind=kind_phys), dimension(nlay), intent(in) :: pavel,      &
     &       coldry, colbrd, fac00, fac01, fac10, fac11, selffac,       &
     &       selffrac, forfac, forfrac, minorfrac, scaleminor,          &
     &       scaleminorn2

      real (kind=kind_phys), dimension(nlay,maxgas), intent(in):: colamt
      real (kind=kind_phys), dimension(nlay,maxxsec),intent(in):: wx

      real (kind=kind_phys), dimension(nbands,nlay), intent(in):: tauaer

      real (kind=kind_phys), dimension(nlay,nrates,2), intent(in) ::    &
     &       rfrate

!  ---  outputs:
      real (kind=kind_phys), dimension(ngptlw,nlay), intent(out) ::     &
     &       fracs, tautot

!  ---  locals
      real (kind=kind_phys), dimension(ngptlw,nlay) :: taug

      integer :: ib, ig, k
!
!===> ...  begin here
!
      call taugb01
      call taugb02
      call taugb03
      call taugb04
      call taugb05
      call taugb06
      call taugb07
      call taugb08
      call taugb09
      call taugb10
      call taugb11
      call taugb12
      call taugb13
      call taugb14
      call taugb15
      call taugb16

!  ---  combine gaseous and aerosol optical depths

      do ig = 1, ngptlw
        ib = ngb(ig)

        do k = 1, nlay
          tautot(ig,k) = taug(ig,k) + tauaer(ib,k)
        enddo
      enddo

! =================
      contains
! =================

!>\ingroup module_radlw_main
!> band 1:  10-350 cm-1 (low key - h2o; low minor - n2);
!!  (high key - h2o; high minor - n2)
! ----------------------------------
      subroutine taugb01
! ..................................

!  ------------------------------------------------------------------  !
!  written by eli j. mlawer, atmospheric & environmental research.     !
!  revised by michael j. iacono, atmospheric & environmental research. !
!                                                                      !
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)             !
!                          (high key - h2o; high minor - n2)           !
!                                                                      !
!  compute the optical depth by interpolating in ln(pressure) and      !
!  temperature.  below laytrop, the water vapor self-continuum and     !
!  foreign continuum is interpolated (in temperature) separately.      !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb01

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &           indm, indmp, ig

      real (kind=kind_phys) :: pp, corradj, scalen2, tauself, taufor,   &
     &       taun2
!
!===> ...  begin here
!
!  ---  minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(1) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(1) + 1
        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)

        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1

        pp = pavel(k)
        scalen2 = colbrd(k) * scaleminorn2(k)
        if (pp < 250.0) then
          corradj = f_one - 0.15 * (250.0-pp) / 154.4
        else
          corradj = f_one
        endif

        do ig = 1, ng01
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) -  forref(ig,indf)))
          taun2   = scalen2 * (ka_mn2(ig,indm) + minorfrac(k)           &
     &            * (ka_mn2(ig,indmp) - ka_mn2(ig,indm)))

          taug(ig,k) = corradj * (colamt(k,1)                           &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor + taun2)

          fracs(ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(1) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(1) + 1
        indf = indfor(k)
        indm = indminor(k)

        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indfp = indf + 1
        indmp = indm + 1

        scalen2 = colbrd(k) * scaleminorn2(k)
        corradj = f_one - 0.15 * (pavel(k) / 95.6)

        do ig = 1, ng01
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))
          taun2  = scalen2 * (kb_mn2(ig,indm) + minorfrac(k)            &
     &           * (kb_mn2(ig,indmp) - kb_mn2(ig,indm)))

          taug(ig,k) = corradj * (colamt(k,1)                           &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))   &
     &           + taufor + taun2)

          fracs(ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb01
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
! ----------------------------------
      subroutine taugb02
! ..................................

!  ------------------------------------------------------------------  !
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)            !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb02

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &           ig

      real (kind=kind_phys) :: corradj, tauself, taufor
!
!===> ...  begin here
!
!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(2) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(2) + 1
        inds = indself(k)
        indf = indfor(k)

        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1

        corradj = f_one - 0.05 * (pavel(k) - 100.0) / 900.0

        do ig = 1, ng02
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns02+ig,k) = corradj * (colamt(k,1)                      &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor)

          fracs(ns02+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(2) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(2) + 1
        indf = indfor(k)

        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indfp = indf + 1

        do ig = 1, ng02
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns02+ig,k) = colamt(k,1)                                 &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))   &
     &           + taufor

          fracs(ns02+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb02
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o);
!!                        (high key - h2o,co2; high minor - n2o)
! ----------------------------------
      subroutine taugb03
! ..................................

!  ------------------------------------------------------------------  !
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)       !
!                           (high key - h2o,co2; high minor - n2o)     !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb03

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmn2o, jmn2op,   &
     &       id001, id011, id101, id111, id201, id211, jpl, jplp,       &
     &       ig, js, js1

      real (kind=kind_phys) ::  absn2o, ratn2o, adjfac, adjcoln2o,      &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mn2o,  specparm_mn2o,  specmult_mn2o,  fmn2o,      &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_planck_b, refrat_m_a, refrat_m_b,   &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      tau_major, tau_major1, tauself, taufor, n2om1, n2om2,       &
     &      p, p4, fk0, fk1, fk2
!
!===> ...  begin here
!
!  --- ...  minor gas mapping levels:
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

      refrat_planck_a = chi_mls(1,9)/chi_mls(2,9)    ! P = 212.725 mb
      refrat_planck_b = chi_mls(1,13)/chi_mls(2,13)  ! P = 95.58   mb
      refrat_m_a      = chi_mls(1,3)/chi_mls(2,3)    ! P = 706.270 mb
      refrat_m_b      = chi_mls(1,13)/chi_mls(2,13)  ! P = 95.58   mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(3) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(3) + js1

        speccomb_mn2o = colamt(k,1) + refrat_m_a*colamt(k,2)
        specparm_mn2o = colamt(k,1) / speccomb_mn2o
        specmult_mn2o = 8.0 * min(specparm_mn2o, oneminus)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = mod(specmult_mn2o, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jmn2op= jmn2o+ 1
        jplp  = jpl  + 1

!  --- ...  in atmospheres where the amount of n2O is too great to be considered
!           a minor species, adjust the column amount of n2O by an empirical factor
!           to obtain the proper contribution.

        p = coldry(k) * chi_mls(4,jp(k)+1)
        ratn2o = colamt(k,4) / p
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o - 0.5)**0.65
          adjcoln2o = adjfac * p
        else
          adjcoln2o = colamt(k,4)
        endif

        if (specparm < 0.125) then
          p = fs - f_one
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        else if (specparm > 0.875) then
          p = -fs
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk0 = f_one - fs
          fk1 = fs
          fk2 = f_zero
          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk0*fac00(k)
        fac100 = fk1*fac00(k)
        fac200 = fk2*fac00(k)
        fac010 = fk0*fac10(k)
        fac110 = fk1*fac10(k)
        fac210 = fk2*fac10(k)

        if (specparm1 < 0.125) then
          p = fs1 - f_one
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p = -fs1
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk0 = f_one - fs1
          fk1 = fs1
          fk2 = f_zero
          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk0*fac01(k)
        fac101 = fk1*fac01(k)
        fac201 = fk2*fac01(k)
        fac011 = fk0*fac11(k)
        fac111 = fk1*fac11(k)
        fac211 = fk2*fac11(k)

        do ig = 1, ng03
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          n2om1   = ka_mn2o(ig,jmn2o,indm) + fmn2o                      &
     &            * (ka_mn2o(ig,jmn2op,indm) - ka_mn2o(ig,jmn2o,indm))
          n2om2   = ka_mn2o(ig,jmn2o,indmp) + fmn2o                     &
     &            * (ka_mn2o(ig,jmn2op,indmp) - ka_mn2o(ig,jmn2o,indmp))
          absn2o  = n2om1 + minorfrac(k) * (n2om2 - n2om1)

          tau_major = speccomb                                          &
     &              * (fac000*absa(ig,id000) + fac010*absa(ig,id010)    &
     &              +  fac100*absa(ig,id100) + fac110*absa(ig,id110)    &
     &              +  fac200*absa(ig,id200) + fac210*absa(ig,id210))

          tau_major1 = speccomb1                                        &
     &              * (fac001*absa(ig,id001) + fac011*absa(ig,id011)    &
     &              +  fac101*absa(ig,id101) + fac111*absa(ig,id111)    &
     &              +  fac201*absa(ig,id201) + fac211*absa(ig,id211))

          taug(ns03+ig,k) = tau_major + tau_major1                      &
     &                    + tauself + taufor + adjcoln2o*absn2o

          fracs(ns03+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo     ! end do_k_loop
      enddo   ! end do_ig_loop

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 4.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-13)*5 + (jt(k)-1)) * nspb(3) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 4.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(3) + js1

        speccomb_mn2o = colamt(k,1) + refrat_m_b*colamt(k,2)
        specparm_mn2o = colamt(k,1) / speccomb_mn2o
        specmult_mn2o = 4.0 * min(specparm_mn2o, oneminus)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = mod(specmult_mn2o, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_b*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 4.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        indf = indfor(k)
        indm = indminor(k)
        indfp = indf + 1
        indmp = indm + 1
        jmn2op= jmn2o+ 1
        jplp  = jpl  + 1

        id000 = ind0
        id010 = ind0 + 5
        id100 = ind0 + 1
        id110 = ind0 + 6
        id001 = ind1
        id011 = ind1 + 5
        id101 = ind1 + 1
        id111 = ind1 + 6

!  --- ...  in atmospheres where the amount of n2o is too great to be considered
!           a minor species, adjust the column amount of N2O by an empirical factor
!           to obtain the proper contribution.

        p = coldry(k) * chi_mls(4,jp(k)+1)
        ratn2o = colamt(k,4) / p
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o - 0.5)**0.65
          adjcoln2o = adjfac * p
        else
          adjcoln2o = colamt(k,4)
        endif

        fk0 = f_one - fs
        fk1 = fs
        fac000 = fk0*fac00(k)
        fac010 = fk0*fac10(k)
        fac100 = fk1*fac00(k)
        fac110 = fk1*fac10(k)

        fk0 = f_one - fs1
        fk1 = fs1
        fac001 = fk0*fac01(k)
        fac011 = fk0*fac11(k)
        fac101 = fk1*fac01(k)
        fac111 = fk1*fac11(k)

        do ig = 1, ng03
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))
          n2om1  = kb_mn2o(ig,jmn2o,indm) + fmn2o                       &
     &           * (kb_mn2o(ig,jmn2op,indm) - kb_mn2o(ig,jmn2o,indm))
          n2om2  = kb_mn2o(ig,jmn2o,indmp) + fmn2o                      &
     &           * (kb_mn2o(ig,jmn2op,indmp) - kb_mn2o(ig,jmn2o,indmp))
          absn2o = n2om1 + minorfrac(k) * (n2om2 - n2om1)

          tau_major = speccomb                                          &
     &              * (fac000*absb(ig,id000) + fac010*absb(ig,id010)    &
     &              +  fac100*absb(ig,id100) + fac110*absb(ig,id110))

          tau_major1 = speccomb1                                        &
     &              * (fac001*absb(ig,id001) + fac011*absb(ig,id011)    &
     &              +  fac101*absb(ig,id101) + fac111*absb(ig,id111))

          taug(ns03+ig,k) = tau_major + tau_major1                      &
     &                    + taufor + adjcoln2o*absn2o

          fracs(ns03+ig,k) = fracrefb(ig,jpl) + fpl                     &
     &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
        enddo
      enddo

! ..................................
      end subroutine taugb03
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
! ----------------------------------
      subroutine taugb04
! ..................................

!  ------------------------------------------------------------------  !
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)     !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb04

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, jpl, jplp,    &
     &       id000, id010, id100, id110, id200, id210, ig, js, js1,     &
     &       id001, id011, id101, id111, id201, id211

      real (kind=kind_phys) :: tauself, taufor, p, p4, fk0, fk1, fk2,   &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      refrat_planck_a, refrat_planck_b, tau_major, tau_major1
!
!===> ...  begin here
!
      refrat_planck_a = chi_mls(1,11)/chi_mls(2,11)     ! P = 142.5940 mb
      refrat_planck_b = chi_mls(3,13)/chi_mls(2,13)     ! P = 95.58350 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(4) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = ( jp(k)*5 + (jt1(k)-1)) * nspa(4) + js1

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, 1.0)

        inds = indself(k)
        indf = indfor(k)
        indsp = inds + 1
        indfp = indf + 1
        jplp  = jpl  + 1

        if (specparm < 0.125) then
          p = fs - f_one
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11
        elseif (specparm > 0.875) then
          p = -fs
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8
        else
          fk0 = f_one - fs
          fk1 = fs
          fk2 = f_zero
          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0
        endif

        fac000 = fk0*fac00(k)
        fac100 = fk1*fac00(k)
        fac200 = fk2*fac00(k)
        fac010 = fk0*fac10(k)
        fac110 = fk1*fac10(k)
        fac210 = fk2*fac10(k)

        if (specparm1 < 0.125) then
          p = fs1 - f_one
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm1 > 0.875) then
          p = -fs1
          p4 = p**4
          fk0 = p4
          fk1 = f_one - p - 2.0*p4
          fk2 = p + p4
          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk0 = f_one - fs1
          fk1 = fs1
          fk2 = f_zero
          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac001 = fk0*fac01(k)
        fac101 = fk1*fac01(k)
        fac201 = fk2*fac01(k)
        fac011 = fk0*fac11(k)
        fac111 = fk1*fac11(k)
        fac211 = fk2*fac11(k)

        do ig = 1, ng04
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          tau_major = speccomb                                          &
     &              * (fac000*absa(ig,id000) + fac010*absa(ig,id010)    &
     &              +  fac100*absa(ig,id100) + fac110*absa(ig,id110)    &
     &              +  fac200*absa(ig,id200) + fac210*absa(ig,id210))

          tau_major1 = speccomb1                                        &
     &              * (fac001*absa(ig,id001) + fac011*absa(ig,id011)    &
     &              +  fac101*absa(ig,id101) + fac111*absa(ig,id111)    &
     &              +  fac201*absa(ig,id201) + fac211*absa(ig,id211))

          taug(ns04+ig,k) = tau_major + tau_major1 + tauself + taufor

          fracs(ns04+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo     ! end do_k_loop
      enddo   ! end do_ig_loop

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        speccomb = colamt(k,3) + rfrate(k,6,1)*colamt(k,2)
        specparm = colamt(k,3) / speccomb
        specmult = 4.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-13)*5 + (jt(k)-1)) * nspb(4) + js

        speccomb1 = colamt(k,3) + rfrate(k,6,2)*colamt(k,2)
        specparm1 = colamt(k,3) / speccomb1
        specmult1 = 4.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(4) + js1

        speccomb_planck = colamt(k,3) + refrat_planck_b*colamt(k,2)
        specparm_planck = colamt(k,3) / speccomb_planck
        specmult_planck = 4.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)
        jplp = jpl + 1

        id000 = ind0
        id010 = ind0 + 5
        id100 = ind0 + 1
        id110 = ind0 + 6
        id001 = ind1
        id011 = ind1 + 5
        id101 = ind1 + 1
        id111 = ind1 + 6

        fk0 = f_one - fs
        fk1 = fs
        fac000 = fk0*fac00(k)
        fac010 = fk0*fac10(k)
        fac100 = fk1*fac00(k)
        fac110 = fk1*fac10(k)

        fk0 = f_one - fs1
        fk1 = fs1
        fac001 = fk0*fac01(k)
        fac011 = fk0*fac11(k)
        fac101 = fk1*fac01(k)
        fac111 = fk1*fac11(k)

        do ig = 1, ng04
          tau_major =  speccomb                                         &
     &              * (fac000*absb(ig,id000) + fac010*absb(ig,id010)    &
     &              +  fac100*absb(ig,id100) + fac110*absb(ig,id110))
          tau_major1 = speccomb1                                        &
     &              * (fac001*absb(ig,id001) + fac011*absb(ig,id011)    &
     &              +  fac101*absb(ig,id101) + fac111*absb(ig,id111))

          taug(ns04+ig,k) =  tau_major + tau_major1

          fracs(ns04+ig,k) = fracrefb(ig,jpl) + fpl                     &
     &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
        enddo

!  --- ...  empirical modification to code to improve stratospheric cooling rates
!           for co2. revised to apply weighting for g-point reduction in this band.

        taug(ns04+ 8,k) = taug(ns04+ 8,k) * 0.92
        taug(ns04+ 9,k) = taug(ns04+ 9,k) * 0.88
        taug(ns04+10,k) = taug(ns04+10,k) * 1.07
        taug(ns04+11,k) = taug(ns04+11,k) * 1.1
        taug(ns04+12,k) = taug(ns04+12,k) * 0.99
        taug(ns04+13,k) = taug(ns04+13,k) * 0.88
        taug(ns04+14,k) = taug(ns04+14,k) * 0.943
      enddo

! ..................................
      end subroutine taugb04
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!!                       (high key - o3,co2)
! ----------------------------------
      subroutine taugb05
! ..................................

!  ------------------------------------------------------------------  !
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)  !
!                           (high key - o3,co2)                        !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb05

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmo3, jmo3p,     &
     &       id001, id011, id101, id111, id201, id211, jpl, jplp,       &
     &       ig, js, js1

      real (kind=kind_phys)  :: tauself, taufor, o3m1, o3m2, abso3,     &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mo3,   specparm_mo3,   specmult_mo3,   fmo3,       &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_planck_b, refrat_m_a,               &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - o3, p = 317.34 mbar, t = 240.77 k
!     lower - ccl4

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower/upper atmosphere.

      refrat_planck_a = chi_mls(1,5)/chi_mls(2,5)      ! P = 473.420 mb
      refrat_planck_b = chi_mls(3,43)/chi_mls(2,43)    ! P = 0.2369  mb
      refrat_m_a = chi_mls(1,7)/chi_mls(2,7)           ! P = 317.348 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(5) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(5) + js1

        speccomb_mo3 = colamt(k,1) + refrat_m_a*colamt(k,2)
        specparm_mo3 = colamt(k,1) / speccomb_mo3
        specmult_mo3 = 8.0 * min(specparm_mo3, oneminus)
        jmo3 = 1 + int(specmult_mo3)
        fmo3 = mod(specmult_mo3, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmo3p = jmo3 + 1

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm > 0.875 .and. specparm1 > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng05
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          o3m1    = ka_mo3(ig,jmo3,indm) + fmo3                         &
     &            * (ka_mo3(ig,jmo3p,indm) -  ka_mo3(ig,jmo3,indm))
          o3m2    = ka_mo3(ig,jmo3,indmp) + fmo3                        &
     &            * (ka_mo3(ig,jmo3p,indmp) - ka_mo3(ig,jmo3,indmp))
          abso3   = o3m1 + minorfrac(k)*(o3m2 - o3m1)

          taug(ns05+ig,k) = speccomb                                    &
     &            * (fac000*absa(ig,id000) + fac010*absa(ig,id010)      &
     &            +  fac100*absa(ig,id100) + fac110*absa(ig,id110)      &
     &            +  fac200*absa(ig,id200) + fac210*absa(ig,id210))     &
     &            +     speccomb1                                       &
     &            * (fac001*absa(ig,id001) + fac011*absa(ig,id011)      &
     &            +  fac101*absa(ig,id101) + fac111*absa(ig,id111)      &
     &            +  fac201*absa(ig,id201) + fac211*absa(ig,id211))     &
     &            + tauself + taufor+abso3*colamt(k,3)+wx(k,1)*ccl4(ig)

          fracs(ns05+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        speccomb = colamt(k,3) + rfrate(k,6,1)*colamt(k,2)
        specparm = colamt(k,3) / speccomb
        specmult = 4.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-13)*5 + (jt(k)-1)) * nspb(5) + js

        speccomb1 = colamt(k,3) + rfrate(k,6,2)*colamt(k,2)
        specparm1 = colamt(k,3) / speccomb1
        specmult1 = 4.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(5) + js1

        speccomb_planck = colamt(k,3) + refrat_planck_b*colamt(k,2)
        specparm_planck = colamt(k,3) / speccomb_planck
        specmult_planck = 4.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)
        jplp= jpl + 1

        id000 = ind0
        id010 = ind0 + 5
        id100 = ind0 + 1
        id110 = ind0 + 6
        id001 = ind1
        id011 = ind1 + 5
        id101 = ind1 + 1
        id111 = ind1 + 6

        fk00 = f_one - fs
        fk10 = fs

        fk01 = f_one - fs1
        fk11 = fs1

        fac000 = fk00 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac100 = fk10 * fac00(k)
        fac110 = fk10 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac101 = fk11 * fac01(k)
        fac111 = fk11 * fac11(k)

        do ig = 1, ng05
          taug(ns05+ig,k) = speccomb                                    &
     &                * (fac000*absb(ig,id000) + fac010*absb(ig,id010)  &
     &                +  fac100*absb(ig,id100) + fac110*absb(ig,id110)) &
     &                +     speccomb1                                   &
     &                * (fac001*absb(ig,id001) + fac011*absb(ig,id011)  &
     &                +  fac101*absb(ig,id101) + fac111*absb(ig,id111)) &
     &                + wx(k,1) * ccl4(ig)

          fracs(ns05+ig,k) = fracrefb(ig,jpl) + fpl                     &
     &                     * (fracrefb(ig,jplp) - fracrefb(ig,jpl))
        enddo
      enddo

! ..................................
      end subroutine taugb05
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!!                       (high key - none; high minor - cfc11, cfc12)
! ----------------------------------
      subroutine taugb06
! ..................................

!  ------------------------------------------------------------------  !
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)           !
!                           (high key - none; high minor - cfc11, cfc12)
!  ------------------------------------------------------------------  !

      use module_radlw_kgb06

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: ratco2, adjfac, adjcolco2, tauself,      &
     &      taufor, absco2, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level:
!     lower - co2, p = 706.2720 mb, t = 294.2 k
!     upper - cfc11, cfc12

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(6) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(6) + 1

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        ind0p = ind0 + 1
        ind1p = ind1 + 1

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.77
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        do ig = 1, ng06
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          absco2  = ka_mco2(ig,indm) + minorfrac(k)                     &
     &            * (ka_mco2(ig,indmp) - ka_mco2(ig,indm))

          taug(ns06+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            +  tauself + taufor + adjcolco2*absco2                &
     &            +  wx(k,2)*cfc11adj(ig) + wx(k,3)*cfc12(ig)

          fracs(ns06+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop
!           nothing important goes on above laytrop in this band.

      do k = laytrop+1, nlay
        do ig = 1, ng06
          taug(ns06+ig,k) = wx(k,2)*cfc11adj(ig) + wx(k,3)*cfc12(ig)

          fracs(ns06+ig,k) = fracrefa(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb06
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!!                        (high key - o3; high minor - co2)
! ----------------------------------
      subroutine taugb07
! ..................................

!  ------------------------------------------------------------------  !
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)       !
!                            (high key - o3; high minor - co2)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb07

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, indm, indmp,     &
     &       id001, id011, id101, id111, id201, id211, jmco2, jmco2p,   &
     &       jpl, jplp, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, co2m1, co2m2, absco2,   &
     &      speccomb,       specparm,       specmult,       fs,         &
     &      speccomb1,      specparm1,      specmult1,      fs1,        &
     &      speccomb_mco2,  specparm_mco2,  specmult_mco2,  fmco2,      &
     &      speccomb_planck,specparm_planck,specmult_planck,fpl,        &
     &      refrat_planck_a, refrat_m_a, ratco2, adjfac, adjcolco2,     &
     &      fac000, fac100, fac200, fac010, fac110, fac210,             &
     &      fac001, fac101, fac201, fac011, fac111, fac211,             &
     &      p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - co2, p = 706.2620 mbar, t= 278.94 k
!     upper - co2, p = 12.9350 mbar, t = 234.01 k

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower atmosphere.

      refrat_planck_a = chi_mls(1,3)/chi_mls(3,3)     ! P = 706.2620 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(3,3)          ! P = 706.2720 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,2,1)*colamt(k,3)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(7) + js

        speccomb1 = colamt(k,1) + rfrate(k,2,2)*colamt(k,3)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(7) + js1

        speccomb_mco2 = colamt(k,1) + refrat_m_a*colamt(k,3)
        specparm_mco2 = colamt(k,1) / speccomb_mco2
        specmult_mco2 = 8.0 * min(specparm_mco2, oneminus)
        jmco2 = 1 + int(specmult_mco2)
        fmco2 = mod(specmult_mco2, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,3)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmco2p= jmco2+ 1
        ind0p = ind0 + 1
        ind1p = ind1 + 1

!  --- ...  in atmospheres where the amount of CO2 is too great to be considered
!           a minor species, adjust the column amount of CO2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 3.0 + (ratco2-3.0)**0.79
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm > 0.875 .and. specparm1 > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng07
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          co2m1   = ka_mco2(ig,jmco2,indm) + fmco2                      &
     &            * (ka_mco2(ig,jmco2p,indm) - ka_mco2(ig,jmco2,indm))
          co2m2   = ka_mco2(ig,jmco2,indmp) + fmco2                     &
     &            * (ka_mco2(ig,jmco2p,indmp) - ka_mco2(ig,jmco2,indmp))
          absco2  = co2m1 + minorfrac(k) * (co2m2 - co2m1)

          taug(ns07+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor + adjcolco2*absco2

          fracs(ns07+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

      do k = laytrop+1, nlay
        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.79
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(7) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(7) + 1

        indm = indminor(k)
        indmp = indm + 1
        ind0p = ind0 + 1
        ind1p = ind1 + 1

        do ig = 1, ng07
          absco2 = kb_mco2(ig,indm) + minorfrac(k)                      &
     &           * (kb_mco2(ig,indmp) - kb_mco2(ig,indm))

          taug(ns07+ig,k) = colamt(k,3)                                 &
     &            * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)   &
     &            +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))  &
     &            + adjcolco2 * absco2

          fracs(ns07+ig,k) = fracrefb(ig)
        enddo

!  --- ...  empirical modification to code to improve stratospheric cooling rates
!           for o3.  revised to apply weighting for g-point reduction in this band.

        taug(ns07+ 6,k) = taug(ns07+ 6,k) * 0.92
        taug(ns07+ 7,k) = taug(ns07+ 7,k) * 0.88
        taug(ns07+ 8,k) = taug(ns07+ 8,k) * 1.07
        taug(ns07+ 9,k) = taug(ns07+ 9,k) * 1.1
        taug(ns07+10,k) = taug(ns07+10,k) * 0.99
        taug(ns07+11,k) = taug(ns07+11,k) * 0.855
      enddo

! ..................................
      end subroutine taugb07
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!!                         (high key - o3; high minor - co2, n2o)
! ----------------------------------
      subroutine taugb08
! ..................................

!  ------------------------------------------------------------------  !
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)  !
!                             (high key - o3; high minor - co2, n2o)   !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb08

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: tauself, taufor, absco2, abso3, absn2o,  &
     &      ratco2, adjfac, adjcolco2, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level:
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - o3,  p = 317.348 mb, t = 240.77 k
!     lower - n2o, p = 706.2720 mb, t= 278.94 k
!     lower - cfc12,cfc11
!     upper - co2, p = 35.1632 mb, t = 223.28 k
!     upper - n2o, p = 8.716e-2 mb, t = 226.03 k

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(8) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(8) + 1

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.65
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        do ig = 1, ng08
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          absco2  = (ka_mco2(ig,indm) + minorfrac(k)                    &
     &            * (ka_mco2(ig,indmp) - ka_mco2(ig,indm)))
          abso3   = (ka_mo3(ig,indm) + minorfrac(k)                     &
     &            * (ka_mo3(ig,indmp) - ka_mo3(ig,indm)))
          absn2o  = (ka_mn2o(ig,indm) + minorfrac(k)                    &
     &            * (ka_mn2o(ig,indmp) - ka_mn2o(ig,indm)))

          taug(ns08+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself+taufor + adjcolco2*absco2                   &
     &            + colamt(k,3)*abso3 + colamt(k,4)*absn2o              &
     &            + wx(k,3)*cfc12(ig) + wx(k,4)*cfc22adj(ig)

          fracs(ns08+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(8) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(8) + 1

        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indmp = indm + 1

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(2,jp(k)+1)
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.65
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        do ig = 1, ng08
          absco2 = (kb_mco2(ig,indm) + minorfrac(k)                     &
     &           * (kb_mco2(ig,indmp) - kb_mco2(ig,indm)))
          absn2o = (kb_mn2o(ig,indm) + minorfrac(k)                     &
     &           * (kb_mn2o(ig,indmp) - kb_mn2o(ig,indm)))

          taug(ns08+ig,k) = colamt(k,3)                                 &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))   &
     &           + adjcolco2*absco2 + colamt(k,4)*absn2o                &
     &           + wx(k,3)*cfc12(ig) + wx(k,4)*cfc22adj(ig)

          fracs(ns08+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb08
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!!                         (high key - ch4; high minor - n2o)
! ----------------------------------
      subroutine taugb09
! ..................................

!  ------------------------------------------------------------------  !
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)     !
!                             (high key - ch4; high minor - n2o)       !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb09

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, indm, indmp,     &
     &       id001, id011, id101, id111, id201, id211, jmn2o, jmn2op,   &
     &       jpl, jplp, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, n2om1, n2om2, absn2o,   &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mn2o,  specparm_mn2o,  specmult_mn2o,  fmn2o,     &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, ratn2o, adjfac, adjcoln2o,    &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - n2o, p = 706.272 mbar, t = 278.94 k
!     upper - n2o, p = 95.58 mbar, t = 215.7 k

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower/upper atmosphere.

      refrat_planck_a = chi_mls(1,9)/chi_mls(6,9)       ! P = 212 mb
      refrat_m_a = chi_mls(1,3)/chi_mls(6,3)            ! P = 706.272 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,4,1)*colamt(k,5)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(9) + js

        speccomb1 = colamt(k,1) + rfrate(k,4,2)*colamt(k,5)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(9) + js1

        speccomb_mn2o = colamt(k,1) + refrat_m_a*colamt(k,5)
        specparm_mn2o = colamt(k,1) / speccomb_mn2o
        specmult_mn2o = 8.0 * min(specparm_mn2o, oneminus)
        jmn2o = 1 + int(specmult_mn2o)
        fmn2o = mod(specmult_mn2o, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,5)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmn2op= jmn2o+ 1

!  --- ...  in atmospheres where the amount of n2o is too great to be considered
!           a minor species, adjust the column amount of n2o by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(4,jp(k)+1)
        ratn2o = colamt(k,4) / temp
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o-0.5)**0.65
          adjcoln2o = adjfac * temp
        else
          adjcoln2o = colamt(k,4)
        endif

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11

        elseif (specparm > 0.875 .and. specparm1 > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng09
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          n2om1   = ka_mn2o(ig,jmn2o,indm) + fmn2o                      &
     &            * (ka_mn2o(ig,jmn2op,indm) - ka_mn2o(ig,jmn2o,indm))
          n2om2   = ka_mn2o(ig,jmn2o,indmp) + fmn2o                     &
     &            * (ka_mn2o(ig,jmn2op,indmp) - ka_mn2o(ig,jmn2o,indmp))
          absn2o  = n2om1 + minorfrac(k) * (n2om2 - n2om1)

          taug(ns09+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor + adjcoln2o*absn2o

          fracs(ns09+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(9) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(9) + 1

        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indmp = indm + 1

!  --- ...  in atmospheres where the amount of n2o is too great to be considered
!           a minor species, adjust the column amount of n2o by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * chi_mls(4,jp(k)+1)
        ratn2o = colamt(k,4) / temp
        if (ratn2o > 1.5) then
          adjfac = 0.5 + (ratn2o - 0.5)**0.65
          adjcoln2o = adjfac * temp
        else
          adjcoln2o = colamt(k,4)
        endif

        do ig = 1, ng09
          absn2o = kb_mn2o(ig,indm) + minorfrac(k)                      &
     &           * (kb_mn2o(ig,indmp) - kb_mn2o(ig,indm))

          taug(ns09+ig,k) = colamt(k,5)                                 &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))   &
     &           + adjcoln2o*absn2o

          fracs(ns09+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb09
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
! ----------------------------------
      subroutine taugb10
! ..................................

!  ------------------------------------------------------------------  !
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb10

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       ig

      real (kind=kind_phys) :: tauself, taufor
!
!===> ...  begin here
!
!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(10) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(10) + 1

        inds = indself(k)
        indf = indfor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1

        do ig = 1, ng10
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns10+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor

          fracs(ns10+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(10) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(10) + 1

        indf = indfor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indfp = indf + 1

        do ig = 1, ng10
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns10+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)   &
     &            +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))  &
     &            + taufor

          fracs(ns10+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb10
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!!                          (high key - h2o; high minor - o2)
! ----------------------------------
      subroutine taugb11
! ..................................

!  ------------------------------------------------------------------  !
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)             !
!                              (high key - h2o; high minor - o2)       !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb11

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       indm, indmp, ig

      real (kind=kind_phys) :: scaleo2, tauself, taufor, tauo2
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - o2, p = 706.2720 mbar, t = 278.94 k
!     upper - o2, p = 4.758820 mbarm t = 250.85 k

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(11) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(11) + 1

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1

        scaleo2 = colamt(k,6) * scaleminor(k)

        do ig = 1, ng11
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          tauo2   = scaleo2 * (ka_mo2(ig,indm) + minorfrac(k)           &
     &            * (ka_mo2(ig,indmp) - ka_mo2(ig,indm)))

          taug(ns11+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor + tauo2

          fracs(ns11+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(11) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(11) + 1

        indf = indfor(k)
        indm = indminor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indfp = indf + 1
        indmp = indm + 1

        scaleo2 = colamt(k,6) * scaleminor(k)

        do ig = 1, ng11
          taufor = forfac(k) * (forref(ig,indf) + forfrac(k)            &
     &           * (forref(ig,indfp) - forref(ig,indf)))
          tauo2  = scaleo2 * (kb_mo2(ig,indm) + minorfrac(k)            &
     &           * (kb_mo2(ig,indmp) - kb_mo2(ig,indm)))

          taug(ns11+ig,k) = colamt(k,1)                                 &
     &            * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)   &
     &            +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))  &
     &            + taufor + tauo2

          fracs(ns11+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb11
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
! ----------------------------------
      subroutine taugb12
! ..................................

!  ------------------------------------------------------------------  !
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)         !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb12

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, jpl, jplp,    &
     &       id000, id010, id100, id110, id200, id210, ig, js, js1,     &
     &       id001, id011, id101, id111, id201, id211

      real (kind=kind_phys) :: tauself, taufor, refrat_planck_a,        &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21
!
!===> ...  begin here
!
!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower/upper atmosphere.

      refrat_planck_a = chi_mls(1,10)/chi_mls(2,10)      ! P =   174.164 mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,1,1)*colamt(k,2)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(12) + js

        speccomb1 = colamt(k,1) + rfrate(k,1,2)*colamt(k,2)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(12) + js1

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,1) / speccomb_planck
        if (specparm_planck >= oneminus) specparm_planck=oneminus
        specmult_planck = 8.0 * specparm_planck
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indsp = inds + 1
        indfp = indf + 1
        jplp  = jpl  + 1

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm > 0.875 .and. specparm1 > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng12
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns12+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor

          fracs(ns12+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     *(fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        do ig = 1, ng12
          taug(ns12+ig,k) = f_zero
          fracs(ns12+ig,k) = f_zero
        enddo
      enddo

! ..................................
      end subroutine taugb12
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 13:  2080-2250 cm-1 (low key-h2o,n2o; high minor-o3 minor)
! ----------------------------------
      subroutine taugb13
! ..................................

!  ------------------------------------------------------------------  !
!     band 13:  2080-2250 cm-1 (low key-h2o,n2o; high minor-o3 minor)  !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb13

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jmco2, jpl,      &
     &       id001, id011, id101, id111, id201, id211, jmco2p, jplp,    &
     &       jmco, jmcop, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, co2m1, co2m2, absco2,   &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mco2,  specparm_mco2,  specmult_mco2,  fmco2,     &
     &       speccomb_mco,   specparm_mco,   specmult_mco,   fmco,      &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, refrat_m_a3, ratco2,          &
     &       adjfac, adjcolco2, com1, com2, absco, abso3,               &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21, temp
!
!===> ...  begin here
!
!  --- ...  minor gas mapping levels :
!     lower - co2, p = 1053.63 mb, t = 294.2 k
!     lower - co, p = 706 mb, t = 278.94 k
!     upper - o3, p = 95.5835 mb, t = 215.7 k

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower/upper atmosphere.

      refrat_planck_a = chi_mls(1,5)/chi_mls(4,5)        ! P = 473.420 mb (Level 5)
      refrat_m_a = chi_mls(1,1)/chi_mls(4,1)             ! P = 1053. (Level 1)
      refrat_m_a3 = chi_mls(1,3)/chi_mls(4,3)            ! P = 706. (Level 3)

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,3,1)*colamt(k,4)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(13) + js

        speccomb1 = colamt(k,1) + rfrate(k,3,2)*colamt(k,4)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(13) + js1

        speccomb_mco2 = colamt(k,1) + refrat_m_a*colamt(k,4)
        specparm_mco2 = colamt(k,1) / speccomb_mco2
        specmult_mco2 = 8.0 * min(specparm_mco2, oneminus)
        jmco2 = 1 + int(specmult_mco2)
        fmco2 = mod(specmult_mco2, f_one)

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        speccomb_mco = colamt(k,1) + refrat_m_a3*colamt(k,4)
        specparm_mco = colamt(k,1) / speccomb_mco
        specmult_mco = 8.0 * min(specparm_mco, oneminus)
        jmco = 1 + int(specmult_mco)
        fmco = mod(specmult_mco, f_one)

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,4)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmco2p= jmco2+ 1
        jmcop = jmco + 1

!  --- ...  in atmospheres where the amount of co2 is too great to be considered
!           a minor species, adjust the column amount of co2 by an empirical factor
!           to obtain the proper contribution.

        temp   = coldry(k) * 3.55e-4
        ratco2 = colamt(k,2) / temp
        if (ratco2 > 3.0) then
          adjfac = 2.0 + (ratco2-2.0)**0.68
          adjcolco2 = adjfac * temp
        else
          adjcolco2 = colamt(k,2)
        endif

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm > 0.875 .and. specparm1 > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng13
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          co2m1   = ka_mco2(ig,jmco2,indm) + fmco2                      &
     &            * (ka_mco2(ig,jmco2p,indm) - ka_mco2(ig,jmco2,indm))
          co2m2   = ka_mco2(ig,jmco2,indmp) + fmco2                     &
     &            * (ka_mco2(ig,jmco2p,indmp) - ka_mco2(ig,jmco2,indmp))
          absco2  = co2m1 + minorfrac(k) * (co2m2 - co2m1)
          com1    = ka_mco(ig,jmco,indm) + fmco                         &
     &            * (ka_mco(ig,jmcop,indm) - ka_mco(ig,jmco,indm))
          com2    = ka_mco(ig,jmco,indmp) + fmco                        &
     &            * (ka_mco(ig,jmcop,indmp) - ka_mco(ig,jmco,indmp))
          absco   = com1 + minorfrac(k) * (com2 - com1)

          taug(ns13+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor + adjcolco2*absco2             &
     &                + colamt(k,7)*absco

          fracs(ns13+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        indm = indminor(k)
        indmp = indm + 1

        do ig = 1, ng13
          abso3 = kb_mo3(ig,indm) + minorfrac(k)                        &
     &          * (kb_mo3(ig,indmp) - kb_mo3(ig,indm))

          taug(ns13+ig,k) = colamt(k,3)*abso3

          fracs(ns13+ig,k) =  fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb13
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 14:  2250-2380 cm-1 (low - co2; high - co2)
! ----------------------------------
      subroutine taugb14
! ..................................

!  ------------------------------------------------------------------  !
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)                 !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb14

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       ig

      real (kind=kind_phys) :: tauself, taufor
!
!===> ...  begin here
!
!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        ind0 = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(14) + 1
        ind1 = ( jp(k)   *5 + (jt1(k)-1)) * nspa(14) + 1

        inds = indself(k)
        indf = indfor(k)
        ind0p = ind0 + 1
        ind1p = ind1 + 1
        indsp = inds + 1
        indfp = indf + 1

        do ig = 1, ng14
          tauself = selffac(k) * (selfref(ig,inds) + selffrac(k)        &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns14+ig,k) = colamt(k,2)                                 &
     &            * (fac00(k)*absa(ig,ind0) + fac10(k)*absa(ig,ind0p)   &
     &            +  fac01(k)*absa(ig,ind1) + fac11(k)*absa(ig,ind1p))  &
     &            + tauself + taufor

          fracs(ns14+ig,k) = fracrefa(ig)
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(14) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(14) + 1

        ind0p = ind0 + 1
        ind1p = ind1 + 1

        do ig = 1, ng14
          taug(ns14+ig,k) = colamt(k,2)                                 &
     &             * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)  &
     &             +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))

          fracs(ns14+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb14
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!!                          (high - nothing)
! ----------------------------------
      subroutine taugb15
! ..................................

!  ------------------------------------------------------------------  !
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)         !
!                              (high - nothing)                        !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb15

!  ---  locals:
      integer :: k, ind0, ind1, inds, indsp, indf, indfp, indm, indmp,  &
     &       id000, id010, id100, id110, id200, id210, jpl, jplp,       &
     &       id001, id011, id101, id111, id201, id211, jmn2, jmn2p,     &
     &       ig, js, js1

      real (kind=kind_phys) :: scalen2, tauself, taufor,                &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_mn2,   specparm_mn2,   specmult_mn2,   fmn2,      &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       refrat_planck_a, refrat_m_a, n2m1, n2m2, taun2,            &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21
!
!===> ...  begin here
!
!  --- ...  minor gas mapping level :
!     lower - nitrogen continuum, P = 1053., T = 294.

!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower atmosphere.

      refrat_planck_a = chi_mls(4,1)/chi_mls(2,1)      ! P = 1053. mb (Level 1)
      refrat_m_a = chi_mls(4,1)/chi_mls(2,1)           ! P = 1053. mb

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,4) + rfrate(k,5,1)*colamt(k,2)
        specparm = colamt(k,4) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(15) + js

        speccomb1 = colamt(k,4) + rfrate(k,5,2)*colamt(k,2)
        specparm1 = colamt(k,4) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(15) + js1

        speccomb_mn2 = colamt(k,4) + refrat_m_a*colamt(k,2)
        specparm_mn2 = colamt(k,4) / speccomb_mn2
        specmult_mn2 = 8.0 * min(specparm_mn2, oneminus)
        jmn2 = 1 + int(specmult_mn2)
        fmn2 = mod(specmult_mn2, f_one)

        speccomb_planck = colamt(k,4) + refrat_planck_a*colamt(k,2)
        specparm_planck = colamt(k,4) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        scalen2 = colbrd(k) * scaleminor(k)

        inds = indself(k)
        indf = indfor(k)
        indm = indminor(k)
        indsp = inds + 1
        indfp = indf + 1
        indmp = indm + 1
        jplp  = jpl  + 1
        jmn2p = jmn2 + 1


        if (specparm < 0.125 .and. specparm1 < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm > 0.875 .and. specparm1 > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng15
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))
          n2m1    = ka_mn2(ig,jmn2,indm) + fmn2                         &
     &            * (ka_mn2(ig,jmn2p,indm) - ka_mn2(ig,jmn2,indm))
          n2m2    = ka_mn2(ig,jmn2,indmp) + fmn2                        &
     &            * (ka_mn2(ig,jmn2p,indmp) - ka_mn2(ig,jmn2,indmp))
          taun2   = scalen2 * (n2m1 + minorfrac(k) * (n2m2 - n2m1))

          taug(ns15+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor + taun2

          fracs(ns15+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        do ig = 1, ng15
          taug(ns15+ig,k) = f_zero

          fracs(ns15+ig,k) = f_zero
        enddo
      enddo

! ..................................
      end subroutine taugb15
! ----------------------------------

!>\ingroup module_radlw_main
!> Band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
! ----------------------------------
      subroutine taugb16
! ..................................

!  ------------------------------------------------------------------  !
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)      !
!  ------------------------------------------------------------------  !

      use module_radlw_kgb16

!  ---  locals:
      integer :: k, ind0, ind0p, ind1, ind1p, inds, indsp, indf, indfp, &
     &       id000, id010, id100, id110, id200, id210, jpl, jplp,       &
     &       id001, id011, id101, id111, id201, id211, ig, js, js1

      real (kind=kind_phys) :: tauself, taufor, refrat_planck_a,        &
     &       speccomb,       specparm,       specmult,       fs,        &
     &       speccomb1,      specparm1,      specmult1,      fs1,       &
     &       speccomb_planck,specparm_planck,specmult_planck,fpl,       &
     &       fac000, fac100, fac200, fac010, fac110, fac210,            &
     &       fac001, fac101, fac201, fac011, fac111, fac211,            &
     &       p0, p40, fk00, fk10, fk20, p1, p41, fk01, fk11, fk21
!
!===> ...  begin here
!
!  --- ...  calculate reference ratio to be used in calculation of Planck
!           fraction in lower atmosphere.

      refrat_planck_a = chi_mls(1,6)/chi_mls(6,6)        ! P = 387. mb (Level 6)

!  --- ...  lower atmosphere loop

      do k = 1, laytrop
        speccomb = colamt(k,1) + rfrate(k,4,1)*colamt(k,5)
        specparm = colamt(k,1) / speccomb
        specmult = 8.0 * min(specparm, oneminus)
        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        ind0 = ((jp(k)-1)*5 + (jt(k)-1)) * nspa(16) + js

        speccomb1 = colamt(k,1) + rfrate(k,4,2)*colamt(k,5)
        specparm1 = colamt(k,1) / speccomb1
        specmult1 = 8.0 * min(specparm1, oneminus)
        js1 = 1 + int(specmult1)
        fs1 = mod(specmult1, f_one)
        ind1 = (jp(k)*5 + (jt1(k)-1)) * nspa(16) + js1

        speccomb_planck = colamt(k,1) + refrat_planck_a*colamt(k,5)
        specparm_planck = colamt(k,1) / speccomb_planck
        specmult_planck = 8.0 * min(specparm_planck, oneminus)
        jpl = 1 + int(specmult_planck)
        fpl = mod(specmult_planck, f_one)

        inds = indself(k)
        indf = indfor(k)
        indsp = inds + 1
        indfp = indf + 1
        jplp  = jpl  + 1

        if (specparm < 0.125 .and. specparm1 < 0.125) then
          p0 = fs - f_one
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = fs1 - f_one
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0 + 2
          id210 = ind0 +11

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1 + 2
          id211 = ind1 +11
        elseif (specparm > 0.875 .and. specparm1 > 0.875) then
          p0 = -fs
          p40 = p0**4
          fk00 = p40
          fk10 = f_one - p0 - 2.0*p40
          fk20 = p0 + p40

          p1 = -fs1
          p41 = p1**4
          fk01 = p41
          fk11 = f_one - p1 - 2.0*p41
          fk21 = p1 + p41

          id000 = ind0 + 1
          id010 = ind0 +10
          id100 = ind0
          id110 = ind0 + 9
          id200 = ind0 - 1
          id210 = ind0 + 8

          id001 = ind1 + 1
          id011 = ind1 +10
          id101 = ind1
          id111 = ind1 + 9
          id201 = ind1 - 1
          id211 = ind1 + 8
        else
          fk00 = f_one - fs
          fk10 = fs
          fk20 = f_zero

          fk01 = f_one - fs1
          fk11 = fs1
          fk21 = f_zero

          id000 = ind0
          id010 = ind0 + 9
          id100 = ind0 + 1
          id110 = ind0 +10
          id200 = ind0
          id210 = ind0

          id001 = ind1
          id011 = ind1 + 9
          id101 = ind1 + 1
          id111 = ind1 +10
          id201 = ind1
          id211 = ind1
        endif

        fac000 = fk00 * fac00(k)
        fac100 = fk10 * fac00(k)
        fac200 = fk20 * fac00(k)
        fac010 = fk00 * fac10(k)
        fac110 = fk10 * fac10(k)
        fac210 = fk20 * fac10(k)

        fac001 = fk01 * fac01(k)
        fac101 = fk11 * fac01(k)
        fac201 = fk21 * fac01(k)
        fac011 = fk01 * fac11(k)
        fac111 = fk11 * fac11(k)
        fac211 = fk21 * fac11(k)

        do ig = 1, ng16
          tauself = selffac(k)* (selfref(ig,inds) + selffrac(k)         &
     &            * (selfref(ig,indsp) - selfref(ig,inds)))
          taufor  = forfac(k) * (forref(ig,indf) + forfrac(k)           &
     &            * (forref(ig,indfp) - forref(ig,indf)))

          taug(ns16+ig,k) = speccomb                                    &
     &                * (fac000*absa(ig,id000) + fac010*absa(ig,id010)  &
     &                +  fac100*absa(ig,id100) + fac110*absa(ig,id110)  &
     &                +  fac200*absa(ig,id200) + fac210*absa(ig,id210)) &
     &                +     speccomb1                                   &
     &                * (fac001*absa(ig,id001) + fac011*absa(ig,id011)  &
     &                +  fac101*absa(ig,id101) + fac111*absa(ig,id111)  &
     &                +  fac201*absa(ig,id201) + fac211*absa(ig,id211)) &
     &                + tauself + taufor

          fracs(ns16+ig,k) = fracrefa(ig,jpl) + fpl                     &
     &                     * (fracrefa(ig,jplp) - fracrefa(ig,jpl))
        enddo
      enddo

!  --- ...  upper atmosphere loop

      do k = laytrop+1, nlay
        ind0 = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(16) + 1
        ind1 = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(16) + 1

        ind0p = ind0 + 1
        ind1p = ind1 + 1

        do ig = 1, ng16
          taug(ns16+ig,k) = colamt(k,5)                                 &
     &           * (fac00(k)*absb(ig,ind0) + fac10(k)*absb(ig,ind0p)    &
     &           +  fac01(k)*absb(ig,ind1) + fac11(k)*absb(ig,ind1p))

          fracs(ns16+ig,k) = fracrefb(ig)
        enddo
      enddo

! ..................................
      end subroutine taugb16
! ----------------------------------

! ..................................
      end subroutine taumol
!! @}
!-----------------------------------


!
!........................................!
      end module rrtmg_lw                !
!========================================!

