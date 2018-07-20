!>  \file radiation_aerosols.f
!!  This file contains climatological atmospheric aerosol schemes for
!!  radiation computations.

!  ==========================================================  !!!!!
!            'module_radiation_aerosols' description           !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!   this module contains climatological atmospheric aerosol schemes for!
!   radiation computations.                                            !
!                                                                      !
!   in the module, the externally callable subroutines are :           !
!                                                                      !
!      'aer_init'   -- initialization                                  !
!         inputs:                                                      !
!           ( NLAY, me )                                               !
!         outputs:                                                     !
!           ( none )                                                   !
!                                                                      !
!      'aer_update' -- updating aerosol data                           !
!         inputs:                                                      !
!           ( iyear, imon, me )                                        !
!         outputs:                                                     !
!           ( none )                                                   !
!                                                                      !
!      'setaer'     -- mapping aeros profile, compute aeros opticals   !
!         inputs:                                                      !
!           (prsi,prsl,prslk,tvly,rhlay,slmsk,tracer,xlon,xlat,        !
!            IMAX,NLAY,NLP1, lsswr,lslwr,                              !
!         outputs:                                                     !
!           (aerosw,aerolw,tau_gocart)                                 !
!!          (aerosw,aerolw,aerodp)                                     !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!       'module physparam'               in 'physparam.f'              !
!       'module physcons'                in 'physcons.f'               !
!       'module module_radsw_parameters' in 'radsw_xxxx#_param.f'      !
!       'module module_radlw_parameters' in 'radlw_xxxx#_param.f'      !
!       'module module_radlw_cntr_para'  in 'radsw_xxxx#_param.f'      !
!                                                                      !
!   output variable definitions:                                       !
!       aerosw(IMAX,NLAY,NBDSW,1) - aerosols optical depth for sw      !
!       aerosw(IMAX,NLAY,NBDSW,2) - aerosols single scat albedo for sw !
!       aerosw(IMAX,NLAY,NBDSW,3) - aerosols asymmetry parameter for sw!
!                                                                      !
!       aerolw(IMAX,NLAY,NBDLW,1) - aerosols optical depth for lw      !
!       aerolw(IMAX,NLAY,NBDLW,2) - aerosols single scattering albedo  !
!       aerolw(IMAX,NLAY,NBDLW,3) - aerosols asymmetry parameter       !
!                                                                      !
!                                                                      !
!   program history:                                                   !
!     apr     2003  ---  y.-t. hou     created                         !
!     nov 04, 2003  ---  y.-t. hou     modified version                !
!     apr 15, 2005  ---  y.-t. hou     modified module structure       !
!     jul     2006  ---  y.-t. hou     add volcanic forcing            !
!     feb     2007  ---  y.-t. hou     add generalized spectral band   !
!                   interpolation for sw aerosol optical properties    !
!     mar     2007  ---  y.-t. hou     add generalized spectral band   !
!                   interpolation for lw aerosol optical properties    !
!     aug     2007  ---  y.-t. hou     change clim-aer vert domain     !
!                   from pressure reference to sigma reference         !
!     sep     2007  ---  y.-t. hou     moving temporary allocatable    !
!                   module variable arrays to subroutine dynamically   !
!                   allocated arrays (eliminate deallocate calls)      !
!     jan     2010  ---  sarah lu      add gocart option               !
!     may     2010  ---  sarah lu      add geos4-gocart climo          !
!     jul     2010  --   s. moorthi - merged NEMS version with new GFS !
!                        version                                       !
!     oct 23, 2010  ---  Hsin-mu lin   modified subr setclimaer to     !
!        interpolate the 5 degree aerosol data to small domain based on!
!        the nearby 4 points instead of previous nearby assignment by  !
!        using the 5 degree data. this process will eliminate the dsw  !
!        jagged edges in the east conus where aerosol effect are lagre.!
!     dec     2010  ---  y.-t. hou     modified and optimized bi-linear!
!        horizontal interpolation in subr setclimaer. added safe guard !
!        measures in lat/lon indexing and added sea/land mask variable !
!        slmsk as input field to help aerosol profile selection.       !
!     jan     2011  ---  y.-t. hou     divided the program into two    !
!        separated interchangeable modules: a climatology aerosol      !
!        module, and a gocart aerosol scheme module. the stratospheric !
!        volcanic aerosol part is still within the two driver modules, !
!        and may also become a separate one in the further development.!
!        unified in/out argument list for both clim and gocart types of!
!        schemes and added vertically integrated aer-opt-dep, aerodp,  !
!        to replace tau_gocart as optional output for various species. !
!     aug     2012  ---  y.-t. hou     changed the initialization subr !
!        'aerinit' into two parts: 'aer_init' is called at the start   !
!        of run to set up module parameters; and 'aer_update' is       !
!        called within the time loop to check and update data sets.    !
!     nov     2012  ---  y.-t. hou     modified control parameters thru!
!        module 'physparam'.                                            !
!     jan     2013  ---  sarah lu and y.-t. hou   reintegrate both     !
!        opac-clim and gocart schemes into one module to make the      !
!        program best utilize common components. added aerosol model   !
!        scheme selection control variable iaer_mdl to the namelist.   !
!      Aug     2013  --- s. moorthi - merge sarah's gocart changes with!
!                                     yutai's changes                  !
!      13Feb2014  --- Sarah lu - compute aod at 550nm                  !
!                                                                      !
!   references for opac climatological aerosols:                       !
!     hou et al. 2002  (ncep office note 441)                          !
!     hess et al. 1998 - bams v79 831-844                              !
!                                                                      !
!   references for gocart interactive aerosols:                        !
!     chin et al., 2000 - jgr, v105, 24671-24687                       !
!                                                                      !
!   references for stratosperic volcanical aerosols:                   !
!     sato et al. 1993 - jgr, v98, d12, 22987-22994                    !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!> \ingroup RRTMG
!! \defgroup module_radiation_aerosols RRTMG Aerosols Module
!! \brief This module contains climatological atmospheric aerosol schemes for
!!  radiation computations.
!!
!!\version NCEP-Radiation_aerosols  v5.2  Jan 2013
!!
!!\n This module has three externally callable subroutines:
!! - aer_init() -- initialization; called at the start of run to set up
!!                 module parameters
!! - aer_update() -- updating aerosol data; called within the time loop
!!                   to check and update data sets
!! - setaer() -- mapping aeros profile, compute aeros opticals
!!
!!\n References:
!! - OPAC climatological aerosols: \cite hou_et_al_2002; \cite hess_et_al_1998
!! - GOCART interactive aerosols: \cite chin_et_al_2000
!! - Stratospheric volcanical aerosols: \cite sato_et_al_1993
!========================================!
      module module_radiation_aerosols   !
!........................................!
!
      use physparam,only : iaermdl, iaerflg, lavoflg, lalwflg, laswflg, &
     &                     lalw1bd, aeros_file, ivflip, kind_phys       &
     &,                    kind_io4, kind_io8
      use physcons, only : con_pi, con_rd, con_g, con_t0c, con_c,       &
     &                     con_boltz, con_plnk, con_amd

      use module_iounitdef,        only : NIAERCM
      use module_radsw_parameters, only : NBDSW,  wvnsw1=>wvnum1,       &
     &                                    NSWSTR, wvnsw2=>wvnum2
      use module_radlw_parameters, only : NBDLW,  wvnlw1, wvnlw2
!
      use funcphys,                     only : fpkap
      use gfs_phy_tracer_config,        only : gfs_phy_tracer, trcindx
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGAER='NCEP-Radiation_aerosols  v5.2  Jan 2013 '
!    &   VTAGAER='NCEP-Radiation_aerosols  v5.1  Nov 2012 '
!    &   VTAGAER='NCEP-Radiation_aerosols  v5.0  Aug 2012 '

! --- general use parameter constants:
! num of output fields for SW rad
      integer, parameter, public :: NF_AESW = 3
! num of output fields for LW rad
      integer, parameter, public :: NF_AELW = 3
! starting band number in ir region
      integer, parameter, public :: NLWSTR  = 1
! num of species for output aod (opnl)
      integer, parameter, public :: NSPC    = 5
! total+species
      integer, parameter, public :: NSPC1   = NSPC + 1

      real (kind=kind_phys), parameter :: f_zero = 0.0
      real (kind=kind_phys), parameter :: f_one  = 1.0

!  ---  module control parameters set in subroutine "aer_init"
! number of actual bands for sw aerosols; calculated according to
!! laswflg setting
      integer, save :: NSWBND  = NBDSW
! number of actual bands for lw aerosols; calculated according to
!! lalwflg and lalw1bd settings
      integer, save :: NLWBND  = NBDLW
! total number of bands for sw+lw aerosols
      integer, save :: NSWLWBD = NBDSW+NBDLW

! --------------------------------------------------------------------- !
!   section-1 : module variables for spectral band interpolation        !
!               similar to gfdl-sw treatment (2000 version)             !
! --------------------------------------------------------------------- !

!  ---  parameter constants:
! num of wvnum regions where solar flux is constant
      integer, parameter, public :: NWVSOL  = 151

! total num of wvnum included
      integer, parameter, public :: NWVTOT  = 57600
! total num of wvnum in ir range
      integer, parameter, public :: NWVTIR  = 4000

! number of wavenumbers in each region where the solar flux is constant
      integer, dimension(NWVSOL), save :: nwvns0

      data nwvns0   / 100,  11,  14,  18,  24,  33,  50,  83,  12,  12, &
     &  13,  15,  15,  17,  18,  20,  21,  24,  26,  30,  32,  37,  42, &
     &  47,  55,  64,  76,  91, 111, 139, 179, 238, 333,  41,  42,  45, &
     &  46,  48,  51,  53,  55,  58,  61,  64,  68,  71,  75,  79,  84, &
     &  89,  95, 101, 107, 115, 123, 133, 142, 154, 167, 181, 197, 217, &
     & 238, 263, 293, 326, 368, 417, 476, 549, 641, 758, 909, 101, 103, &
     & 105, 108, 109, 112, 115, 117, 119, 122, 125, 128, 130, 134, 137, &
     & 140, 143, 147, 151, 154, 158, 163, 166, 171, 175, 181, 185, 190, &
     & 196, 201, 207, 213, 219, 227, 233, 240, 248, 256, 264, 274, 282, &
     & 292, 303, 313, 325, 337, 349, 363, 377, 392, 408, 425, 444, 462, &
     & 483, 505, 529, 554, 580, 610, 641, 675, 711, 751, 793, 841, 891, &
     & 947,1008,1075,1150,1231,1323,1425,1538,1667,1633,14300 /

! solar flux \f$w/m^2\f$ in each wvnumb region where it is constant
      real (kind=kind_phys), dimension(NWVSOL), save :: s0intv

      data  s0intv(  1: 50)       /                                     &
     &     1.60000E-6, 2.88000E-5, 3.60000E-5, 4.59200E-5, 6.13200E-5,  &
     &     8.55000E-5, 1.28600E-4, 2.16000E-4, 2.90580E-4, 3.10184E-4,  &
     &     3.34152E-4, 3.58722E-4, 3.88050E-4, 4.20000E-4, 4.57056E-4,  &
     &     4.96892E-4, 5.45160E-4, 6.00600E-4, 6.53600E-4, 7.25040E-4,  &
     &     7.98660E-4, 9.11200E-4, 1.03680E-3, 1.18440E-3, 1.36682E-3,  &
     &     1.57560E-3, 1.87440E-3, 2.25500E-3, 2.74500E-3, 3.39840E-3,  &
     &     4.34000E-3, 5.75400E-3, 7.74000E-3, 9.53050E-3, 9.90192E-3,  &
     &     1.02874E-2, 1.06803E-2, 1.11366E-2, 1.15830E-2, 1.21088E-2,  &
     &     1.26420E-2, 1.32250E-2, 1.38088E-2, 1.44612E-2, 1.51164E-2,  &
     &     1.58878E-2, 1.66500E-2, 1.75140E-2, 1.84450E-2, 1.94106E-2 /
      data  s0intv( 51:100)       /                                     &
     &     2.04864E-2, 2.17248E-2, 2.30640E-2, 2.44470E-2, 2.59840E-2,  &
     &     2.75940E-2, 2.94138E-2, 3.13950E-2, 3.34800E-2, 3.57696E-2,  &
     &     3.84054E-2, 4.13490E-2, 4.46880E-2, 4.82220E-2, 5.22918E-2,  &
     &     5.70078E-2, 6.19888E-2, 6.54720E-2, 6.69060E-2, 6.81226E-2,  &
     &     6.97788E-2, 7.12668E-2, 7.27100E-2, 7.31610E-2, 7.33471E-2,  &
     &     7.34814E-2, 7.34717E-2, 7.35072E-2, 7.34939E-2, 7.35202E-2,  &
     &     7.33249E-2, 7.31713E-2, 7.35462E-2, 7.36920E-2, 7.23677E-2,  &
     &     7.25023E-2, 7.24258E-2, 7.20766E-2, 7.18284E-2, 7.32757E-2,  &
     &     7.31645E-2, 7.33277E-2, 7.36128E-2, 7.33752E-2, 7.28965E-2,  &
     &     7.24924E-2, 7.23307E-2, 7.21050E-2, 7.12620E-2, 7.10903E-2 /
      data  s0intv(101:151)       /                        7.12714E-2,  &
     &     7.08012E-2, 7.03752E-2, 7.00350E-2, 6.98639E-2, 6.90690E-2,  &
     &     6.87621E-2, 6.52080E-2, 6.65184E-2, 6.60038E-2, 6.47615E-2,  &
     &     6.44831E-2, 6.37206E-2, 6.24102E-2, 6.18698E-2, 6.06320E-2,  &
     &     5.83498E-2, 5.67028E-2, 5.51232E-2, 5.48645E-2, 5.12340E-2,  &
     &     4.85581E-2, 4.85010E-2, 4.79220E-2, 4.44058E-2, 4.48718E-2,  &
     &     4.29373E-2, 4.15242E-2, 3.81744E-2, 3.16342E-2, 2.99615E-2,  &
     &     2.92740E-2, 2.67484E-2, 1.76904E-2, 1.40049E-2, 1.46224E-2,  &
     &     1.39993E-2, 1.19574E-2, 1.06386E-2, 1.00980E-2, 8.63808E-3,  &
     &     6.52736E-3, 4.99410E-3, 4.39350E-3, 2.21676E-3, 1.33812E-3,  &
     &     1.12320E-3, 5.59000E-4, 3.60000E-4, 2.98080E-4, 7.46294E-5  /

! --------------------------------------------------------------------- !
!   section-2 : module variables for stratospheric volcanic aerosols    !
!               from historical data (sato et al. 1993)                 !
! --------------------------------------------------------------------- !

!  ---  parameter constants:
! lower limit (year) data available
      integer, parameter :: MINVYR = 1850
! upper limit (year) data available
      integer, parameter :: MAXVYR = 1999

! monthly, 45-deg lat-zone aerosols data set in subroutine 'aer_init'
      integer, allocatable, save :: ivolae(:,:,:)

!  ---  static control variables:
! starting year of data in the input file
      integer :: kyrstr
! ending year of data in the input file
      integer :: kyrend
! the year of data in use in the input file
      integer :: kyrsav
! the month of data in use in the input file
      integer :: kmonsav

! --------------------------------------------------------------------- !
!   section-3 : module variables for opac climatological aerosols       !
!               optical properties (hess et al. 1989)                   !
! --------------------------------------------------------------------- !

!  ---  parameters and constants:
! num of max componets in a profile
      integer, parameter :: NXC = 5     !< num of max componets in a profile
! num of aerosols profile structures
      integer, parameter :: NAE = 7
! num of atmos aerosols domains
      integer, parameter :: NDM = 5
! num of lon-points in glb aeros data set
      integer, parameter :: IMXAE = 72
! num of lat-points in glb aeros data set
      integer, parameter :: JMXAE = 37
! num of bands for clim aer data (opac)
      integer, parameter :: NAERBND=61
! num of rh levels for rh-dep components
      integer, parameter :: NRHLEV =8
! num of rh independent aeros species
      integer, parameter :: NCM1 = 6
! num of rh dependent aeros species
      integer, parameter :: NCM2 = 4
      integer, parameter :: NCM  = NCM1+NCM2

! predefined relative humidity levels
      real (kind=kind_phys), dimension(NRHLEV), save :: rhlev
      data  rhlev (:) / 0.0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99 /

!  ---  the following arrays are for climatological data that are
!           allocated and read in subroutine 'clim_aerinit'.
!   - global aerosol distribution:
!      haer  (NDM,NAE)  - scale height of aerosols (km)
!      prsref(NDM,NAE)  - ref pressure lev (sfc to toa) in mb (100Pa)
!      sigref(NDM,NAE)  - ref sigma lev (sfc to toa)

! scale height of aerosols (km)
      real (kind=kind_phys), save, dimension(NDM,NAE) :: haer
! ref pressure lev (sfc to toa) in mb (100Pa)
      real (kind=kind_phys), save, dimension(NDM,NAE) :: prsref
! ref sigma lev (sfc to toa)
      real (kind=kind_phys), save, dimension(NDM,NAE) :: sigref

!  ---  the following arrays are allocate and setup in subr 'clim_aerinit'
!   - for relative humidity independent aerosol optical properties:
!      species : insoluble        (inso); soot             (soot);
!                mineral nuc mode (minm); mineral acc mode (miam);
!                mineral coa mode (micm); mineral transport(mitr).
!      extrhi(NCM1,NSWLWBD) - extinction coefficient for sw+lw spectral band
!      scarhi(NCM1,NSWLWBD) - scattering coefficient for sw+lw spectral band
!      ssarhi(NCM1,NSWLWBD) - single scattering albedo for sw+lw spectral band
!      asyrhi(NCM1,NSWLWBD) - asymmetry parameter for sw+lw spectral band
!   - for relative humidity dependent aerosol optical properties:
!      species : water soluble    (waso); sea salt acc mode(ssam);
!                sea salt coa mode(sscm); sulfate droplets (suso).
!      rh level: 00%, 50%, 70%, 80%, 90%, 95%, 98%, 99%
!      extrhd(NRHLEV,NCM2,NSWLWBD) - extinction coefficient for sw+lw band
!      scarhd(NRHLEV,NCM2,NSWLWBD) - scattering coefficient for sw+lw band
!      ssarhd(NRHLEV,NCM2,NSWLWBD) - single scattering albedo for sw+lw band
!      asyrhd(NRHLEV,NCM2,NSWLWBD) - asymmetry parameter for sw+lw band
!   - for stratospheric aerosols optical properties:
!      extstra(NSWLWBD)            - extinction coefficient for sw+lw band

      real (kind=kind_phys), allocatable, save, dimension(:,:)   ::     &
     &       extrhi, scarhi, ssarhi, asyrhi
      real (kind=kind_phys), allocatable, save, dimension(:,:,:) ::     &
     &       extrhd, scarhd, ssarhd, asyrhd
      real (kind=kind_phys), allocatable, save, dimension(:)     ::     &
     &       extstra

!  ---  the following arrays are calculated in subr 'clim_aerinit'
!   - for topospheric aerosol profile distibution:
!      kprfg (    IMXAE*JMXAE)   - aeros profile index
!      idxcg (NXC*IMXAE*JMXAE)   - aeros component index
!      cmixg (NXC*IMXAE*JMXAE)   - aeros component mixing ratio
!      denng ( 2 *IMXAE*JMXAE)   - aerosols number density

! \name topospheric aerosol profile distribution

! aeros component mixing ratio
      real (kind=kind_phys), dimension(NXC,IMXAE,JMXAE), save :: cmixg
! aeros number density
      real (kind=kind_phys), dimension( 2 ,IMXAE,JMXAE), save :: denng
! aeros component index
      integer,               dimension(NXC,IMXAE,JMXAE), save :: idxcg
! aeros profile index
      integer,               dimension(    IMXAE,JMXAE), save :: kprfg

! --------------------------------------------------------------------- !
!   section-4 : module variables for gocart aerosol optical properties  !
! --------------------------------------------------------------------- !

! \name module variables for gocart aerosol optical properties

!  ---  parameters and constants:
!   - KCM, KCM1, KCM2 are determined from subroutine 'set_aerspc'
! num of bands for aer data (gocart)
      integer, parameter :: KAERBND=61
! num of rh levels for rh-dep components
      integer, parameter :: KRHLEV =36
!*    integer, parameter :: KCM1 = 8   ! num of rh independent aer !species
!*    integer, parameter :: KCM2 = 5   ! num of rh dependent aer species
!*    integer, parameter :: KCM  = KCM1 + KCM2
! num of rh indep aerosols (set in subr set_aerspc)
      integer, save      :: KCM1 = 0
! num of rh dep aerosols (set in subr set_aerspc)
      integer, save      :: KCM2 = 0
! =KCM1+KCM2 (set in subr set_aerspc)
      integer, save      :: KCM

      real (kind=kind_phys), dimension(KRHLEV) :: rhlev_grt             
      data  rhlev_grt (:)/ .00, .05, .10, .15, .20, .25, .30, .35,      &
     &      .40, .45, .50, .55, .60, .65, .70, .75, .80, .81, .82,      &
     &      .83, .84, .85, .86, .87, .88, .89, .90, .91, .92, .93,      &
     &      .94, .95, .96, .97, .98, .99 /

!  ---  the following arrays are allocate and setup in subr 'gocrt_aerinit'
!  ------  gocart aerosol specification    ------
!  =>  transported aerosol species:
!      DU (5-bins)
!      SS (4 bins for climo mode and 5 bins for fcst mode)
!      SU (dms, so2, so4, msa)
!      OC (phobic, philic) and BC (phobic, philic)
!  =>  species and lumped species for aerosol optical properties
!      DU (5-bins, with 4 sub-groups in the submicron bin )
!      SS (ssam for submicron, sscm for coarse mode)
!      SU (so4)
!      OC (phobic, philic) and BC (phobic, philic)
!  =>  specification used for aerosol optical properties luts
!      DU (8 bins)
!      SS (ssam, sscm)
!      SU (suso)
!      OC (waso) and BC (soot)
!
!   - spectral band structure:
!      iendwv_grt(KAERBND)      - ending wavenumber (cm**-1) for each band
!   - relative humidity independent aerosol optical properties:
!   ===> species : dust (8 bins)
!      rhidext0_grt(KAERBND,KCM1) - extinction coefficient
!      rhidssa0_grt(KAERBND,KCM1) - single scattering albedo
!      rhidasy0_grt(KAERBND,KCM1) - asymmetry parameter
!   - relative humidity dependent aerosol optical properties:
!   ===> species : soot, suso, waso, ssam, sscm
!      rhdpext0_grt(KAERBND,KRHLEV,KCM2) - extinction coefficient
!      rhdpssa0_grt(KAERBND,KRHLEV,KCM2) - single scattering albedo
!      rhdpasy0_grt(KAERBND,KRHLEV,KCM2) - asymmetry parameter

! spectral band structure: ending wavenumber (\f$cm^-1\f$) for each band
      integer,               allocatable, dimension(:) :: iendwv_grt
! relative humidity independent aerosol optical properties:
!! species : dust (8 bins)

! \name relative humidity independent aerosol optical properties:
! species : dust (8 bins)

! extinction coefficient
      real (kind=kind_phys),allocatable, dimension(:,:)  :: rhidext0_grt
! single scattering albedo
      real (kind=kind_phys),allocatable, dimension(:,:)  :: rhidssa0_grt
! asymmetry parameter
      real (kind=kind_phys), allocatable, dimension(:,:) :: rhidasy0_grt
!
! relative humidity dependent aerosol optical properties:
! species : soot, suso, waso, ssam, sscm

! \name relative humidity dependent aerosol optical properties:
! species : soot, suso, waso, ssam, sscm

! extinction coefficient
      real (kind=kind_phys),allocatable,dimension(:,:,:) :: rhdpext0_grt
! single scattering albedo
      real (kind=kind_phys),allocatable,dimension(:,:,:) :: rhdpssa0_grt
! asymmetry parameter
      real (kind=kind_phys),allocatable,dimension(:,:,:) :: rhdpasy0_grt

!   - relative humidity independent aerosol optical properties:
!      extrhi_grt(KCM1,NSWLWBD) - extinction coefficient for sw+lw spectral band
!      ssarhi_grt(KCM1,NSWLWBD) - single scattering albedo for sw+lw spectral band
!      asyrhi_grt(KCM1,NSWLWBD) - asymmetry parameter for sw+lw spectral band
!   - relative humidity dependent aerosol optical properties:
!      extrhd_grt(KRHLEV,KCM2,NSWLWBD) - extinction coefficient for sw+lw band
!      ssarhd_grt(KRHLEV,KCM2,NSWLWBD) - single scattering albedo for sw+lw band
!      asyrhd_grt(KRHLEV,KCM2,NSWLWBD) - asymmetry parameter for sw+lw band

!\name relative humidity independent aerosol optical properties

! extinction coefficient for SW+LW spectral band
      real (kind=kind_phys),allocatable,save,dimension(:,:) ::          &
     &      extrhi_grt
! single scattering albedo for SW+LW spectral band
      real (kind=kind_phys),allocatable,save,dimension(:,:) ::          &
     &      ssarhi_grt
! asymmetry parameter for SW+LW spectral band
      real (kind=kind_phys),allocatable,save,dimension(:,:) ::          &
     &      asyrhi_grt

! \name relative humidity dependent aerosol optical properties

! extinction coefficient for SW+LW spectral band
      real (kind=kind_phys),allocatable,save,dimension(:,:,:) ::        &
     &      extrhd_grt
! single scattering albedo for SW+LW band
      real (kind=kind_phys),allocatable,save,dimension(:,:,:) ::        &
     &      ssarhd_grt
! asymmetry parameter for SW+LW band
      real (kind=kind_phys),allocatable,save,dimension(:,:,:) ::        &
     &      asyrhd_grt

! \name module variables for gocart aerosol clim data set

! --------------------------------------------------------------------- !
!   section-5 : module variables for gocart aerosol climo data set      !
! --------------------------------------------------------------------- !
!     This version only supports geos3-gocart data set (Jan 2010)
!     Modified to support geos4-gocart data set        (May 2010)
!
!  geos3-gocart vs geos4-gocart
!  (1) Use the same module variables
!      IMXG,JMXG,KMXG,NMXG,psclmg,dmclmg,geos_rlon,geos_rlat
!  (2) Similarity between geos3 and geos 4:
!      identical lat/lon grids and aerosol specification;
!      direction of vertical index is bottom-up (sfc to toa)
!  (3) Difference between geos3 and geos4
!      vertical coordinate (sigma for geos3/hybrid_sigma_pressure for geos4)
!      aerosol units (mass concentration for geos3/mixing ratio for geos4)

! num of lon-points in geos dataset
      integer, parameter :: IMXG = 144
! num of lat-points in geos dataset
      integer, parameter :: JMXG = 91
! num of vertical layers in geos dataset
      integer, parameter :: KMXG = 30
!*    integer, parameter :: NMXG = 12
! to be determined by set_aerspc
      integer, save      :: NMXG

      real (kind=kind_phys), parameter :: dltx = 360.0 / float(IMXG)
      real (kind=kind_phys), parameter :: dlty = 180.0 / float(JMXG-1)

!  --- the following arrays are allocated and setup in 'rd_gocart_clim'
!   - geos-gocart climo data (input dataset)
!     psclmg  - pressure in cb                   IMXG*JMXG*KMXG
!     dmclmg  - aerosol dry mass in g/m3         IMXG*JMXG*KMXG*NMXG
!               or aerosol mixing ratio in mol/mol or Kg/Kg

! pressure in cb
      real (kind=kind_phys),allocatable, save:: psclmg(:,:,:)
! aerosol dry mass in g/m3 or aerosol mixing ration in mol/mol or Kg/Kg
      real (kind=kind_phys),allocatable, save:: dmclmg(:,:,:,:)

!   - geos-gocart lat/lon arrays
      real (kind=kind_phys), allocatable, save, dimension(:):: geos_rlon
      real (kind=kind_phys), allocatable, save, dimension(:):: geos_rlat

! control flag for gocart climo data set: xxxx as default; ver3 for geos3;
!! ver4 for geos4; 0000 for unknown data
      character*4, save  :: gocart_climo = 'xxxx'

! molecular wght of gocart aerosol species
      real (kind=kind_io4), allocatable :: molwgt(:)

! ---------------------------------------------------------------------
! !
!   section-6 : module variables for gocart aerosol scheme options
!   !
! ---------------------------------------------------------------------
! !

! logical parameter for gocart initialization control
      logical, save :: lgrtint = .true.

! logical parameter for gocart debug print control
!     logical, save :: lckprnt = .true.
      logical, save :: lckprnt = .false.

!  --- the following index/flag/weight are set up in 'set_aerspc'

! merging coefficients for fcst/clim; determined from fdaer
      real (kind=kind_phys), save :: ctaer = f_zero   ! user specified wgt

! option to get fcst gocart aerosol field
      logical, save :: get_fcst = .true.
! option to get clim gocart aerosol field
      logical, save :: get_clim = .true.

!  ------  gocart aerosol specification    ------
!  =>  transported aerosol species:
!      DU (5-bins)
!      SS (4 bins for climo mode and 5 bins for fcst mode)
!      SU (dms, so2, so4, msa)
!      OC (phobic, philic) and BC (phobic, philic)
!  =>  species and lumped species for aerosol optical properties
!      DU (5-bins, with 4 sub-groups in the submicron bin )
!      SS (ssam for submicron, sscm for coarse mode)
!      SU (so4)
!      OC (phobic, philic) and BC (phobic, philic)
!  =>  specification used for aerosol optical properties luts
!      DU (8 bins)
!      SS (ssam, sscm)
!      SU (suso)
!      OC (waso) and BC (soot)
!

! index for rh dependent aerosol optical properties (2nd
! dimension for extrhd_grt, ssarhd_grt, and asyrhd_grt)
      integer, save :: isoot, iwaso, isuso, issam, isscm

!   - index for rh independent aerosol optical properties (1st
!     dimension for extrhi_grt, ssarhi_grt, and asyrhi_grt) is
!     not needed ===> hardwired to 8-bin dust

!   - index for gocart aerosol species to be included in the
!     calculations of aerosol optical properties (ext, ssa, asy)
      type  gocart_index_type
! dust
         integer :: dust1, dust2, dust3, dust4, dust5
! sea salt
         integer :: ssam,  sscm
! sulfate
         integer :: suso
! oc
         integer :: waso_phobic, waso_philic
! bc
         integer :: soot_phobic, soot_philic
      endtype
      type (gocart_index_type), save :: dm_indx

!  index for gocart aerosols from prognostic tracer fields
      type  tracer_index_type
! dust
         integer :: du001, du002, du003, du004, du005
! sea salt
         integer :: ss001, ss002, ss003, ss004, ss005
! sulfate
         integer :: so4
! oc
         integer :: ocphobic, ocphilic
! bc
         integer :: bcphobic, bcphilic
      endtype
      type (tracer_index_type), save :: dmfcs_indx

!   - grid components to be included in the aeropt calculations
! number of aerosol grid components
      integer, save                  :: num_gridcomp = 0
! aerosol grid components
      character, allocatable , save  :: gridcomp(:)*2

! default full-package setting
      integer, parameter          :: max_num_gridcomp = 5
! data max_gridcomp  /'DU', 'BC', 'OC', 'SU', 'SS'/
      character*2                 :: max_gridcomp(max_num_gridcomp)
      data max_gridcomp  /'DU', 'BC', 'OC', 'SU', 'SS'/

! GOCART code modification end here (Sarah Lu)
! ------------------------!
! =======================================================================

!! ---  the following are for diagnostic purpose to output aerosol optical depth
!       aod from 10 components are grouped into 5 major different species:
!      1:dust (inso,minm,miam,micm,mitr); 2:black carbon (soot)
!      3:water soluble (waso);            4:sulfate (suso);      5:sea salt (ssam,sscm)
!
!      idxspc (NCM)         - index conversion array
!      lspcaod              - logical flag for aod from individual species
!
! index conversion array:data  idxspc / 1, 2, 1, 1, 1, 1, 3, 5, 5, 4 /
      integer, dimension(NCM) :: idxspc
      data  idxspc / 1, 2, 1, 1, 1, 1, 3, 5, 5, 4 /
!
!   - wvn550 is the wavenumber (1/cm) of wavelenth 550nm for diagnostic aod output
!     nv_aod is the sw spectral band covering wvn550 (comp in aer_init)
!
! the wavenumber (\f$cm^-1\f$) of wavelength 550nm for diagnostic aod output
      real (kind=kind_phys), parameter :: wvn550 = 1.0e4/0.55
! the sw spectral band covering wvn550 (comp in aer_init)
      integer, save      :: nv_aod = 1

!  ---  public interfaces

      public aer_init, aer_update, setaer


! =================
      contains
! =================

!>\ingroup module_radiation_aerosols
!> The initialization program is to set up necessary parameters and
!! working arrays.
!!
!>\param NLAY    number of model vertical layers (not used)
!>\param me      print message control flag
!>\section aer_init_gen_al aer_init General Algorithm
!! @{
!-----------------------------------
      subroutine aer_init                                               &
     &     ( NLAY, me ) !  ---  inputs
!  ---  outputs: ( to module variables )

!  ==================================================================  !
!                                                                      !
!  aer_init is the initialization program to set up necessary          !
!    parameters and working arrays.                                    !
!                                                                      !
!  inputs:                                                             !
!     NLAY    - number of model vertical layers  (not used)            !
!     me      - print message control flag                             !
!                                                                      !
!  outputs: (to module variables)                                      !
!                                                                      !
!  external module variables: (in physparam)                           !
!     iaermdl - tropospheric aerosol model scheme flag                 !
!               =0 opac-clim; =1 gocart-clim, =2 gocart-prognostic     !
!     lalwflg - logical lw aerosols effect control flag                !
!               =t compute lw aerosol optical prop                     !
!     laswflg - logical sw aerosols effect control flag                !
!               =t compute sw aerosol optical prop                     !
!     lavoflg - logical stratosphere volcanic aerosol control flag     !
!               =t include volcanic aerosol effect                     !
!     lalw1bd = logical lw aeros propty 1 band vs multi-band cntl flag !
!               =t use 1 broad band optical property                   !
!               =f use multi bands optical property                    !
!                                                                      !
!  module constants:                                                   !
!     NWVSOL  - num of wvnum regions where solar flux is constant      !
!     NWVTOT  - total num of wave numbers used in sw spectrum          !
!     NWVTIR  - total num of wave numbers used in the ir region        !
!     NSWBND  - total number of sw spectral bands                      !
!     NLWBND  - total number of lw spectral bands                      !
!                                                                      !
!  usage:    call aer_init                                             !
!                                                                      !
!  subprograms called:  clim_aerinit, gcrt_aerinit,                    !
!                       wrt_aerlog, set_volcaer, set_spectrum,         !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      integer,  intent(in) :: NLAY, me

!  ---  output: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(NWVTOT) :: solfwv        ! one wvn sol flux
      real (kind=kind_phys), dimension(NWVTIR) :: eirfwv        ! one wvn ir flux
!
!===>  ...  begin here
!
      kyrstr  = 1
      kyrend  = 1
      kyrsav  = 1
      kmonsav = 1

!> -# Call wrt_aerlog() to write aerosol parameter configuration to output logs.

      if ( me == 0 ) then

        call wrt_aerlog      ! write aerosol param info to log file
!  ---  inputs:   (in scope variables)
!  ---  outputs:  ( none )

      endif

      if ( iaerflg == 0 ) return      ! return without any aerosol calculations

!  --- ...  in sw, aerosols optical properties are computed for each radiation
!           spectral band; while in lw, optical properties can be calculated
!           for either only one broad band or for each of the lw radiation bands

      if ( laswflg ) then
        NSWBND = NBDSW
      else
        NSWBND = 0
      endif

      if ( lalwflg ) then
        if ( lalw1bd ) then
          NLWBND = 1
        else
          NLWBND = NBDLW
        endif
      else
        NLWBND = 0
      endif

      NSWLWBD = NSWBND + NLWBND

      if ( iaerflg /= 100 ) then

!> -# Call set_spectrum() to set up spectral one wavenumber solar/IR
!! fluxes.

        call set_spectrum
!  ---  inputs:   (module constants)
!  ---  outputs:  (in-scope variables)

!> -# Call clim_aerinit() to invoke tropospheric aerosol initialization.

        if ( iaermdl == 0 ) then                    ! opac-climatology scheme

          call clim_aerinit                                             &
!  ---  inputs:
     &     ( solfwv, eirfwv, me                                         &
!  ---  outputs:
     &     )

!       elseif ( iaermdl == 1 ) then                ! gocart-climatology scheme
!       elseif ( iaermdl==1 .or. iaermdl==2 ) then  ! gocart-clim/prog scheme

!         call gcrt_climinit

!       elseif ( iaermdl == 2 ) then                ! gocart-prognostic scheme

!         call gcrt_aerinit

        else
          if ( me == 0 ) then
            print *,'  !!! ERROR in aerosol model scheme selection',    &
     &              ' iaermdl =',iaermdl
            stop
          endif
        endif

      endif    ! end if_iaerflg_block

!> -# Call set_volcaer() to invoke stratospheric volcanic aerosol
!! initialization.

      if ( lavoflg ) then

        call set_volcaer
!  ---  inputs:  (module variables)
!  ---  outputs: (module variables)

      endif    ! end if_lavoflg_block


! =================
      contains
! =================

!>\ingroup module_radiation_aerosols
!> This subroutine writes aerosol parameter configuration to run log file.
!>\section wrt_aerlog_gen wrt_aerlog General Algorithm
!! @{
!--------------------------------
      subroutine wrt_aerlog
!................................
!  ---  inputs:    (in scope variables)
!  ---  outputs:   ( none )

!  ==================================================================  !
!                                                                      !
!  subprogram : wrt_aerlog                                             !
!                                                                      !
!    write aerosol parameter configuration to run log file.            !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  external module variables:  (in physparam)                          !
!   iaermdl  - aerosol scheme flag: 0:opac-clm; 1:gocart-clim;         !
!                                   2:gocart-prog                      !
!   iaerflg  - aerosol effect control flag: 3-digits (volc,lw,sw)      !
!   lalwflg  - toposphere lw aerosol effect: =f:no; =t:yes             !
!   laswflg  - toposphere sw aerosol effect: =f:no; =t:yes             !
!   lavoflg  - stratospherer volcanic aeros effect: =f:no; =t:yes      !
!                                                                      !
!  outputs: ( none )                                                   !
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call wrt_aerlog                                           !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: ( none )
!  ---  output: ( none )
!  ---  locals:

!
!===>  ...  begin here
!
      print *, VTAGAER    ! print out version tag

      if ( iaermdl == 0 ) then
        print *,' - Using OPAC-seasonal climatology for tropospheric',  &
     &          ' aerosol effect'
      elseif ( iaermdl == 1 ) then
        print *,' - Using GOCART-climatology for tropospheric',         &
     &          ' aerosol effect'
      elseif ( iaermdl == 2 ) then
        print *,' - Using GOCART-prognostic aerosols for tropospheric', &
     &          ' aerosol effect'
      else
        print *,' !!! ERROR in selection of aerosol model scheme',      &
     &          ' IAER_MDL =',iaermdl
        stop
      endif   ! end_if_iaermdl_block

      print *,'   IAER=',iaerflg,'  LW-trop-aer=',lalwflg,              &
     &        '  SW-trop-aer=',laswflg,'  Volc-aer=',lavoflg

      if ( iaerflg <= 0 ) then        ! turn off all aerosol effects
        print *,' - No tropospheric/volcanic aerosol effect included'
        print *,'      Input values of aerosol optical properties to'   &
     &         ,' both SW and LW radiations are set to zeros'
      else
        if ( iaerflg >= 100 ) then    ! incl stratospheric volcanic aerosols
          print *,' - Include stratospheric volcanic aerosol effect'
        else                       ! no stratospheric volcanic aerosols
          print *,' - No stratospheric volcanic aerosol effect'
        endif

        if ( laswflg ) then          ! chcek for sw effect
          print *,'   - Compute multi-band aerosol optical'             &
     &           ,' properties for SW input parameters'
        else
          print *,'   - No SW radiation aerosol effect, values of'      &
     &           ,' aerosol properties to SW input are set to zeros'
        endif

        if ( lalwflg ) then          ! check for lw effect
          if ( lalw1bd ) then
            print *,'   - Compute 1 broad-band aerosol optical'         &
     &           ,' properties for LW input parameters'
          else
            print *,'   - Compute multi-band aerosol optical'           &
     &           ,' properties for LW input parameters'
          endif
        else
          print *,'   - No LW radiation aerosol effect, values of'      &
     &           ,' aerosol properties to LW input are set to zeros'
        endif
      endif     ! end if_iaerflg_block
!
      return
!................................
      end subroutine wrt_aerlog
!! @}
!--------------------------------

!>\ingroup module_radiation_aerosols
!> This subroutine defines the one wavenumber solar fluxes based on toa
!! solar spectral distribution, and define the one wavenumber IR fluxes
!! based on black-body emission distribution at a predefined temperature.
!>\section gel_set_spec set_spectrum General Algorithm
!! @{
!--------------------------------
      subroutine set_spectrum
!................................
!  ---  inputs:   (module constants)
!  ---  outputs:  (in-scope variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : set_spectrum                                           !
!                                                                      !
!    define the one wavenumber solar fluxes based on toa solar spectral!
!    distrobution, and define the one wavenumber ir fluxes based on    !
!    black-body emission distribution at a predefined temperature.     !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!> -  inputs:  (module constants)
!!  -   NWVTOT:  total num of wave numbers used in sw spectrum
!!  -   NWVTIR:  total num of wave numbers used in the ir region
!!
!> -  outputs: (in-scope variables)
!!  -   solfwv(NWVTOT):   solar flux for each individual wavenumber
!!                        (\f$W/m^2\f$)
!!  -   eirfwv(NWVTIR):   ir flux(273k) for each individual wavenumber
!!                        (\f$W/m^2\f$)
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call set_spectrum                                         !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: (module constants)
!     integer :: NWVTOT, NWVTIR

!  ---  output: (in-scope variables)
!     real (kind=kind_phys), dimension(NWVTOT) :: solfwv        ! one wvn sol flux
!     real (kind=kind_phys), dimension(NWVTIR) :: eirfwv        ! one wvn ir flux

!  ---  locals:
      real (kind=kind_phys) :: soltot, tmp1, tmp2, tmp3

      integer :: nb, nw, nw1, nw2, nmax, nmin
!
!===>  ...  begin here
!
!     nmax = min( NWVTOT, nint( maxval(wvnsw2) ))
!     nmin = max( 1,      nint( minval(wvnsw1) ))

!  ---  check print
!     print *,' MINWVN, MAXWVN = ',nmin, nmax
!  --- ...  define the one wavenumber solar fluxes based on toa solar
!           spectral distribution

!     soltot1 = f_zero
!     soltot  = f_zero
      do nb = 1, NWVSOL
        if ( nb == 1 ) then
          nw1 = 1
        else
          nw1 = nw1 + nwvns0(nb-1)
        endif

        nw2 = nw1 + nwvns0(nb) - 1

        do nw = nw1, nw2
          solfwv(nw) = s0intv(nb)
!         soltot1 = soltot1 + s0intv(nb)
!         if ( nw >= nmin .and. nw <= nmax ) then
!           soltot = soltot + s0intv(nb)
!         endif
        enddo
      enddo

!  --- ...  define the one wavenumber ir fluxes based on black-body
!           emission distribution at a predefined temperature

      tmp1 = (con_pi + con_pi) * con_plnk * con_c* con_c
      tmp2 = con_plnk * con_c / (con_boltz * con_t0c)

!$omp parallel do private(nw,tmp3)
      do nw = 1, NWVTIR
        tmp3 = 100.0 * nw
        eirfwv(nw) = (tmp1 * tmp3**3) / (exp(tmp2*tmp3) - 1.0)
      enddo
!
      return
!................................
      end subroutine set_spectrum
!--------------------------------
!! @}

!>\ingroup module_radiation_aerosols
!> The initialization program for stratospheric volcanic aerosols.
!>\section set_volcaer_gen set_volcaer General Algorithm
!! @{
!-----------------------------
      subroutine set_volcaer
!.............................
!  ---  inputs:   ( none )
!  ---  outputs:  (module variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : set_volcaer                                            !
!                                                                      !
!    this is the initialization progrmam for stratospheric volcanic    !
!    aerosols.                                                         !
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call set_volcaer                                          !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: (none)

!  ---  output: (module variables)
!     integer :: ivolae(:,:,:)

!  ---  locals:
!
!===>  ...  begin here
!
!  ---  allocate data space

      if ( .not. allocated(ivolae) ) then
        allocate ( ivolae(12,4,10) )   ! for 12-mon,4-lat_zone,10-year
      endif
!
      return
!................................
      end subroutine set_volcaer
!! @}
!--------------------------------
!
!...................................
      end subroutine aer_init
!-----------------------------------
!!@}

!>\ingroup module_radiation_aerosols
!> This subroutine is the opac-climatology aerosol initialization
!! program to set up necessary parameters and working arrays.
!>\param solfwv     (NWVTOT), solar flux for each individual wavenumber
!!                  \f$(w/m^2)\f$
!!\param eirfwv     (NWVTIR), IR flux(273k) for each individual wavenumber
!!                  \f$(w/m^2)\f$
!!\param me         print message control flag
!!
!!\section gen_clim_aerinit clim_aerinit General Algorithm
!!@{
!-----------------------------------
      subroutine clim_aerinit                                           &
     &     ( solfwv, eirfwv, me                                         &          ! ---  inputs
     &     )                                                           !  ---  outputs

!  ==================================================================  !
!                                                                      !
!  clim_aerinit is the opac-climatology aerosol initialization program !
!  to set up necessary parameters and working arrays.                  !
!                                                                      !
!  inputs:                                                             !
!   solfwv(NWVTOT)   - solar flux for each individual wavenumber (w/m2)!
!   eirfwv(NWVTIR)   - ir flux(273k) for each individual wavenum (w/m2)!
!   me               - print message control flag                      !
!                                                                      !
!  outputs: (to module variables)                                      !
!                                                                      !
!  external module variables: (in physparam)                           !
!     iaerflg - abc 3-digit integer aerosol flag (abc:volc,lw,sw)      !
!               a: =0 use background stratospheric aerosol             !
!                  =1 incl stratospheric vocanic aeros (MINVYR-MAXVYR) !
!               b: =0 no topospheric aerosol in lw radiation           !
!                  =1 include tropspheric aerosols for lw radiation    !
!               c: =0 no topospheric aerosol in sw radiation           !
!                  =1 include tropspheric aerosols for sw radiation    !
!     lalwflg - logical lw aerosols effect control flag                !
!               =t compute lw aerosol optical prop                     !
!     laswflg - logical sw aerosols effect control flag                !
!               =t compute sw aerosol optical prop                     !
!     lalw1bd = logical lw aeros propty 1 band vs multi-band cntl flag !
!               =t use 1 broad band optical property                   !
!               =f use multi bands optical property                    !
!                                                                      !
!  module constants:                                                   !
!     NWVSOL  - num of wvnum regions where solar flux is constant      !
!     NWVTOT  - total num of wave numbers used in sw spectrum          !
!     NWVTIR  - total num of wave numbers used in the ir region        !
!     NSWBND  - total number of sw spectral bands                      !
!     NLWBND  - total number of lw spectral bands                      !
!     NAERBND - number of bands for climatology aerosol data           !
!     NCM1    - number of rh independent aeros species                 !
!     NCM2    - number of rh dependent aeros species                   !
!                                                                      !
!  usage:    call clim_aerinit                                         !
!                                                                      !
!  subprograms called:  set_aercoef, optavg                            !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      real (kind=kind_phys), dimension(:) :: solfwv        ! one wvn sol flux
      real (kind=kind_phys), dimension(:) :: eirfwv        ! one wvn ir flux

      integer,  intent(in) :: me

!  ---  output: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(NAERBND,NCM1)       ::           &
     &       rhidext0, rhidsca0, rhidssa0, rhidasy0
      real (kind=kind_phys), dimension(NAERBND,NRHLEV,NCM2)::           &
     &       rhdpext0, rhdpsca0, rhdpssa0, rhdpasy0
      real (kind=kind_phys), dimension(NAERBND)            :: straext0

      real (kind=kind_phys), dimension(NSWBND,NAERBND) :: solwaer
      real (kind=kind_phys), dimension(NSWBND)         :: solbnd
      real (kind=kind_phys), dimension(NLWBND,NAERBND) :: eirwaer
      real (kind=kind_phys), dimension(NLWBND)         :: eirbnd

      integer, dimension(NSWBND) :: nv1, nv2
      integer, dimension(NLWBND) :: nr1, nr2
!
!===>  ...  begin here
!
!  --- ...  invoke tropospheric aerosol initialization

!> - call set_aercoef() to invoke tropospheric aerosol initialization.
      call set_aercoef
!  ---  inputs:   (in-scope variables, module constants)
!  ---  outputs:  (module variables)


! =================
      contains
! =================

!>\ingroup module_radiation_aerosols
!> The initialization program for climatological aerosols. The program
!! reads and maps the pre-tabulated aerosol optical spectral data onto
!! corresponding SW radiation spectral bands.
!!\section det_set_aercoef set_aercoef General Algorithm
!! @{
!--------------------------------
      subroutine set_aercoef
!................................
!  ---  inputs:   (in-scope variables, module constants)
!  ---  outputs:  (module variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : set_aercoef                                            !
!                                                                      !
!    this is the initialization progrmam for climatological aerosols   !
!                                                                      !
!    the program reads and maps the pre-tabulated aerosol optical      !
!    spectral data onto corresponding sw radiation spectral bands.     !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  inputs:  (in-scope variables, module constants)                     !
!   solfwv(:)    - real, solar flux for individual wavenumber (w/m2)   !
!   eirfwv(:)    - real, lw flux(273k) for individual wavenum (w/m2)   !
!   me           - integer, select cpu number as print control flag    !
!                                                                      !
!  outputs: (to the module variables)                                  !
!                                                                      !
!  external module variables:  (in physparam)                          !
!   lalwflg   - module control flag for lw trop-aer: =f:no; =t:yes     !
!   laswflg   - module control flag for sw trop-aer: =f:no; =t:yes     !
!   aeros_file- external aerosol data file name                        !
!                                                                      !
!  internal module variables:                                          !
!     IMXAE   - number of longitude points in global aeros data set    !
!     JMXAE   - number of latitude points in global aeros data set     !
!     wvnsw1,wvnsw2 (NSWSTR:NSWEND)                                    !
!             - start/end wavenumbers for each of sw bands             !
!     wvnlw1,wvnlw2 (     1:NBDLW)                                     !
!             - start/end wavenumbers for each of lw bands             !
!     NSWLWBD - total num of bands (sw+lw) for aeros optical properties!
!     NSWBND  - number of sw spectral bands actually invloved          !
!     NLWBND  - number of lw spectral bands actually invloved          !
!     NIAERCM - unit number for reading input data set                 !
!     extrhi  - extinction coef for rh-indep aeros         NCM1*NSWLWBD!
!     scarhi  - scattering coef for rh-indep aeros         NCM1*NSWLWBD!
!     ssarhi  - single-scat-alb for rh-indep aeros         NCM1*NSWLWBD!
!     asyrhi  - asymmetry factor for rh-indep aeros        NCM1*NSWLWBD!
!     extrhd  - extinction coef for rh-dep aeros    NRHLEV*NCM2*NSWLWBD!
!     scarhd  - scattering coef for rh-dep aeros    NRHLEV*NCM2*NSWLWBD!
!     ssarhd  - single-scat-alb for rh-dep aeros    NRHLEV*NCM2*NSWLWBD!
!     asyrhd  - asymmetry factor for rh-dep aeros   NRHLEV*NCM2*NSWLWBD!
!                                                                      !
!  major local variables:                                              !
!   for handling spectral band structures                              !
!     iendwv   - ending wvnum (cm**-1) for each band  NAERBND          !
!   for handling optical properties of rh independent species (NCM1)   !
!         1. insoluble        (inso); 2. soot             (soot);      !
!         3. mineral nuc mode (minm); 4. mineral acc mode (miam);      !
!         5. mineral coa mode (micm); 6. mineral transport(mitr).      !
!     rhidext0 - extinction coefficient             NAERBND*NCM1       !
!     rhidsca0 - scattering coefficient             NAERBND*NCM1       !
!     rhidssa0 - single scattering albedo           NAERBND*NCM1       !
!     rhidasy0 - asymmetry parameter                NAERBND*NCM1       !
!   for handling optical properties of rh ndependent species (NCM2)    !
!         1. water soluble    (waso); 2. sea salt acc mode(ssam);      !
!         3. sea salt coa mode(sscm); 4. sulfate droplets (suso).      !
!         rh level (NRHLEV): 00%, 50%, 70%, 80%, 90%, 95%, 98%, 99%    !
!     rhdpext0 - extinction coefficient             NAERBND,NRHLEV,NCM2!
!     rhdpsca0 - scattering coefficient             NAERBND,NRHLEV,NCM2!
!     rhdpssa0 - single scattering albedo           NAERBND,NRHLEV,NCM2!
!     rhdpasy0 - asymmetry parameter                NAERBND,NRHLEV,NCM2!
!   for handling optical properties of stratospheric bkgrnd aerosols   !
!     straext0 - extingction coefficients             NAERBND          !
!                                                                      !
!  usage:    call set_aercoef                                          !
!                                                                      !
!  subprograms called:  optavg                                         !
!                                                                      !
!  ==================================================================  !
!
!  ---  inputs:  ( none )
!  ---  output: ( none )

!  ---  locals:
      integer, dimension(NAERBND) :: iendwv

      integer :: i, j, k, m, mb, ib, ii, id, iw, iw1, iw2

      real (kind=kind_phys) :: sumsol, sumir

      logical :: file_exist
      character :: cline*80
!
!===>  ...  begin here
!
!> -# Reading climatological aerosols optical data from aeros_file,
!! including:

      inquire (file=aeros_file, exist=file_exist)

      if ( file_exist ) then
        close (NIAERCM)
        open  (unit=NIAERCM,file=aeros_file,status='OLD',               &
     &        action='read',form='FORMATTED')
        rewind (NIAERCM)
      else
        print *,'    Requested aerosol data file "',aeros_file,         &
     &          '" not found!'
        print *,'    *** Stopped in subroutine aero_init !!'
        stop
      endif     ! end if_file_exist_block

!  --- ...  skip monthly global distribution

      do m = 1, 12
        read (NIAERCM,12) cline
  12    format(a80/)

        do j = 1, JMXAE
          do i = 1, IMXAE
            read(NIAERCM,*) id
          enddo
        enddo
      enddo   ! end do_m_block

!  --- ...  aloocate and input aerosol optical data

      if ( .not. allocated( extrhi ) ) then
        allocate ( extrhi (       NCM1,NSWLWBD) )
        allocate ( scarhi (       NCM1,NSWLWBD) )
        allocate ( ssarhi (       NCM1,NSWLWBD) )
        allocate ( asyrhi (       NCM1,NSWLWBD) )
        allocate ( extrhd (NRHLEV,NCM2,NSWLWBD) )
        allocate ( scarhd (NRHLEV,NCM2,NSWLWBD) )
        allocate ( ssarhd (NRHLEV,NCM2,NSWLWBD) )
        allocate ( asyrhd (NRHLEV,NCM2,NSWLWBD) )
        allocate ( extstra(            NSWLWBD) )
      endif

!>  - ending wave num for 61 aerosol spectral bands
      read(NIAERCM,21) cline
  21  format(a80)
      read(NIAERCM,22) iendwv(:)
  22  format(13i6)

!>  - atmos scale height for 5 domains, 7 profs
      read(NIAERCM,21) cline
      read(NIAERCM,24) haer(:,:)
  24  format(20f4.1)

!>  - reference pressure for 5 domains, 7 profs
      read(NIAERCM,21) cline
      read(NIAERCM,26) prsref(:,:)
  26  format(10f7.2)

!>  - rh independent ext coef for 61 bands, 6 species
      read(NIAERCM,21) cline
      read(NIAERCM,28) rhidext0(:,:)
  28  format(8e10.3)

!>  - rh independent sca coef for 61 bands, 6 species
      read(NIAERCM,21) cline
      read(NIAERCM,28) rhidsca0(:,:)

!>  - rh independent ssa coef for 61 bands, 6 species
      read(NIAERCM,21) cline
      read(NIAERCM,28) rhidssa0(:,:)

!>  - rh independent asy coef for 61 bands, 6 species
      read(NIAERCM,21) cline
      read(NIAERCM,28) rhidasy0(:,:)

!>  - rh dependent ext coef for 61 bands, 8 rh lev, 4 species
      read(NIAERCM,21) cline
      read(NIAERCM,28) rhdpext0(:,:,:)

!>  - rh dependent sca coef for 61 bands, 8 rh lev, 4 species
      read(NIAERCM,21) cline
      read(NIAERCM,28) rhdpsca0(:,:,:)

!>  - rh dependent ssa coef for 61 bands, 8 rh lev, 4 species
      read(NIAERCM,21) cline
      read(NIAERCM,28) rhdpssa0(:,:,:)

!>  - rh dependent asy coef for 61 bands, 8 rh lev, 4 species
      read(NIAERCM,21) cline
      read(NIAERCM,28) rhdpasy0(:,:,:)

!>  - stratospheric background aeros for 61 bands
      read(NIAERCM,21) cline
      read(NIAERCM,28) straext0(:)

      close (NIAERCM)

!> -# Convert pressure reference level (in mb) to sigma reference level
!!    assume an 1000mb reference surface pressure.

      sigref(:,:) = 0.001 * prsref(:,:)

!> -# Compute solar flux weights and interval indices for mapping
!!    spectral bands between SW radiation and aerosol data.

      if ( laswflg ) then
        solbnd (:)   = f_zero
!$omp parallel do private(i,j)
        do j=1,naerbnd
          do i=1,nswbnd
            solwaer(i,j) = f_zero
          enddo
        enddo

!$omp parallel do private(ib,mb,ii,iw1,iw2,iw,nv_aod,sumsol)
        do ib = 1, NSWBND
          mb = ib + NSWSTR - 1
          ii = 1
          iw1 = nint(wvnsw1(mb))
          iw2 = nint(wvnsw2(mb))

          if ( wvnsw2(mb) >= wvn550 .and. wvn550 >= wvnsw1(mb) ) then
            nv_aod = ib                  ! sw band number covering 550nm wavelenth
          endif

          Lab_swdowhile : do while ( iw1 > iendwv(ii) )
            if ( ii == NAERBND ) exit Lab_swdowhile
            ii = ii + 1
          enddo  Lab_swdowhile

          sumsol = f_zero
          nv1(ib) = ii

          do iw = iw1, iw2
            solbnd(ib) = solbnd(ib) + solfwv(iw)
            sumsol     = sumsol     + solfwv(iw)

            if ( iw == iendwv(ii) ) then
              solwaer(ib,ii) = sumsol

              if ( ii < NAERBND ) then
                sumsol = f_zero
                ii = ii + 1
              endif
            endif
          enddo

          if ( iw2 /= iendwv(ii) ) then
            solwaer(ib,ii) = sumsol
          endif

          nv2(ib) = ii
!         frcbnd(ib) = solbnd(ib) / soltot
        enddo     ! end do_ib_block for sw
      endif    ! end if_laswflg_block

!> -# Compute LW flux weights and interval indices for mapping
!!    spectral bands between lw radiation and aerosol data.

      if ( lalwflg ) then
        eirbnd (:)   = f_zero
!$omp parallel do private(i,j)
        do j=1,naerbnd
          do i=1,nlwbnd
            eirwaer(i,j) = f_zero
          enddo
        enddo

!$omp parallel do private(ib,ii,iw1,iw2,iw,mb,sumir)
        do ib = 1, NLWBND
          ii = 1
          if ( NLWBND == 1 ) then
!           iw1 = 250                   ! corresponding 40 mu
            iw1 = 400                   ! corresponding 25 mu
            iw2 = 2500                  ! corresponding 4  mu
          else
            mb = ib + NLWSTR - 1
            iw1 = nint(wvnlw1(mb))
            iw2 = nint(wvnlw2(mb))
          endif

          Lab_lwdowhile : do while ( iw1 > iendwv(ii) )
            if ( ii == NAERBND ) exit Lab_lwdowhile
            ii = ii + 1
          enddo  Lab_lwdowhile

          sumir = f_zero
          nr1(ib) = ii

          do iw = iw1, iw2
            eirbnd(ib) = eirbnd(ib) + eirfwv(iw)
            sumir  = sumir  + eirfwv(iw)

            if ( iw == iendwv(ii) ) then
              eirwaer(ib,ii) = sumir

              if ( ii < NAERBND ) then
                sumir = f_zero
                ii = ii + 1
              endif
            endif
          enddo

          if ( iw2 /= iendwv(ii) ) then
            eirwaer(ib,ii) = sumir
          endif

          nr2(ib) = ii
        enddo     ! end do_ib_block for lw
      endif    ! end if_lalwflg_block

!> -# Call optavg() to compute spectral band mean properties for each
!! species.

      call optavg
!  ---  inputs:  (in-scope variables, module variables)
!  ---  outputs: (module variables)

!  ---  check print
!     do ib = 1, NSWBND
!       print *,' After optavg, for sw band:',ib
!       print *,'  extrhi:', extrhi(:,ib)
!       print *,'  scarhi:', scarhi(:,ib)
!       print *,'  ssarhi:', ssarhi(:,ib)
!       print *,'  asyrhi:', asyrhi(:,ib)
!       mb = ib + NSWSTR - 1
!       print *,'  wvnsw1,wvnsw2 :',wvnsw1(mb),wvnsw2(mb)
!       do i = 1, NRHLEV
!         print *,'  extrhd for rhlev:',i
!         print *,extrhd(i,:,ib)
!         print *,'  scarhd for rhlev:',i
!         print *,scarhd(i,:,ib)
!         print *,'  ssarhd for rhlev:',i
!         print *,ssarhd(i,:,ib)
!         print *,'  asyrhd for rhlev:',i
!         print *,asyrhd(i,:,ib)
!       enddo
!       print *,' extstra:', extstra(ib)
!     enddo
!     print *,'  wvnlw1 :',wvnlw1
!     print *,'  wvnlw2 :',wvnlw2
!     do ib = 1, NLWBND
!       ii = NSWBND + ib
!       print *,' After optavg, for lw band:',ib
!       print *,'  extrhi:', extrhi(:,ii)
!       print *,'  scarhi:', scarhi(:,ii)
!       print *,'  ssarhi:', ssarhi(:,ii)
!       print *,'  asyrhi:', asyrhi(:,ii)
!       do i = 1, NRHLEV
!         print *,'  extrhd for rhlev:',i
!         print *,extrhd(i,:,ii)
!         print *,'  scarhd for rhlev:',i
!         print *,scarhd(i,:,ii)
!         print *,'  ssarhd for rhlev:',i
!         print *,ssarhd(i,:,ii)
!         print *,'  asyrhd for rhlev:',i
!         print *,asyrhd(i,:,ii)
!       enddo
!       print *,' extstra:', extstra(ii)
!     enddo
!
      return
!................................
      end subroutine set_aercoef
!--------------------------------
!! @}

!>\ingroup module_radiation_aerosols
!> This subroutine computes mean aerosols optical properties over each
!! SW radiation spectral band for each of the species components. This
!! program follows GFDL's approach for thick cloud optical property in
!! SW radiation scheme (2000).
!>\section optave_gen optavg General Algorithm
!! @{
!--------------------------------
      subroutine optavg
!................................
!  ---  inputs:  (in-scope variables, module variables
!  ---  outputs: (module variables)

! ==================================================================== !
!                                                                      !
! subprogram: optavg                                                   !
!                                                                      !
!   compute mean aerosols optical properties over each sw radiation    !
!   spectral band for each of the species components.  This program    !
!   follows gfdl's approach for thick cloud opertical property in      !
!   sw radiation scheme (2000).                                        !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
! major input variables:                                               !
!   nv1,nv2 (NSWBND) - start/end spectral band indices of aerosol data !
!                      for each sw radiation spectral band             !
!   nr1,nr2 (NLWBND) - start/end spectral band indices of aerosol data !
!                      for each ir radiation spectral band             !
!   solwaer (NSWBND,NAERBND)                                           !
!                    - solar flux weight over each sw radiation band   !
!                      vs each aerosol data spectral band              !
!   eirwaer (NLWBND,NAERBND)                                           !
!                    - ir flux weight over each lw radiation band      !
!                      vs each aerosol data spectral band              !
!   solbnd  (NSWBND) - solar flux weight over each sw radiation band   !
!   eirbnd  (NLWBND) - ir flux weight over each lw radiation band      !
!   NSWBND           - total number of sw spectral bands               !
!   NLWBND           - total number of lw spectral bands               !
!                                                                      !
! external module variables:  (in physparam)                           !
!   laswflg          - control flag for sw spectral region             !
!   lalwflg          - control flag for lw spectral region             !
!                                                                      !
! output variables: (to module variables)                              !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
!  ---  output:

!  ---  locals:
      real (kind=kind_phys) :: sumk, sums, sumok, sumokg, sumreft,      &
     &       sp, refb, reft, rsolbd, rirbd

      integer :: ib, nb, ni, nh, nc
!
!===> ...  begin here
!
!  --- ...  loop for each sw radiation spectral band

      if ( laswflg ) then

!$omp parallel do private(nb,nc,sumk,sums,sumok,sumokg,sumreft)
!$omp+private(ni,nh,sp,reft,refb,rsolbd)
        do nb = 1, NSWBND
          rsolbd = f_one / solbnd(nb)

!  ---  for rh independent aerosol species

          do nc = 1, NCM1        !  ---  for rh independent aerosol species
            sumk    = f_zero
            sums    = f_zero
            sumok   = f_zero
            sumokg  = f_zero
            sumreft = f_zero

            do ni = nv1(nb), nv2(nb)
              sp   = sqrt( (f_one - rhidssa0(ni,nc))                    &
     &             / (f_one - rhidssa0(ni,nc)*rhidasy0(ni,nc)) )
              reft = (f_one - sp) / (f_one + sp)
              sumreft = sumreft + reft*solwaer(nb,ni)

              sumk    = sumk    + rhidext0(ni,nc)*solwaer(nb,ni)
              sums    = sums    + rhidsca0(ni,nc)*solwaer(nb,ni)
              sumok   = sumok   + rhidssa0(ni,nc)*solwaer(nb,ni)        &
     &                * rhidext0(ni,nc)
              sumokg  = sumokg  + rhidssa0(ni,nc)*solwaer(nb,ni)        &
     &                * rhidext0(ni,nc)*rhidasy0(ni,nc)
            enddo

            refb = sumreft * rsolbd

            extrhi(nc,nb) = sumk   * rsolbd
            scarhi(nc,nb) = sums   * rsolbd
            asyrhi(nc,nb) = sumokg / (sumok + 1.0e-10)
            ssarhi(nc,nb) = 4.0*refb                                    &
     &         / ( (f_one+refb)**2 - asyrhi(nc,nb)*(f_one-refb)**2 )
          enddo   ! end do_nc_block for rh-ind aeros


          do nc = 1, NCM2        !  ---  for rh dependent aerosols species
            do nh = 1, NRHLEV
              sumk    = f_zero
              sums    = f_zero
              sumok   = f_zero
              sumokg  = f_zero
              sumreft = f_zero

              do ni = nv1(nb), nv2(nb)
                sp   = sqrt( (f_one - rhdpssa0(ni,nh,nc))               &
     &               / (f_one - rhdpssa0(ni,nh,nc)*rhdpasy0(ni,nh,nc)) )
                reft = (f_one - sp) / (f_one + sp)
                sumreft = sumreft + reft*solwaer(nb,ni)

                sumk    = sumk    + rhdpext0(ni,nh,nc)*solwaer(nb,ni)
                sums    = sums    + rhdpsca0(ni,nh,nc)*solwaer(nb,ni)
                sumok   = sumok   + rhdpssa0(ni,nh,nc)*solwaer(nb,ni)   &
     &                  * rhdpext0(ni,nh,nc)
                sumokg  = sumokg  + rhdpssa0(ni,nh,nc)*solwaer(nb,ni)   &
     &                  * rhdpext0(ni,nh,nc)*rhdpasy0(ni,nh,nc)
              enddo

              refb = sumreft * rsolbd

              extrhd(nh,nc,nb) = sumk   * rsolbd
              scarhd(nh,nc,nb) = sums   * rsolbd
              asyrhd(nh,nc,nb) = sumokg / (sumok + 1.0e-10)
              ssarhd(nh,nc,nb) = 4.0*refb                               &
     &         / ( (f_one+refb)**2 - asyrhd(nh,nc,nb)*(f_one-refb)**2 )
            enddo   ! end do_nh_block
          enddo   ! end do_nc_block for rh-dep aeros

!  ---  for stratospheric background aerosols

          sumk = f_zero
          do ni = nv1(nb), nv2(nb)
            sumk = sumk + straext0(ni)*solwaer(nb,ni)
          enddo

          extstra(nb) = sumk * rsolbd

!  ---  check print
!         if ( nb > 6 .and. nb < 10) then
!           print *,' in optavg for sw band',nb
!           print *,'  nv1, nv2:',nv1(nb),nv2(nb)
!           print *,'  solwaer:',solwaer(nb,nv1(nb):nv2(nb))
!           print *,'  extrhi:', extrhi(:,nb)
!           do i = 1, NRHLEV
!             print *,'  extrhd for rhlev:',i
!             print *,extrhd(i,:,nb)
!           enddo
!           print *,'  sumk, rsolbd, extstra:',sumk,rsolbd,extstra(nb)
!         endif

        enddo   !  end do_nb_block for sw
      endif   !  end if_laswflg_block

!  --- ...  loop for each lw radiation spectral band

      if ( lalwflg ) then

!$omp parallel do private(nb,ib,nc,rirbd,sumk,sums,sumok,sumokg,sumreft)
!$omp+private(ni,nh,sp,reft,refb,rsolbd)
        do nb = 1, NLWBND

          ib = NSWBND + nb
          rirbd = f_one / eirbnd(nb)

          do nc = 1, NCM1        !  ---  for rh independent aerosol species
            sumk    = f_zero
            sums    = f_zero
            sumok   = f_zero
            sumokg  = f_zero
            sumreft = f_zero

            do ni = nr1(nb), nr2(nb)
              sp   = sqrt( (f_one - rhidssa0(ni,nc))                    &
     &             / (f_one - rhidssa0(ni,nc)*rhidasy0(ni,nc)) )
              reft = (f_one - sp) / (f_one + sp)
              sumreft = sumreft + reft*eirwaer(nb,ni)

              sumk    = sumk    + rhidext0(ni,nc)*eirwaer(nb,ni)
              sums    = sums    + rhidsca0(ni,nc)*eirwaer(nb,ni)
              sumok   = sumok   + rhidssa0(ni,nc)*eirwaer(nb,ni)        &
     &                * rhidext0(ni,nc)
              sumokg  = sumokg  + rhidssa0(ni,nc)*eirwaer(nb,ni)        &
     &                * rhidext0(ni,nc)*rhidasy0(ni,nc)
            enddo

            refb = sumreft * rirbd

            extrhi(nc,ib) = sumk   * rirbd
            scarhi(nc,ib) = sums   * rirbd
            asyrhi(nc,ib) = sumokg / (sumok + 1.0e-10)
            ssarhi(nc,ib) = 4.0*refb                                       &
     &         / ( (f_one+refb)**2 - asyrhi(nc,ib)*(f_one-refb)**2 )
          enddo   ! end do_nc_block for rh-ind aeros

          do nc = 1, NCM2        !  ---  for rh dependent aerosols species
            do nh = 1, NRHLEV
              sumk    = f_zero
              sums    = f_zero
              sumok   = f_zero
              sumokg  = f_zero
              sumreft = f_zero

              do ni = nr1(nb), nr2(nb)
                sp   = sqrt( (f_one - rhdpssa0(ni,nh,nc))               &
     &             / (f_one - rhdpssa0(ni,nh,nc)*rhdpasy0(ni,nh,nc)) )
                reft = (f_one - sp) / (f_one + sp)
                sumreft = sumreft + reft*eirwaer(nb,ni)

                sumk    = sumk    + rhdpext0(ni,nh,nc)*eirwaer(nb,ni)
                sums    = sums    + rhdpsca0(ni,nh,nc)*eirwaer(nb,ni)
                sumok   = sumok   + rhdpssa0(ni,nh,nc)*eirwaer(nb,ni)   &
     &                  * rhdpext0(ni,nh,nc)
                sumokg  = sumokg  + rhdpssa0(ni,nh,nc)*eirwaer(nb,ni)   &
     &                  * rhdpext0(ni,nh,nc)*rhdpasy0(ni,nh,nc)
              enddo

              refb = sumreft * rirbd

              extrhd(nh,nc,ib) = sumk   * rirbd
              scarhd(nh,nc,ib) = sums   * rirbd
              asyrhd(nh,nc,ib) = sumokg / (sumok + 1.0e-10)
              ssarhd(nh,nc,ib) = 4.0*refb                               &
     &         / ( (f_one+refb)**2 - asyrhd(nh,nc,ib)*(f_one-refb)**2 )
            enddo   ! end do_nh_block
          enddo   ! end do_nc_block for rh-dep aeros

!  ---  for stratospheric background aerosols

          sumk = f_zero
          do ni = nr1(nb), nr2(nb)
            sumk = sumk + straext0(ni)*eirwaer(nb,ni)
          enddo

          extstra(ib) = sumk * rirbd

!  ---  check print
!         if ( nb >= 1 .and. nb < 5) then
!           print *,' in optavg for ir band:',nb
!           print *,'  nr1, nr2:',nr1(nb),nr2(nb)
!           print *,'  eirwaer:',eirwaer(nb,nr1(nb):nr2(nb))
!           print *,'  extrhi:', extrhi(:,ib)
!           do i = 1, NRHLEV
!             print *,'  extrhd for rhlev:',i
!             print *,extrhd(i,:,ib)
!           enddo
!           print *,'  sumk, rirbd, extstra:',sumk,rirbd,extstra(ib)
!         endif

        enddo   !  end do_nb_block for lw
      endif   !  end if_lalwflg_block
!
      return
!................................
      end subroutine optavg
!! @}
!--------------------------------
!
!...................................
      end subroutine clim_aerinit
!-----------------------------------
!!@}

!>\ingroup module_radiation_aerosols
!> This subroutine checks and updates time varying climatology aerosol
!! data sets.
!>\param iyear    4-digit calender year
!!\param imon     month of the year
!!\param me       print message control flag
!>\section gen_aer_upd aer_update General Algorithm
!! @{
!-----------------------------------
      subroutine aer_update                                             &
     &     ( iyear, imon, me ) !  ---  inputs:
!  ---  outputs: ( to module variables )

!  ==================================================================  !
!                                                                      !
!  aer_update checks and update time varying climatology aerosol       !
!    data sets.                                                        !
!                                                                      !
!  inputs:                                          size               !
!     iyear   - 4-digit calender year                 1                !
!     imon    - month of the year                     1                !
!     me      - print message control flag            1                !
!                                                                      !
!  outputs: ( none )                                                   !
!                                                                      !
!  external module variables: (in physparam)                           !
!     lalwflg     - control flag for tropospheric lw aerosol           !
!     laswflg     - control flag for tropospheric sw aerosol           !
!     lavoflg     - control flag for stratospheric volcanic aerosol    !
!                                                                      !
!  usage:    call aero_update                                          !
!                                                                      !
!  subprograms called:  trop_update, volc_update                       !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      integer,  intent(in) :: iyear, imon, me

!  ---  output: ( none )

!  ---  locals: ( none )
!
!===> ...  begin here
!
      if ( imon < 1 .or. imon > 12 ) then
        print *,' ***** ERROR in specifying requested month !!! ',      &
     &          'imon=', imon
        print *,' ***** STOPPED in subroutinte aer_update !!!'
        stop
      endif

!> -# Call trop_update() to update monthly tropospheric aerosol data.
      if ( lalwflg .or. laswflg ) then
        call trop_update
      endif

!> -# Call volc_update() to update yearly stratospheric volcanic aerosol data.
      if ( lavoflg ) then
        call volc_update
      endif


! =================
      contains
! =================

!>\ingroup module_radiation_aerosols
!> This subroutine updates the monthly global distribution of aerosol
!! profiles in five degree horizontal resolution.
!>\section trop_update_gen trop_update General Algorithm
!! @{
!--------------------------------
      subroutine trop_update
!................................
!  ---  inputs:    (in scope variables, module variables)
!  ---  outputs:   (module variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : trop_update                                            !
!                                                                      !
!    updates the  monthly global distribution of aerosol profiles in   !
!    five degree horizontal resolution.                                !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  inputs:  (in-scope variables, module constants)                     !
!   imon     - integer, month of the year                              !
!   me       - integer, print message control flag                     !
!                                                                      !
!  outputs: (module variables)                                         !
!                                                                      !
!  external module variables: (in physparam)                           !
!    aeros_file   - external aerosol data file name                    !
!                                                                      !
!  internal module variables:                                          !
!    kprfg (    IMXAE*JMXAE)   - aeros profile index                   !
!    idxcg (NXC*IMXAE*JMXAE)   - aeros component index                 !
!    cmixg (NXC*IMXAE*JMXAE)   - aeros component mixing ratio          !
!    denng ( 2 *IMXAE*JMXAE)   - aerosols number density               !
!                                                                      !
!    NIAERCM      - unit number for input data set                     !
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call trop_update                                          !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: ( none )
!  ---  output: ( none )

!  ---  locals:
!     real (kind=kind_io8)  :: cmix(NXC), denn, tem
      real (kind=kind_phys) :: cmix(NXC), denn, tem
      integer               :: idxc(NXC), kprf

      integer :: i, id, j, k, m, nc
      logical :: file_exist

      character :: cline*80, ctyp*3
!
!===>  ...  begin here
!
!  --- ...  reading climatological aerosols data

      inquire (file=aeros_file, exist=file_exist)

      if ( file_exist ) then
        close(NIAERCM)
        open (unit=NIAERCM,file=aeros_file,status='OLD',                &
     &        action='read',form='FORMATTED')
        rewind (NIAERCM)

        if ( me == 0 ) then
          print *,'   Opened aerosol data file: ',aeros_file
        endif
      else
        print *,'    Requested aerosol data file "',aeros_file,         &
     &          '" not found!'
        print *,'    *** Stopped in subroutine trop_update !!'
        stop
      endif      ! end if_file_exist_block

!$omp parallel do private(i,j,m)
      do j = 1, JMXAE
        do i = 1, IMXAE
          do m = 1, NXC
            idxcg(m,i,j) = 0
            cmixg(m,i,j) = f_zero
          enddo
        enddo
      enddo

!$omp parallel do private(i,j)
      do j = 1, JMXAE
        do i = 1, IMXAE
          denng(1,i,j) = f_zero
          denng(2,i,j) = f_zero
        enddo
      enddo

!  --- ...  loop over 12 month global distribution

      Lab_do_12mon : do m = 1, 12

        read(NIAERCM,12) cline
  12    format(a80/)

        if ( m /= imon ) then
!         if ( me == 0 ) print *,'  *** Skipped ',cline

          do j = 1, JMXAE
            do i = 1, IMXAE
              read(NIAERCM,*) id
            enddo
          enddo
        else
          if ( me == 0 ) print *,'  --- Reading ',cline

          do j = 1, JMXAE
            do i = 1, IMXAE
              read(NIAERCM,14) (idxc(k),cmix(k),k=1,NXC),kprf,denn,nc,  &
     &                         ctyp
  14          format(5(i2,e11.4),i2,f8.2,i3,1x,a3)

              kprfg(i,j)     = kprf
              denng(1,i,j)   = denn       ! num density of 1st layer
              if ( kprf >= 6 ) then
                denng(2,i,j) = cmix(NXC)  ! num density of 2dn layer
              else
                denng(2,i,j) = f_zero
              endif

              tem = f_one
              do k = 1, NXC-1
                idxcg(k,i,j) = idxc(k)    ! component index
                cmixg(k,i,j) = cmix(k)    ! component mixing ratio
                tem          = tem - cmix(k)
              enddo
              idxcg(NXC,i,j) = idxc(NXC)
              cmixg(NXC,i,j) = tem        ! to make sure all add to 1.
            enddo
          enddo

          close (NIAERCM)
          exit  Lab_do_12mon
        endif     ! end if_m_block

      enddo  Lab_do_12mon

!  --  check print

!     print *,'  IDXCG :'
!     print 16,idxcg
! 16  format(40i3)
!     print *,'  CMIXG :'
!     print 17,cmixg
!     print *,'  DENNG :'
!     print 17,denng
!     print *,'  KPRFG :'
!     print 17,kprfg
! 17  format(8e16.9)
!
      return
!................................
      end subroutine trop_update
!! @}
!--------------------------------

!>\ingroup module_radiation_aerosols
!> This subroutine searches historical volcanic data sets to find and
!! read in monthly 45-degree lat-zone band of optical depth.
!>\section volc_update_gen volc_update General Algorithm
!! @{
!--------------------------------
      subroutine volc_update
!................................
!  ---  inputs:    (in scope variables, module variables)
!  ---  outputs:   (module variables)

!  ==================================================================  !
!                                                                      !
!  subprogram : volc_update                                            !
!                                                                      !
!    searches historical volcanic data sets to find and read in        !
!    monthly 45-degree lat-zone band data of optical depth.            !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  inputs:  (in-scope variables, module constants)                     !
!   iyear    - integer, 4-digit calender year                 1        !
!   imon     - integer, month of the year                     1        !
!   me       - integer, print message control flag            1        !
!   NIAERCM  - integer, unit number for input data set        1        !
!                                                                      !
!  outputs: (module variables)                                         !
!   ivolae   - integer, monthly, 45-deg lat-zone volc odp      12*4*10 !
!   kyrstr   - integer, starting year of data in the input file        !
!   kyrend   - integer, ending   year of data in the input file        !
!   kyrsav   - integer, the year of data in use in the input file      !
!   kmonsav  - integer, the month of data in use in the input file     !
!                                                                      !
!  subroutines called: none                                            !
!                                                                      !
!  usage:    call volc_aerinit                                         !
!                                                                      !
!  ==================================================================  !

!  ---  inputs: (in-scope variables, module constants)
!     integer :: iyear, imon, me, NIAERCM

!  ---  output: (module variables)
!     integer :: ivolae(:,:,:), kyrstr, kyrend, kyrsav, kmonsav

!  ---  locals:
      integer :: i, j, k
      logical :: file_exist

      character :: cline*80, volcano_file*32
      data volcano_file / 'volcanic_aerosols_1850-1859.txt ' /
!
!===>  ...  begin here
!
      kmonsav = imon

      if ( kyrstr<=iyear .and. iyear<=kyrend ) then   ! use previously input data
        kyrsav = iyear
        return
      else                                            ! need to input new data
        kyrsav = iyear
        kyrstr = iyear - mod(iyear,10)
        kyrend = kyrstr + 9

!  ---  check print
!       print *,'  kyrstr, kyrend, kyrsav, kmonsav =',                  &
!    &          kyrstr,kyrend,kyrsav,kmonsav

        if ( iyear < MINVYR .or. iyear > MAXVYR ) then
!         if ( .not. allocated(ivolae) ) then
!           allocate ( ivolae(12,4,10) )   ! for 12-mon,4-lat_zone,10-year
!         endif
          ivolae(:,:,:) = 1            ! set as lowest value
          if ( me == 0 ) then
            print *,'   Request volcanic date out of range,',           &
     &              ' optical depth set to lowest value'
          endif
        else
          write(volcano_file(19:27),60) kyrstr,kyrend
  60      format(i4.4,'-',i4.4)

          inquire (file=volcano_file, exist=file_exist)
          if ( file_exist ) then
            close(NIAERCM)
            open (unit=NIAERCM,file=volcano_file,status='OLD',          &
     &            action='read',form='FORMATTED')

            read(NIAERCM,62) cline
  62        format(a80)

!  ---  check print
            if ( me == 0 ) then
              print *,'   Opened volcanic data file: ',volcano_file
              print *, cline
            endif

            do k = 1, 10
              do j = 1, 4
                read(NIAERCM,64) (ivolae(i,j,k),i=1,12)
  64            format(12i5)
              enddo
            enddo

            close (NIAERCM)
          else
            print *,'   Requested volcanic data file "',                &
     &              volcano_file,'" not found!'
            print *,'   *** Stopped in subroutine VOLC_AERINIT !!'
            stop
          endif   ! end if_file_exist_block

        endif   ! end if_iyear_block
      endif   ! end if_kyrstr_block

!  ---  check print
      if ( me == 0 ) then
        k = mod(kyrsav,10) + 1
        print *,' CHECK: Sample Volcanic data used for month, year:',   &
     &           imon, iyear
        print *,  ivolae(kmonsav,:,k)
      endif
!
      return
!................................
      end subroutine volc_update
!! @}
!--------------------------------
!
!...................................
      end subroutine aer_update
!-----------------------------------
!! @}


!>\ingroup module_radiation_aerosols
!> This subroutine computes aerosols optical properties.
!>\param prsi    (IMAX,NLP1), pressure at interface in mb
!!\param prsl    (IMAX,NLAY), layer mean pressure in mb
!!\param prslk   (IMAX,NLAY), exner function = \f$(p/p0)^{rocp}\f$
!!\param tvly    (IMAX,NLAY), layer virtual temperature in K
!!\param rhlay   (IMAX,NLAY), layer mean relative humidity
!!\param slmsk   (IMAX), sea/land mask (sea:0,land:1,sea-ice:2)
!!\param tracer  (IMAX,NLAY,NTRAC), aerosol tracer concentration
!!\param xlon    (IMAX), longitude of given points in radiance, ok for
!!               both 0->2pi or -pi->+pi ranges
!!\param xlat    (IMAX), latitude of given points in radiance, default
!!               to pi/2 -> -pi/2, otherwise see in-line comment
!!\param IMAX           1, horizontal dimension of arrays
!!\param NLAY,NLP1      1, vertical dimensions of arrays
!!\param lsswr,lslwr    logical flags for sw/lw radiation calls
!!\param aerosw    (IMAX,NLAY,NBDSW,NF_AESW), aeros opt properties for sw
!!\n                    (:,:,:,1): optical depth
!!\n                    (:,:,:,2): single scattering albedo
!!\n                    (:,:,:,3): asymmetry parameter
!!\param aerolw    (IMAX,NLAY,NBDLW,NF_AELW), aeros opt properties for lw
!!\n                    (:,:,:,1): optical depth
!!\n                    (:,:,:,2): single scattering albedo
!!\n                    (:,:,:,3): asymmetry parameter
!!\param aerodp    (IMAX,NSPC1), vertically integrated optical depth
!>\section general_setaer setaer General Algorithm
!> @{
!-----------------------------------
      subroutine setaer                                                 &
     &     ( prsi,prsl,prslk,tvly,rhlay,slmsk,tracer,xlon,xlat,         &   !  ---  inputs
     &       IMAX,NLAY,NLP1, lsswr,lslwr,                               &
     &       aerosw,aerolw                                              &   !  ---  outputs
     &,      aerodp                                                     &
     &     )

!  ==================================================================  !
!                                                                      !
!  setaer computes aerosols optical properties                         !
!                                                                      !
!  inputs:                                                   size      !
!     prsi    - pressure at interface              mb      IMAX*NLP1   !
!     prsl    - layer mean pressure                mb      IMAX*NLAY   !
!     prslk   - exner function = (p/p0)**rocp              IMAX*NLAY   !
!     tvly    - layer virtual temperature          k       IMAX*NLAY   !
!     rhlay   - layer mean relative humidity               IMAX*NLAY   !
!     slmsk   - sea/land mask (sea:0,land:1,sea-ice:2)       IMAX      !
!     tracer  - aerosol tracer concentration           IMAX*NLAY*NTRAC !
!     xlon    - longitude of given points in radiance        IMAX      !
!               ok for both 0->2pi or -pi->+pi ranges                  !
!     xlat    - latitude of given points in radiance         IMAX      !
!               default to pi/2 -> -pi/2, otherwise see in-line comment!
!     IMAX    - horizontal dimension of arrays                  1      !
!     NLAY,NLP1-vertical dimensions of arrays                   1      !
!     lsswr,lslwr                                                      !
!             - logical flags for sw/lw radiation calls         1      !
!                                                                      !
!  outputs:                                                            !
!     aerosw - aeros opt properties for sw      IMAX*NLAY*NBDSW*NF_AESW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!     aerolw - aeros opt properties for lw      IMAX*NLAY*NBDLW*NF_AELW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!     tau_gocart - 550nm aeros opt depth     IMAX*NLAY*MAX_NUM_GRIDCOMP!
!!    aerodp - vertically integrated optical depth         IMAX*NSPC1  !
!                                                                      !
!  external module variable: (in physparam)                            !
!     iaerflg - aerosol effect control flag (volc,lw,sw, 3-dig)        !
!     laswflg - tropospheric aerosol control flag for sw radiation     !
!               =f: no sw aeros calc.  =t: do sw aeros calc.           !
!     lalwflg - tropospheric aerosol control flag for lw radiation     !
!               =f: no lw aeros calc.  =t: do lw aeros calc.           !
!     lavoflg - control flag for stratospheric vocanic aerosols        !
!               =t: add volcanic aerosols to the background aerosols   !
!     ivflip  - control flag for direction of vertical index           !
!               =0: index from toa to surface                          !
!               =1: index from surface to toa                          !
!                                                                      !
!  internal module variable: (set by subroutine aer_init)              !
!     ivolae  - stratosphere volcanic aerosol optical depth (fac 1.e4) !
!                                                     12*4*10          !
!  usage:    call setaer                                               !
!                                                                      !
!  subprograms called:  aer_property                                   !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      integer, intent(in) :: IMAX, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: prsi, prsl,  &
     &       prslk, tvly, rhlay
      real (kind=kind_phys), dimension(:),   intent(in) :: xlon, xlat,  &
     &       slmsk
      real (kind=kind_phys), dimension(:,:,:),intent(in):: tracer

      logical, intent(in) :: lsswr, lslwr


!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:,:), intent(out) ::         &
     &       aerosw, aerolw

      real (kind=kind_phys), dimension(:,:)    , intent(out) :: aerodp

!  ---  locals:
      real (kind=kind_phys), parameter :: psrfh = 5.0    ! ref press (mb) for upper bound

      real (kind=kind_phys), dimension(IMAX) :: alon,alat,volcae,rdelp
!     real (kind=kind_phys), dimension(IMAX) :: sumodp
      real (kind=kind_phys) :: prsln(NLP1),hz(IMAX,NLP1),dz(IMAX,NLAY)
      real (kind=kind_phys) :: tmp1, tmp2, psrfl

      integer               :: kcutl(IMAX), kcuth(IMAX)
      integer               :: i, i1, j, k, m, mb, kh, kl

      logical               :: laddsw=.false.,  laersw=.false.
      logical               :: laddlw=.false.,  laerlw=.false.

!  ---  conversion constants
      real (kind=kind_phys), parameter :: rdg  = 180.0 / con_pi
      real (kind=kind_phys), parameter :: rovg = 0.001 * con_rd / con_g

!===>  ...  begin here

      do m = 1, NF_AESW
        do j = 1, NBDSW
          do k = 1, NLAY
            do i = 1, IMAX
              aerosw(i,k,j,m) = f_zero
            enddo
          enddo
        enddo
      enddo

      do m = 1, NF_AELW
        do j = 1, NBDLW
          do k = 1, NLAY
            do i = 1, IMAX
              aerolw(i,k,j,m) = f_zero
            enddo
          enddo
        enddo
      enddo

!     sumodp = f_zero
      do i = 1, IMAX
       do k = 1, NSPC1
         aerodp(i,k) = f_zero
       enddo
      enddo


      if ( .not. (lsswr .or. lslwr) ) then
        return
      endif

      if ( iaerflg <= 0 ) then
        return
      endif

      laersw = lsswr .and. laswflg
      laerlw = lslwr .and. lalwflg

!> -# Convert lat/lon from radiance to degree.

      do i = 1, IMAX
        alon(i) = xlon(i) * rdg
        if (alon(i) < f_zero) alon(i) = alon(i) + 360.0
        alat(i) = xlat(i) * rdg          ! if xlat in pi/2 -> -pi/2 range
!       alat(i) = 90.0 - xlat(i)*rdg     ! if xlat in 0 -> pi range
      enddo

!> -# Compute level height and layer thickness.

      if ( laswflg .or. lalwflg ) then

        lab_do_IMAX : do i = 1, IMAX

          lab_if_flip : if (ivflip == 1) then       ! input from sfc to toa

            do k = 1, NLAY
              prsln(k) = log(prsi(i,k))
            enddo
            prsln(NLP1)= log(prsl(i,NLAY))

            do k = NLAY, 1, -1
              dz(i,k) = rovg * (prsln(k) - prsln(k+1)) * tvly(i,k)
            enddo
            dz(i,NLAY)  = 2.0 * dz(i,NLAY)

            hz(i,1) = f_zero
            do k = 1, NLAY
              hz(i,k+1) = hz(i,k) + dz(i,k)
            enddo

          else  lab_if_flip                         ! input from toa to sfc

            prsln(1) = log(prsl(i,1))
            do k = 2, NLP1
              prsln(k) = log(prsi(i,k))
            enddo

            do k = 1, NLAY
              dz(i,k) = rovg * (prsln(k+1) - prsln(k)) * tvly(i,k)
            enddo
            dz(i,1) = 2.0 * dz(i,1)

            hz(i,NLP1) = f_zero
            do k = NLAY, 1, -1
              hz(i,k) = hz(i,k+1) + dz(i,k)
            enddo

          endif  lab_if_flip

        enddo  lab_do_IMAX


!> -# Calculate SW aerosol optical properties for the corresponding
!!    frequency bands:
!!    - if opac aerosol climatology is used, call aer_property(): this
!!      subroutine maps the 5 degree global climatological aerosol data
!!      set onto model grids, and compute aerosol optical properties for
!!      SW and LW radiations.
!!    - if gocart aerosol scheme is used, call setgocartaer(): this
!!      subroutine computes sw + lw aerosol optical properties for gocart
!!      aerosol species (merged from fcst and clim fields).

!SARAH
!         if ( iaerflg == 1 ) then      ! use opac aerosol climatology
          if ( iaermdl == 0 ) then      ! use opac aerosol climatology

          call aer_property                                               &
!  ---  inputs:
     &       ( prsi,prsl,prslk,tvly,rhlay,dz,hz,tracer,                   &
     &         alon,alat,slmsk, laersw,laerlw,                            &
     &         IMAX,NLAY,NLP1,                                            &
!    &         IMAX,NLAY,NLP1,NSPC1,                                      &
!  ---  outputs:
     &         aerosw,aerolw,aerodp                                       &
     &       )

!  ---  check print
!       do m = 1, NBDSW
!         print *,'  ***  CHECK AEROSOLS PROPERTIES FOR SW BAND =',m,   &
!    &            ' ***'
!         do k = 1, 10
!           print *,'  LEVEL :',k
!           print *,'  TAUAER:',aerosw(:,k,m,1)
!           print *,'  SSAAER:',aerosw(:,k,m,2)
!           print *,'  ASYAER:',aerosw(:,k,m,3)
!         enddo
!       enddo
!       print *,'  ***  CHECK AEROSOLS OPTICAL DEPTH FOR 550nm REGION'
!       print *, aerodp(:,1)
!       if ( laod_out ) then
!         do m = 1, NSPC1
!           print *,'  ***  CHECK AEROSOLS OPTICAL DEPTH FOR SPECIES:', &
!    &              m
!           print *, aerodp(:,m)
!           sumodp(:) = sumodp(:) + aerodp(:,m)
!         enddo

!
!         print *,'  ***  CHECK AEROSOLS OPTICAL DEPTH FOR ALL SPECIES:'
!         print *, sumodp(:)
!       endif
!       do m = 1, NBDLW
!         print *,'  ***  CHECK AEROSOLS PROPERTIES FOR LW BAND =',m,   &
!    &            ' ***'
!         do k = 1, 10
!           print *,'  LEVEL :',k
!           print *,'  TAUAER:',aerolw(:,k,m,1)
!           print *,'  SSAAER:',aerolw(:,k,m,2)
!           print *,'  ASYAER:',aerolw(:,k,m,3)
!         enddo
!       enddo
! SARAH
!        elseif ( iaerflg == 2 )   then    ! use gocart aerosol scheme
         elseif ( iaermdl == 1 )   then    ! use gocart aerosol scheme

          call setgocartaer                                               &

!  ---  inputs:
     &       ( alon,alat,prslk,rhlay,dz,hz,NSWLWBD,                       &
     &         prsl,tvly,tracer,                                          &
     &         IMAX,NLAY,NLP1, ivflip, lsswr,lslwr,                       &
!  ---  outputs:
     &         aerosw,aerolw                                              &
     &     )

        endif     ! end if_iaerflg_block

      endif   ! end if_laswflg_or_lalwflg_block

!> -# Compute stratosphere volcanic forcing:
!!    - select data in 4 lat bands, interpolation at the boundaries
!!    - Find lower boundary of stratosphere: polar, fixed at 25000pa
!!      (250mb); tropic, fixed at 15000pa (150mb); mid-lat, interpolation
!!    - SW: add volcanic aerosol optical depth to the background value
!!    - Smoothing profile at boundary if needed
!!    - LW: add volcanic aerosol optical depth to the background value
!  ---  ...  stratosphere volcanic forcing

      if ( lavoflg ) then

        if ( iaerflg == 100 ) then
          laddsw = lsswr
          laddlw = lslwr
        else
          laddsw = lsswr .and. laswflg
          laddlw = lslwr .and. lalwflg
        endif

        i1 = mod(kyrsav, 10) + 1

!  ---  select data in 4 lat bands, interpolation at the boundaires

        do i = 1, IMAX
          if      ( alat(i) > 46.0 ) then
            volcae(i) = 1.0e-4 * ivolae(kmonsav,1,i1)
          else if ( alat(i) > 44.0 ) then
            volcae(i) = 5.0e-5                                          &
     &                * (ivolae(kmonsav,1,i1) + ivolae(kmonsav,2,i1))
          else if ( alat(i) >  1.0 ) then
            volcae(i) = 1.0e-4 * ivolae(kmonsav,2,i1)
          else if ( alat(i) > -1.0 ) then
            volcae(i) = 5.0e-5                                          &
     &                * (ivolae(kmonsav,2,i1) + ivolae(kmonsav,3,i1))
          else if ( alat(i) >-44.0 ) then
            volcae(i) = 1.0e-4 * ivolae(kmonsav,3,i1)
          else if ( alat(i) >-46.0 ) then
            volcae(i) = 5.0e-5                                          &
     &                * (ivolae(kmonsav,3,i1) + ivolae(kmonsav,4,i1))
          else
            volcae(i) = 1.0e-4 * ivolae(kmonsav,4,i1)
          endif
        enddo

        if ( ivflip == 0 ) then         ! input data from toa to sfc

!  ---  find lower boundary of stratosphere

          do i = 1, IMAX

            tmp1 = abs( alat(i) )
            if ( tmp1 > 70.0 ) then          ! polar, fixed at 25000pa (250mb)
              psrfl = 250.0
            elseif ( tmp1 < 20.0 ) then      ! tropic, fixed at 15000pa (150mb)
              psrfl = 150.0
            else                             ! mid-lat, interpolation
              psrfl = 110.0 + 2.0*tmp1
            endif

            kcuth(i) = NLAY - 1
            kcutl(i) = 2
            rdelp(i) = f_one / prsi(i,2)

            lab_do_kcuth0 : do k = 2, NLAY-2
              if ( prsi(i,k) >= psrfh ) then
                kcuth(i) = k - 1
                exit lab_do_kcuth0
              endif
            enddo  lab_do_kcuth0

            lab_do_kcutl0 : do k = 2, NLAY-2
              if ( prsi(i,k) >= psrfl ) then
                kcutl(i) = k - 1
                rdelp(i) = f_one / (prsi(i,k) - prsi(i,kcuth(i)))
                exit lab_do_kcutl0
              endif
            enddo  lab_do_kcutl0
          enddo

!  ---  sw: add volcanic aerosol optical depth to the background value

          if ( laddsw ) then
            do m = 1, NBDSW
              mb = NSWSTR + m - 1

              if     ( wvnsw1(mb) > 20000 ) then   ! range of wvlth < 0.5mu
                tmp2 = 0.74
              elseif ( wvnsw2(mb) < 20000 ) then   ! range of wvlth > 0.5mu
                tmp2 = 1.14
              else                                 ! range of wvlth in btwn
                tmp2 = 0.94
              endif
              tmp1 = (0.275e-4 * (wvnsw2(mb)+wvnsw1(mb))) ** tmp2

              do i = 1, IMAX
                kh = kcuth(i)
                kl = kcutl(i)
                do k = kh, kl
                  tmp2 = tmp1 * ((prsi(i,k+1) - prsi(i,k)) * rdelp(i))
                  aerosw(i,k,m,1) = aerosw(i,k,m,1) + tmp2*volcae(i)
                enddo

!  ---  smoothing profile at boundary if needed

                if ( aerosw(i,kl,m,1) > 10.*aerosw(i,kl+1,m,1) ) then
                  tmp2 = aerosw(i,kl,m,1) + aerosw(i,kl+1,m,1)
                  aerosw(i,kl  ,m,1) = 0.8 * tmp2
                  aerosw(i,kl+1,m,1) = 0.2 * tmp2
                endif
              enddo    ! end do_i_block
            enddo      ! end do_m_block

!  ---  check print
!           do i = 1, IMAX
!             print *,' LEV  PRESS      TAU      FOR PROFILE:',i,       &
!    &                '  KCUTH, KCUTL =',kcuth(i),kcutl(i)
!             kh = kcuth(i) - 1
!             kl = kcutl(i) + 10
!             do k = kh, kl
!               write(6,71) k, prsl(i,k), aerosw(i,k,1,1)
! 71            format(i3,2e11.4)
!             enddo
!           enddo

          endif        ! end if_laddsw_block

!  ---  lw: add volcanic aerosol optical depth to the background value

          if ( laddlw ) then
            if ( NLWBND == 1 ) then

              tmp1 = (0.55 / 11.0) ** 1.2
              do i = 1, IMAX
                kh = kcuth(i)
                kl = kcutl(i)
                do k = kh, kl
                  tmp2 = tmp1 * ((prsi(i,k+1) - prsi(i,k)) * rdelp(i))  &
     &                 * volcae(i)
                  do m = 1, NBDLW
                    aerolw(i,k,m,1) = aerolw(i,k,m,1) + tmp2
                  enddo
                enddo
              enddo    ! end do_i_block

            else

              do m = 1, NBDLW
                tmp1 = (0.275e-4 * (wvnlw2(m) + wvnlw1(m))) ** 1.2

                do i = 1, IMAX
                  kh = kcuth(i)
                  kl = kcutl(i)
                  do k = kh, kl
                    tmp2 = tmp1 * ((prsi(i,k+1)-prsi(i,k)) * rdelp(i))
                    aerolw(i,k,m,1) = aerolw(i,k,m,1) + tmp2*volcae(i)
                  enddo
                enddo    ! end do_i_block
              enddo      ! end do_m_block

            endif      ! end if_NLWBND_block
          endif        ! end if_laddlw_block

        else                            ! input data from sfc to toa

!  ---  find lower boundary of stratosphere

          do i = 1, IMAX

            tmp1 = abs( alat(i) )
            if ( tmp1 > 70.0 ) then          ! polar, fixed at 25000pa (250mb)
              psrfl = 250.0
            elseif ( tmp1 < 20.0 ) then      ! tropic, fixed at 15000pa (150mb)
              psrfl = 150.0
            else                             ! mid-lat, interpolation
              psrfl = 110.0 + 2.0*tmp1
            endif

            kcuth(i) = 2
            kcutl(i) = NLAY - 1
            rdelp(i) = f_one / prsi(i,NLAY-1)

            lab_do_kcuth1 : do k = NLAY-1, 2, -1
              if ( prsi(i,k) >= psrfh ) then
                kcuth(i) = k
                exit lab_do_kcuth1
              endif
            enddo  lab_do_kcuth1

            lab_do_kcutl1 : do k = NLAY, 2, -1
              if ( prsi(i,k) >= psrfl ) then
                kcutl(i) = k
                rdelp(i) = f_one / (prsi(i,k) - prsi(i,kcuth(i)+1))
                exit lab_do_kcutl1
              endif
            enddo  lab_do_kcutl1
          enddo

!  ---  sw: add volcanic aerosol optical depth to the background value

          if ( laddsw ) then
            do m = 1, NBDSW
              mb = NSWSTR + m - 1

              if     ( wvnsw1(mb) > 20000 ) then   ! range of wvlth < 0.5mu
                tmp2 = 0.74
              elseif ( wvnsw2(mb) < 20000 ) then   ! range of wvlth > 0.5mu
                tmp2 = 1.14
              else                                 ! range of wvlth in btwn
                tmp2 = 0.94
              endif
              tmp1 = (0.275e-4 * (wvnsw2(mb)+wvnsw1(mb))) ** tmp2

              do i = 1, IMAX
                kh = kcuth(i)
                kl = kcutl(i)
                do k = kl, kh
                  tmp2 = tmp1 * ((prsi(i,k) - prsi(i,k+1)) * rdelp(i))
                  aerosw(i,k,m,1) = aerosw(i,k,m,1) + tmp2*volcae(i)
                enddo

!  ---  smoothing profile at boundary if needed

                if ( aerosw(i,kl,m,1) > 10.*aerosw(i,kl-1,m,1) ) then
                  tmp2 = aerosw(i,kl,m,1) + aerosw(i,kl-1,m,1)
                  aerosw(i,kl  ,m,1) = 0.8 * tmp2
                  aerosw(i,kl-1,m,1) = 0.2 * tmp2
                endif
              enddo    ! end do_i_block
            enddo      ! end do_m_block

!  ---  check print
!           do i = 1, IMAX
!             print *,' LEV  PRESS      TAU      FOR PROFILE:',i,       &
!    &                '  KCUTH, KCUTL =',kcuth(i),kcutl(i)
!             kh = kcuth(i) + 1
!             kl = kcutl(i) - 10
!             do k = kh, kl, -1
!               write(6,71) NLP1-k,prsl(i,k),aerosw(i,k,1,1)
!             enddo
!           enddo

          endif        ! end if_laddsw_block

!  ---  lw: add volcanic aerosol optical depth to the background value

          if ( laddlw ) then
            if ( NLWBND == 1 ) then

              tmp1 = (0.55 / 11.0) ** 1.2
              do i = 1, IMAX
                kh = kcuth(i)
                kl = kcutl(i)
                do k = kl, kh
                  tmp2 = tmp1 * ((prsi(i,k) - prsi(i,k+1)) * rdelp(i))  &
     &                 * volcae(i)
                  do m = 1, NBDLW
                    aerolw(i,k,m,1) = aerolw(i,k,m,1) + tmp2
                  enddo
                enddo
              enddo    ! end do_i_block

            else

              do m = 1, NBDLW
                tmp1 = (0.275e-4 * (wvnlw2(m) + wvnlw1(m))) ** 1.2

                do i = 1, IMAX
                  kh = kcuth(i)
                  kl = kcutl(i)
                  do k = kl, kh
                    tmp2 = tmp1 * ((prsi(i,k)-prsi(i,k+1)) * rdelp(i))
                    aerolw(i,k,m,1) = aerolw(i,k,m,1) + tmp2*volcae(i)
                  enddo
                enddo    ! end do_i_block
              enddo      ! end do_m_block

            endif      ! end if_NLWBND_block
          endif        ! end if_laddlw_block

        endif                           ! end if_ivflip_block

      endif   ! end if_lavoflg_block
!
      return
!...................................
      end subroutine setaer
!-----------------------------------
!> @}


!>\ingroup module_radiation_aerosols
!> This subroutine maps the 5 degree global climatological aerosol data
!! set onto model grids, and compute aerosol optical properties for SW
!! and LW radiations.
!!\param prsi           (IMAX,NLP1), pressure at interface in mb
!!\param prsl           (IMAX,NLAY), layer mean pressure(not used)
!!\param prslk          (IMAX,NLAY), exner function=\f$(p/p0)^{rocp}\f$ (not used)
!!\param tvly           (IMAX,NLAY), layer virtual temperature (not used)
!!\param rhlay          (IMAX,NLAY), layer mean relative humidity
!!\param dz             (IMAX,NLAY), layer thickness in m
!!\param hz             (IMAX,NLP1), level high in m
!!\param tracer         (IMAX,NLAY,NTRAC), aer tracer concentrations (not used)
!!\param alon, alat     (IMAX), longitude and latitude of given points in degree
!!\param slmsk          (IMAX), sea/land mask (sea:0,land:1,sea-ice:2)
!!\param laersw,laerlw  logical flag for sw/lw aerosol calculations
!!\param IMAX           horizontal dimension of arrays
!!\param NLAY,NLP1      vertical dimensions of arrays
!!\param NSPC           num of species for optional aod output fields
!!\param aerosw        (IMAX,NLAY,NBDSW,NF_AESW), aeros opt properties for sw
!!\n                              (:,:,:,1): optical depth
!!\n                              (:,:,:,2): single scattering albedo
!!\n                              (:,:,:,3): asymmetry parameter
!!\param aerolw        (IMAX,NLAY,NBDLW,NF_AELW), aeros opt properties for lw
!!\n                              (:,:,:,1): optical depth
!!\n                              (:,:,:,2): single scattering albedo
!!\n                              (:,:,:,3): asymmetry parameter
!!\param aerodp        (IMAX,NSPC+1), vertically integrated aer-opt-depth
!!\section gel_aer_pro aer_property General Algorithm
!> @{
!-----------------------------------
      subroutine aer_property                                           &
     &     ( prsi,prsl,prslk,tvly,rhlay,dz,hz,tracer,                   &     !  ---  inputs:
     &       alon,alat,slmsk, laersw,laerlw,                            &
     &       IMAX,NLAY,NLP1,                                            &
     &       aerosw,aerolw,aerodp                                       &     !  ---  outputs:
     &     )

!  ==================================================================  !
!                                                                      !
!  aer_property maps the 5 degree global climatological aerosol data   !
!  set onto model grids, and compute aerosol optical properties for sw !
!  and lw radiations.                                                  !
!                                                                      !
!  inputs:                                                             !
!     prsi    - pressure at interface              mb      IMAX*NLP1   !
!     prsl    - layer mean pressure         (not used)     IMAX*NLAY   !
!     prslk   - exner function=(p/p0)**rocp (not used)     IMAX*NLAY   !
!     tvly    - layer virtual temperature   (not used)     IMAX*NLAY   !
!     rhlay   - layer mean relative humidity               IMAX*NLAY   !
!     dz      - layer thickness                    m       IMAX*NLAY   !
!     hz      - level high                         m       IMAX*NLP1   !
!     tracer  - aer tracer concentrations   (not used)  IMAX*NLAY*NTRAC!
!     alon, alat                                             IMAX      !
!             - longitude and latitude of given points in degree       !
!     slmsk   - sea/land mask (sea:0,land:1,sea-ice:2)       IMAX      !
!     laersw,laerlw                                             1      !
!             - logical flag for sw/lw aerosol calculations            !
!     IMAX    - horizontal dimension of arrays                  1      !
!     NLAY,NLP1-vertical dimensions of arrays                   1      !
!!    NSPC    - num of species for optional aod output fields   1      !
!                                                                      !
!  outputs:                                                            !
!     aerosw - aeros opt properties for sw      IMAX*NLAY*NBDSW*NF_AESW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!     aerolw - aeros opt properties for lw      IMAX*NLAY*NBDLW*NF_AELW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!!    aerodp - vertically integrated aer-opt-depth         IMAX*NSPC+1 !
!                                                                      !
!  module parameters and constants:                                    !
!     NSWBND  - total number of actual sw spectral bands computed      !
!     NLWBND  - total number of actual lw spectral bands computed      !
!     NSWLWBD - total number of sw+lw bands computed                   !
!                                                                      !
!  external module variables: (in physparam)                           !
!     ivflip  - control flag for direction of vertical index           !
!               =0: index from toa to surface                          !
!               =1: index from surface to toa                          !
!                                                                      !
!  module variable: (set by subroutine aer_init)                       !
!     kprfg   - aerosols profile index                IMXAE*JMXAE      !
!               1:ant  2:arc  3:cnt  4:mar  5:des  6:marme 7:cntme     !
!     idxcg   - aerosols component index              NXC*IMXAE*JMXAE  !
!               1:inso    2:soot    3:minm    4:miam    5:micm         !
!               6:mitr    7:waso    8:ssam    9:sscm   10:suso         !
!     cmixg   - aerosols component mixing ratio       NXC*IMXAE*JMXAE  !
!     denng   - aerosols number density                2 *IMXAE*JMXAE  !
!               1:for domain-1   2:domain-2 (prof marme/cntme only)    !
!                                                                      !
!  usage:    call aer_property                                         !
!                                                                      !
!  subprograms called:  radclimaer                                     !
!                                                                      !
!  ==================================================================  !

!  ---  inputs:
      integer, intent(in) :: IMAX, NLAY, NLP1
!     integer, intent(in) :: IMAX, NLAY, NLP1, NSPC
      logical, intent(in) :: laersw, laerlw

      real (kind=kind_phys), dimension(:,:), intent(in) :: prsi, prsl,  &
     &       prslk, tvly, rhlay, dz, hz
      real (kind=kind_phys), dimension(:),   intent(in) :: alon, alat,  &
     &       slmsk
      real (kind=kind_phys), dimension(:,:,:),intent(in):: tracer

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:,:), intent(out) ::         &
     &       aerosw, aerolw
      real (kind=kind_phys), dimension(:,:)    , intent(out) :: aerodp

!  ---  locals:
      real (kind=kind_phys), dimension(NCM) :: cmix
      real (kind=kind_phys), dimension(  2) :: denn
      real (kind=kind_phys), dimension(NSPC) :: spcodp

      real (kind=kind_phys), dimension(NLAY) :: delz, rh1, dz1
      integer,               dimension(NLAY) :: idmaer

      real (kind=kind_phys), dimension(NLAY,NSWLWBD):: tauae,ssaae,asyae
!test real (kind=kind_phys), dimension(IMAX,NLAY) :: aersav

      real (kind=kind_phys) :: tmp1, tmp2, rps, dtmp, h1
      real (kind=kind_phys) :: wi, wj, w11, w12, w21, w22

      integer :: i, ii, i1, i2, i3,  j1, j2, j3,  k, m, m1,             &
     &           kp, kpa, kpi, kpj

!  ---  conversion constants
      real (kind=kind_phys), parameter :: dltg = 360.0 / float(IMXAE)
      real (kind=kind_phys), parameter :: hdlt = 0.5 * dltg
      real (kind=kind_phys), parameter :: rdlt = 1.0 / dltg

!
!===>  ...  begin here
!
!> -# Map aerosol data to model grids
!!    - Map grid in longitude direction, lon from 0 to 355 deg resolution
!!    - Map grid in latitude direction, lat from 90n to 90s in 5 deg resolution

      i1 = 1
      i2 = 2
      j1 = 1
      j2 = 2

      lab_do_IMAX : do i = 1, IMAX

!  ---  map grid in longitude direction, lon from 0 to 355 deg resolution

!       print *,' Seeking lon index for point i =',i
        i3 = i1
        lab_do_IMXAE : do while ( i3 <= IMXAE )
          tmp1 = dltg * (i3 - 1)
          dtmp = alon(i) - tmp1
!         print *,'   alon, i3, tlon, dlon =',alon(i),i3,tmp1,dtmp

          if ( dtmp > dltg ) then
            i3 = i3 + 1
            if ( i3 > IMXAE ) then
              print *,' ERROR! In setclimaer alon>360. ipt =',i,        &
     &           ',  dltg,alon,tlon,dlon =',dltg,alon(i),tmp1,dtmp
              stop
            endif
          elseif ( dtmp >= f_zero ) then
            i1 = i3
            i2 = mod(i3,IMXAE) + 1
            wi = dtmp * rdlt
            if ( dtmp <= hdlt ) then
              kpi = i3
            else
              kpi = i2
            endif
!           print *,'   found i1, i2, wi =',i1,i2,wi
            exit lab_do_IMXAE
          else
            i3 = i3 - 1
            if ( i3 < 1 ) then
              print *,' ERROR! In setclimaer alon< 0. ipt =',i,         &
     &           ',  dltg,alon,tlon,dlon =',dltg,alon(i),tmp1,dtmp
              stop
            endif
          endif
        enddo  lab_do_IMXAE

!  ---  map grid in latitude direction, lat from 90n to 90s in 5 deg resolution

!       print *,' Seeking lat index for point i =',i
        j3 = j1
        lab_do_JMXAE : do while ( j3 <= JMXAE )
          tmp2 = 90.0 - dltg * (j3 - 1)
          dtmp = tmp2 - alat(i)
!         print *,'   alat, j3, tlat, dlat =',alat(i),j3,tmp2,dtmp

          if ( dtmp > dltg ) then
            j3 = j3 + 1
            if ( j3 >= JMXAE ) then
              print *,' ERROR! In setclimaer alat<-90. ipt =',i,        &
     &           ',  dltg,alat,tlat,dlat =',dltg,alat(i),tmp2,dtmp
              stop
            endif
          elseif ( dtmp >= f_zero ) then
            j1 = j3
            j2 = j3 + 1
            wj = dtmp * rdlt
            if ( dtmp <= hdlt ) then
              kpj = j3
            else
              kpj = j2
            endif
!           print *,'   found j1, j2, wj =',j1,j2,wj
            exit lab_do_JMXAE
          else
            j3 = j3 - 1
            if ( j3 < 1 ) then
              print *,' ERROR! In setclimaer alat>90. ipt =',i,         &
     &           ',  dltg,alat,tlat,dlat =',dltg,alat(i),tmp2,dtmp
              stop
            endif
          endif
        enddo  lab_do_JMXAE

!> -# Determin the type of aerosol profile (kp) and scale hight for
!!    domain 1 (h1) to be used at this grid point.

        kp = kprfg(kpi,kpj)                     ! nearest typical aeros profile as default
        kpa = max( kprfg(i1,j1),kprfg(i1,j2),kprfg(i2,j1),kprfg(i2,j2) )
        h1 = haer(1,kp)
        denn(2) = f_zero
        ii = 1

        if ( kp /= kpa ) then
          if ( kpa == 6 ) then                  ! if ocean prof with mineral aeros overlay
            ii = 2                              ! need 2 types of densities
            if ( slmsk(i) > f_zero ) then       ! but actually a land/sea-ice point
              kp = 7                            ! reset prof index to land
              h1 = 0.5*(haer(1,6) + haer(1,7))  ! use a transition scale hight
            else
              kp = kpa
              h1 = haer(1,6)
            endif
          elseif ( kpa == 7 ) then              ! if land prof with mineral aeros overlay
            ii = 2                              ! need 2 types of densities
            if ( slmsk(i) <= f_zero ) then      ! but actually an ocean point
              kp = 6                            ! reset prof index to ocean
              h1 = 0.5*(haer(1,6) + haer(1,7))  ! use a transition scale hight
            else
              kp = kpa
              h1 = haer(1,7)
            endif
          else                                  ! lower atmos without mineral aeros overlay
!           h1 = 0.5*(haer(1,kp) + haer(1,kpa)) ! use a transition scale hight
            h1 = haer(1,kpa)
            kp = kpa
          endif
        endif

!> -# Compute horizontal bi-linear interpolation weights

        w11 = (f_one-wi) * (f_one-wj)
        w12 = (f_one-wi) *       wj
        w21 =        wi  * (f_one-wj)
        w22 =        wi  * wj

!  ---  check print
!       print *,'  Grid pt', i,',   alon, alat =',alon(i),alat(i),      &
!    &                       ',   tlon, tlat =',tmp1,tmp2
!       print *,'   lon grid index i1, i2 =',i1,i2,',  weight wi =',wi
!       print *,'   lat grid index j1, j2 =',j1,j2,',  weight wj =',wj
!       print *,'   bi-linear weights w11,w21,w12,w22 =',w11,w21,w12,w22
!       print *,'   kp,kpa,slmsk,h1 =',kp,m1,slmsk(i),h1

!> -# Do horizontal bi-linear interpolation on aerosol partical density
!!   (denn)

        do m = 1, ii                            ! ii=1 for domain 1; =2 for domain 2.
          denn(m) = w11*denng(m,i1,j1) + w12*denng(m,i1,j2)             &
     &            + w21*denng(m,i2,j1) + w22*denng(m,i2,j2)
        enddo  ! end_do_m_loop

!> -# Do horizontal bi-linear interpolation on mixing ratios

        cmix(:) = f_zero
        do m = 1, NXC
          ii = idxcg(m,i1,j1)
          if ( ii > 0 ) then
            cmix(ii) = cmix(ii) + w11*cmixg(m,i1,j1)
          endif
          ii = idxcg(m,i1,j2)
          if ( ii > 0 ) then
            cmix(ii) = cmix(ii) + w12*cmixg(m,i1,j2)
          endif
          ii = idxcg(m,i2,j1)
          if ( ii > 0 ) then
            cmix(ii) = cmix(ii) + w21*cmixg(m,i2,j1)
          endif
          ii = idxcg(m,i2,j2)
          if ( ii > 0 ) then
            cmix(ii) = cmix(ii) + w22*cmixg(m,i2,j2)
          endif
        enddo  ! end_do_m_loop

!  ---  check print
!       print *,'   denn =',denn(:)
!       print *,'   cmix =',cmix(:)

!> -# Prepare to setup domain index array and effective layer thickness,
!!    also convert pressure level to sigma level to follow the terrain.

        do k = 1, NLAY
          rh1(k) = rhlay(i,k)
          dz1(k) = dz   (i,k)
        enddo

        lab_if_flip : if (ivflip == 1) then       ! input from sfc to toa

          if ( prsi(i,1) > 100.0 ) then
            rps = f_one / prsi(i,1)
          else
            print *,' !!! (1) Error in subr radiation_aerosols:',       &
     &              ' unrealistic surface pressure =', i,prsi(i,1)
            stop
          endif

          ii = 1
          do k = 1, NLAY
            if (prsi(i,k+1)*rps < sigref(ii,kp)) then
              ii = ii + 1
              if (ii == 2 .and. prsref(2,kp) == prsref(3,kp)) then
                ii = 3
              endif
            endif
            idmaer(k) = ii

            if ( ii > 1 ) then
              tmp1 = haer(ii,kp)
            else
              tmp1 = h1
            endif

            if (tmp1 > f_zero) then
              tmp2 = f_one / tmp1
              delz(k) = tmp1 * (exp(-hz(i,k)*tmp2)-exp(-hz(i,k+1)*tmp2))
            else
              delz(k) = dz1(k)
            endif
          enddo

        else  lab_if_flip                         ! input from toa to sfc

          if ( prsi(i,NLP1) > 100.0 ) then
            rps =  1.0 / prsi(i,NLP1)
          else
            print *,' !!! (2) Error in subr radiation_aerosols:',       &
     &              ' unrealistic surface pressure =', i,prsi(i,NLP1)
          endif

          ii = 1
          do k = NLAY, 1, -1
            if (prsi(i,k)*rps < sigref(ii,kp)) then
              ii = ii + 1
              if (ii == 2 .and. prsref(2,kp) == prsref(3,kp)) then
                ii = 3
              endif
            endif
            idmaer(k) = ii

            if ( ii > 1 ) then
              tmp1 = haer(ii,kp)
            else
              tmp1 = h1
            endif

            if (tmp1 > f_zero) then
              tmp2   = f_one / tmp1
              delz(k) = tmp1 * (exp(-hz(i,k+1)*tmp2)-exp(-hz(i,k)*tmp2))
            else
              delz(k) = dz1(k)
            endif
          enddo

        endif  lab_if_flip

!  ---  check print

!       print *,' in setclimaer, profile:',i
!       print *,'  rh   :',rh1
!       print *,'  dz   :',dz1
!       print *,'  delz :',delz
!       print *,'  idmaer:',idmaer

!> -# Call radclimaer() to calculate SW/LW aerosol optical properties
!!    for the corresponding frequency bands.

        call radclimaer
!  ---  inputs:  (in-scope variables)
!  ---  outputs: (in-scope variables)

        if ( laersw ) then

          do m = 1, NBDSW
            do k = 1, NLAY
              aerosw(i,k,m,1) = tauae(k,m)
              aerosw(i,k,m,2) = ssaae(k,m)
              aerosw(i,k,m,3) = asyae(k,m)
            enddo
          enddo

!  ---  total aod (optional)
         do k = 1, NLAY
           aerodp(i,1) = aerodp(i,1) + tauae(k,nv_aod)
         enddo

!  ---  for diagnostic output (optional)
!         if ( lspcaod ) then
           do m = 1, NSPC
             aerodp(i,m+1) = spcodp(m)
           enddo
!         endif

        endif     ! end if_larsw_block

        if ( laerlw ) then

          if ( NLWBND == 1 ) then
            m1 = NSWBND + 1
            do m = 1, NBDLW
              do k = 1, NLAY
                aerolw(i,k,m,1) = tauae(k,m1)
                aerolw(i,k,m,2) = ssaae(k,m1)
                aerolw(i,k,m,3) = asyae(k,m1)
              enddo
            enddo
          else
            do m = 1, NBDLW
              m1 = NSWBND + m
              do k = 1, NLAY
                aerolw(i,k,m,1) = tauae(k,m1)
                aerolw(i,k,m,2) = ssaae(k,m1)
                aerolw(i,k,m,3) = asyae(k,m1)
              enddo
            enddo
          endif

        endif     ! end if_laerlw_block

      enddo  lab_do_IMAX

! =================
      contains
! =================

!>\ingroup module_radiation_aerosols
!> This subroutine computes aerosols optical properties in NSWLWBD
!! bands. there are seven different vertical profile structures. in the
!! troposphere, aerosol distribution at each grid point is composed
!! from up to six components out of ten different substances.
!!\section radclimaer_gen radclimaer General Algorithm
!--------------------------------
      subroutine radclimaer
!................................

!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

!  ==================================================================  !
!                                                                      !
!  compute aerosols optical properties in NSWLWBD bands. there are     !
!  seven different vertical profile structures. in the troposphere,    !
!  aerosol distribution at each grid point is composed from up to      !
!  six components out of a total of ten different substances.          !
!                                                                      !
!  ref: wmo report wcp-112 (1986)                                      !
!                                                                      !
!  input variables:                                                    !
!     cmix   - mixing ratioes of aerosol components  -     NCM         !
!     denn   - aerosol number densities              -     2           !
!     rh1    - relative humidity                     -     NLAY        !
!     delz   - effective layer thickness             km    NLAY        !
!     idmaer - aerosol domain index                  -     NLAY        !
!     NXC    - number of different aerosol components-     1           !
!     NLAY   - vertical dimensions                   -     1           !
!                                                                      !
!  output variables:                                                   !
!     tauae  - optical depth                         -     NLAY*NSWLWBD!
!     ssaae  - single scattering albedo              -     NLAY*NSWLWBD!
!     asyae  - asymmetry parameter                   -     NLAY*NSWLWBD!
!!    aerodp - vertically integrated aer-opt-depth   -     IMAX*NSPC+1 !
!                                                                      !
!  ==================================================================  !
!
      real (kind=kind_phys) :: crt1, crt2
      parameter (crt1=30.0, crt2=0.03333)

!  ---  inputs:
!  ---  outputs:

!  ---  locals:
      real (kind=kind_phys) :: cm, hd, hdi, sig0u, sig0l, ratio, tt0,   &
     &      ex00, sc00, ss00, as00, ex01, sc01, ss01, as01,     tt1,    &
     &      ex02, sc02, ss02, as02, ex03, sc03, ss03, as03,     tt2,    &
     &      ext1, sca1, ssa1, asy1, drh0, drh1, rdrh

      integer :: ih1, ih2, kk, idom, icmp, ib, ii, ic, ic1
      integer :: idx

!===> ...  begin here

       spcodp = f_zero

!===> ... loop over vertical layers from top to surface

      lab_do_layer : do kk = 1, NLAY

! --- linear interp coeffs for rh-dep species

        ih2 = 1
        do while ( rh1(kk) > rhlev(ih2) )
          ih2 = ih2 + 1
          if ( ih2 > NRHLEV ) exit
        enddo
        ih1 = max( 1, ih2-1 )
        ih2 = min( NRHLEV, ih2 )

        drh0 = rhlev(ih2) - rhlev(ih1)
        drh1 = rh1(kk) - rhlev(ih1)
        if ( ih1 == ih2 ) then
          rdrh = f_zero
        else
          rdrh = drh1 / drh0
        endif

! --- assign optical properties in each domain

        idom = idmaer(kk)

        lab_if_idom : if (idom == 5) then
! --- 5th domain - upper stratosphere assume no aerosol

          do ib = 1, NSWLWBD
            tauae(kk,ib) = f_zero
            if ( ib <= NSWBND ) then
              ssaae(kk,ib) = 0.99
              asyae(kk,ib) = 0.696
            else
              ssaae(kk,ib) = 0.5
              asyae(kk,ib) = 0.3
            endif
          enddo

        elseif (idom == 4) then    lab_if_idom
! --- 4th domain - stratospheric layers

          do ib = 1, NSWLWBD
            tauae(kk,ib) = extstra(ib) * delz(kk)
            if ( ib <= NSWBND ) then
              ssaae(kk,ib) = 0.99
              asyae(kk,ib) = 0.696
            else
              ssaae(kk,ib) = 0.5
              asyae(kk,ib) = 0.3
            endif
          enddo

! --- compute aod from individual species' contribution (optional)
          idx = idxspc(10)             ! for sulfate
          spcodp(idx) = spcodp(idx) + tauae(kk,nv_aod)

        elseif (idom == 3) then    lab_if_idom
! --- 3rd domain - free tropospheric layers
!   1:inso 0.17e-3; 2:soot 0.4; 7:waso 0.59983; n:730

          do ib = 1, NSWLWBD
            ex01 = extrhi(1,ib)
            sc01 = scarhi(1,ib)
            ss01 = ssarhi(1,ib)
            as01 = asyrhi(1,ib)

            ex02 = extrhi(2,ib)
            sc02 = scarhi(2,ib)
            ss02 = ssarhi(2,ib)
            as02 = asyrhi(2,ib)

            ex03 = extrhd(ih1,1,ib)                                     &
     &           + rdrh * (extrhd(ih2,1,ib) - extrhd(ih1,1,ib))
            sc03 = scarhd(ih1,1,ib)                                     &
     &           + rdrh * (scarhd(ih2,1,ib) - scarhd(ih1,1,ib))
            ss03 = ssarhd(ih1,1,ib)                                     &
     &           + rdrh * (ssarhd(ih2,1,ib) - ssarhd(ih1,1,ib))
            as03 = asyrhd(ih1,1,ib)                                     &
     &           + rdrh * (asyrhd(ih2,1,ib) - asyrhd(ih1,1,ib))

            ext1 = 0.17e-3*ex01 + 0.4*ex02 + 0.59983*ex03
            sca1 = 0.17e-3*sc01 + 0.4*sc02 + 0.59983*sc03
            ssa1 = 0.17e-3*ss01*ex01 + 0.4*ss02*ex02 + 0.59983*ss03*ex03
            asy1 = 0.17e-3*as01*sc01 + 0.4*as02*sc02 + 0.59983*as03*sc03

            tauae(kk,ib) = ext1 * 730.0 * delz(kk)
            ssaae(kk,ib) = min(f_one, ssa1/ext1)
            asyae(kk,ib) = min(f_one, asy1/sca1)

! --- compute aod from individual species' contribution (optional)
            if ( ib==nv_aod ) then
             spcodp(1) = spcodp(1) + 0.17e-3*ex01*730.0*delz(kk)   ! dust (inso)   #1
             spcodp(2) = spcodp(2) + 0.4    *ex02*730.0*delz(kk)   ! black carbon  #2
             spcodp(3) = spcodp(3) + 0.59983*ex03*730.0*delz(kk)   ! water soluble #7
            endif

          enddo

        elseif (idom == 1) then    lab_if_idom
! --- 1st domain - mixing layer

          lab_do_ib : do ib = 1, NSWLWBD
            ext1 = f_zero
            sca1 = f_zero
            ssa1 = f_zero
            asy1 = f_zero

            lab_do_icmp : do icmp = 1, NCM
              ic = icmp
              idx = idxspc(icmp)

              cm = cmix(icmp)
              lab_if_cm : if ( cm > f_zero ) then

                lab_if_ic : if ( ic <= NCM1 ) then        ! component withour rh dep
                  tt0  = cm * extrhi(ic,ib)
                  ext1 = ext1 + tt0
                  sca1 = sca1 + cm * scarhi(ic,ib)
                  ssa1 = ssa1 + cm * ssarhi(ic,ib) * extrhi(ic,ib)
                  asy1 = asy1 + cm * asyrhi(ic,ib) * scarhi(ic,ib)
                else  lab_if_ic                           ! component with rh dep
                  ic1 = ic - NCM1

                  ex00 = extrhd(ih1,ic1,ib)                             &
     &               + rdrh * (extrhd(ih2,ic1,ib) - extrhd(ih1,ic1,ib))
                  sc00 = scarhd(ih1,ic1,ib)                             &
     &               + rdrh * (scarhd(ih2,ic1,ib) - scarhd(ih1,ic1,ib))
                  ss00 = ssarhd(ih1,ic1,ib)                             &
     &               + rdrh * (ssarhd(ih2,ic1,ib) - ssarhd(ih1,ic1,ib))
                  as00 = asyrhd(ih1,ic1,ib)                             &
     &               + rdrh * (asyrhd(ih2,ic1,ib) - asyrhd(ih1,ic1,ib))

                  tt0  = cm * ex00
                  ext1 = ext1 + tt0
                  sca1 = sca1 + cm * sc00
                  ssa1 = ssa1 + cm * ss00 * ex00
                  asy1 = asy1 + cm * as00 * sc00
                endif  lab_if_ic

! --- compute aod from individual species' contribution (optional)
                if ( ib==nv_aod ) then
                 spcodp(idx) = spcodp(idx) + tt0*denn(1)*delz(kk)   ! idx for dif species
                endif

              endif  lab_if_cm
            enddo  lab_do_icmp

            tauae(kk,ib) = ext1 * denn(1) * delz(kk)
            ssaae(kk,ib) = min(f_one, ssa1/ext1)
            asyae(kk,ib) = min(f_one, asy1/sca1)
          enddo  lab_do_ib

        elseif (idom == 2) then    lab_if_idom
! --- 2nd domain - mineral transport layers

          do ib = 1, NSWLWBD
            tauae(kk,ib) = extrhi(6,ib) * denn(2) * delz(kk)
            ssaae(kk,ib) = ssarhi(6,ib)
            asyae(kk,ib) = asyrhi(6,ib)
          enddo

! --- compute aod from individual species' contribution (optional)
          spcodp(1) = spcodp(1) + tauae(kk,nv_aod)            ! dust

        else  lab_if_idom
! --- domain index out off range, assume no aerosol

          do ib = 1, NSWLWBD
            tauae(kk,ib) = f_zero
            ssaae(kk,ib) = f_one
            asyae(kk,ib) = f_zero
          enddo

!         write(6,19) kk,idom
! 19      format(/'  ***  ERROR in sub AEROS: domain index out'         &
!    &,            ' of range!  K, IDOM =',3i5,' ***')
!         stop 19

        endif  lab_if_idom

      enddo  lab_do_layer

!
!===> ... smooth profile at domain boundaries
!
      if ( ivflip == 0 ) then    ! input from toa to sfc

        do ib = 1, NSWLWBD
        do kk = 2, NLAY
          if ( tauae(kk,ib) > f_zero ) then
            ratio = tauae(kk-1,ib) / tauae(kk,ib)
          else
            ratio = f_one
          endif

          tt0 = tauae(kk,ib) + tauae(kk-1,ib)
          tt1 = 0.2 * tt0
          tt2 = tt0 - tt1

          if ( ratio > crt1 ) then
            tauae(kk,ib)   = tt1
            tauae(kk-1,ib) = tt2
          endif

          if ( ratio < crt2 ) then
            tauae(kk,ib)   = tt2
            tauae(kk-1,ib) = tt1
          endif
        enddo   ! do_kk_loop
        enddo   ! do_ib_loop

      else                      ! input from sfc to toa

        do ib = 1, NSWLWBD
        do kk = NLAY-1, 1, -1
          if ( tauae(kk,ib) > f_zero ) then
            ratio = tauae(kk+1,ib) / tauae(kk,ib)
          else
            ratio = f_one
          endif

          tt0 = tauae(kk,ib) + tauae(kk+1,ib)
          tt1 = 0.2 * tt0
          tt2 = tt0 - tt1

          if ( ratio > crt1 ) then
            tauae(kk,ib)   = tt1
            tauae(kk+1,ib) = tt2
          endif

          if ( ratio < crt2 ) then
            tauae(kk,ib)   = tt2
            tauae(kk+1,ib) = tt1
          endif
        enddo   ! do_kk_loop
        enddo   ! do_ib_loop

      endif

!
      return
!................................
      end subroutine radclimaer
!--------------------------------
!
!...................................
      end subroutine aer_property
!-----------------------------------

!> @}
! =======================================================================
! GOCART code modification starts here (Sarah lu)  ---------------------!
!!
!! gocart_init : set_aerspc, rd_gocart_clim, rd_gocart_luts, optavg_grt
!! setgocartaer: aeropt_grt, map_aermr

!>\ingroup module_radiation_aerosols
!> The initialization program for gocart aerosols
!! - determine weight and index for aerosol composition/luts
!! - read in monthly global distribution of gocart aerosols
!! - read and map the tabulated aerosol optical spectral data onto
!!   corresponding SW/LW radiation spectral bands.
!!
!>\param NWVTOT           total num of wave numbers used in sw spectrum
!!\param solfwv           (NWVTOT), solar flux for each individual
!!                        wavenumber (w/m2)
!!\param soltot           total solar flux for the spectrual range (w/m2)
!!\param NWVTIR           total num of wave numbers used in the ir region
!!\param eirfwv           (NWVTIR), ir flux(273k) for each individual
!!                        wavenumber (w/m2)
!!\param NBDSW            num of bands calculated for sw aeros opt prop
!!\param NLWBND           num of bands calculated for lw aeros opt prop
!!\param NSWLWBD          total num of bands calc for sw+lw aeros opt prop
!!\param imon             month of the year
!!\param me               print message control flag
!!\param raddt            radiation time step
!!\param fdaer
!>\section gel_go_ini gocart_init General Algorithm
!! @{
!-----------------------------------
      subroutine gocart_init                                            &
     &     ( NWVTOT,solfwv,soltot,NWVTIR,eirfwv,                        &   !  ---  inputs:
     &       NBDSW,NLWBND,NSWLWBD,imon,me,raddt,fdaer                   &   !  ---  outputs: ( none )
     &     )

!  ==================================================================  !
!                                                                      !
!  subprogram : gocart_init                                            !
!                                                                      !
!    this is the initialization program for gocart aerosols            !
!                                                                      !
!    - determine weight and index for aerosol composition/luts         !
!    - read in monthly global distribution of gocart aerosols          !
!    - read and map the tabulated aerosol optical spectral data        !
!        onto corresponding sw/lw radiation spectral bands.            !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
!  inputs:                                                             !
!   NWVTOT           - total num of wave numbers used in sw spectrum   !
!   solfwv(NWVTOT)   - solar flux for each individual wavenumber (w/m2)!
!   soltot           - total solar flux for the spectrual range  (w/m2)!
!   NWVTIR           - total num of wave numbers used in the ir region !
!   eirfwv(NWVTIR)   - ir flux(273k) for each individual wavenum (w/m2)!
!   NBDSW            - num of bands calculated for sw aeros opt prop   !
!   NLWBND           - num of bands calculated for lw aeros opt prop   !
!   NSWLWBD          - total num of bands calc for sw+lw aeros opt prop!
!   imon             - month of the year                               !
!   me               - print message control flag                      !
!                                                                      !
!  outputs: (to the module variables)                                  !
!                                                                      !
!  module variables:                                                   !
!     NBDSW   - total number of sw spectral bands                      !
!     wvnum1,wvnum2 (NSWSTR:NSWEND)                                    !
!             - start/end wavenumbers for each of sw bands             !
!     NBDLW   - total number of lw spectral bands                      !
!     wvnlw1,wvnlw2 (NBDLW)                                            !
!             - start/end wavenumbers for each of lw bands             !
!     NSWLWBD - total number of sw+lw bands used in this version       !
!     extrhi_grt  - extinction coef for rh-indep aeros  KCM1*NSWLWBD   !
!     ssarhi_grt  - single-scat-alb for rh-indep aeros  KCM1*NSWLWBD   !
!     asyrhi_grt  - asymmetry factor for rh-indep aeros KCM1*NSWLWBD   !
!     extrhd_grt  - extinction coef for rh-dep aeros KRHLEV*KCM2*NSWLWBD!
!     ssarhd_grt  - single-scat-alb for rh-dep aeros KRHLEV*KCM2*NSWLWBD!
!     asyrhd_grt  - asymmetry factor for rh-dep aerosKRHLEV*KCM2*NSWLWBD!
!     ctaer       - merging coefficients for fcst/clim fields          !
!     get_fcst    - option to get fcst aerosol fields                  !
!     get_clim    - option to get clim aerosol fields                  !
!     dm_indx  - index for aer spec to be included in aeropt calculations  !
!     dmfcs_indx  - index for prognostic aerosol fields                !
!     psclmg      - geos3/4-gocart pressure            IMXG*JMXG*KMXG  !
!     dmclmg      - geos3-gocart aerosol dry mass   IMXG*JMXG*KMXG*NMXG!
!                   or geos4-gocart aerosol mixing ratio               !
!                                                                      !
!  usage:    call gocart_init                                          !
!                                                                      !
!  subprograms called:  set_aerspc, rd_gocart_clim,                    !
!                       rd_gocart_luts, optavg_grt                     !
!                                                                      !
!  ==================================================================  !

      implicit none

!  ---  inputs:
      integer, intent(in) :: NWVTOT,NWVTIR,NBDSW,NLWBND,NSWLWBD,imon,me

      real (kind=kind_phys), intent(in) :: raddt, fdaer

      real (kind=kind_phys), intent(in) :: solfwv(:),soltot, eirfwv(:)

!  ---  output: ( none )

!  ---  locals:

      real (kind=kind_phys), dimension(NBDSW,KAERBND)  :: solwaer
      real (kind=kind_phys), dimension(NBDSW)          :: solbnd
      real (kind=kind_phys), dimension(NLWBND,KAERBND) :: eirwaer
      real (kind=kind_phys), dimension(NLWBND)         :: eirbnd
      real (kind=kind_phys) :: sumsol, sumir

      integer, dimension(NBDSW)  :: nv1, nv2
      integer, dimension(NLWBND) :: nr1, nr2

      integer :: i, mb, ib, ii, iw, iw1, iw2

!===>  ...  begin here

!--------------------------------------------------------------------------
!  (1) determine aerosol specification index and merging coefficients
!--------------------------------------------------------------------------

      if ( .not. lgrtint ) then

!  --- ...  already done aerspc initialization, continue

        continue

      else

!  --- ...  set aerosol specification index and merging coefficients

        call set_aerspc(raddt,fdaer)
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

      endif  ! end if_lgrtinit_block

!
!--------------------------------------------------------------------------
!  (2) read gocart climatological data
!--------------------------------------------------------------------------

!  --- ...  read gocart climatological data, if needed

      if ( get_clim ) then

        call rd_gocart_clim
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

      endif

!
!--------------------------------------------------------------------------
!  (3) read and map the tabulated aerosol optical spectral data
!           onto corresponding radiation spectral bands
!--------------------------------------------------------------------------

      if ( .not. lgrtint ) then

!  --- ...  already done optical property interpolation, exit

        return

      else

!  --- ...  reset lgrtint

        lgrtint = .false.

!  --- ...  read tabulated aerosol optical input data
        call rd_gocart_luts
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

!  --- ...  compute solar flux weights and interval indices for mapping
!           spectral bands between sw radiation and aerosol data

        solbnd (:)   = f_zero
        solwaer(:,:) = f_zero

        nv_aod = 1

        do ib = 1, NBDSW
          mb = ib + NSWSTR - 1
          ii = 1
          iw1 = nint(wvnsw1(mb))
          iw2 = nint(wvnsw2(mb))
!
! ---  locate the spectral band for 550nm (for aod diag)
!
          if (10000./iw1 >= 0.55 .and.                                  &
     &        10000./iw2 <= 0.55 )  then
              nv_aod =  ib
          endif

          Lab_swdowhile : do while ( iw1 > iendwv_grt(ii) )
            if ( ii == KAERBND ) exit Lab_swdowhile
            ii = ii + 1
          enddo  Lab_swdowhile

          sumsol = f_zero
          nv1(ib) = ii

          do iw = iw1, iw2
            solbnd(ib) = solbnd(ib) + solfwv(iw)
            sumsol = sumsol + solfwv(iw)

            if ( iw == iendwv_grt(ii) ) then
              solwaer(ib,ii) = sumsol

              if ( ii < KAERBND ) then
                sumsol = f_zero
                ii = ii + 1
              endif
            endif
          enddo

          if ( iw2 /= iendwv_grt(ii) ) then
            solwaer(ib,ii) = sumsol
          endif

          nv2(ib) = ii

          if((me==0) .and. lckprnt) print *,'RAD-nv1,nv2:',             &
     &        ib,iw1,iw2,nv1(ib),iendwv_grt(nv1(ib)),                   &
     &        nv2(ib),iendwv_grt(nv2(ib)),                              &
     &        10000./iw1, 10000./iw2
        enddo     ! end do_ib_block for sw

! --- check the spectral range for the nv_550 band
        if((me==0) .and. lckprnt) then
          mb = nv_aod + NSWSTR - 1
          iw1 = nint(wvnsw1(mb))
          iw2 = nint(wvnsw2(mb))
           print *,'RAD-nv_aod:',                                       &
     &       nv_aod, iw1, iw2, 10000./iw1, 10000./iw2
        endif
!
!  --- ...  compute ir flux weights and interval indices for mapping
!           spectral bands between lw radiation and aerosol data

        eirbnd (:)   = f_zero
        eirwaer(:,:) = f_zero

        do ib = 1, NLWBND
          ii = 1
          if ( NLWBND == 1 ) then
            iw1 = 400                   ! corresponding 25 mu
            iw2 = 2500                  ! corresponding 4  mu
          else
            iw1 = nint(wvnlw1(ib))
            iw2 = nint(wvnlw2(ib))
          endif

          Lab_lwdowhile : do while ( iw1 > iendwv_grt(ii) )
            if ( ii == KAERBND ) exit Lab_lwdowhile
            ii = ii + 1
          enddo  Lab_lwdowhile

          sumir = f_zero
          nr1(ib) = ii

          do iw = iw1, iw2
            eirbnd(ib) = eirbnd(ib) + eirfwv(iw)
            sumir  = sumir  + eirfwv(iw)

            if ( iw == iendwv_grt(ii) ) then
              eirwaer(ib,ii) = sumir

              if ( ii < KAERBND ) then
                sumir = f_zero
                ii = ii + 1
              endif
            endif
          enddo

          if ( iw2 /= iendwv_grt(ii) ) then
            eirwaer(ib,ii) = sumir
          endif

          nr2(ib) = ii

          if(me==0 .and. lckprnt) print *,'RAD-nr1,nr2:',                &
     &        ib,iw1,iw2,nr1(ib),iendwv_grt(nr1(ib)),                    &
     &        nr2(ib),iendwv_grt(nr2(ib)),                               &
     &        10000./iw1, 10000./iw2
        enddo     ! end do_ib_block for lw

!  ---  compute spectral band mean properties for each species

        call optavg_grt
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

        if(me==0 .and. lckprnt) then
          print *, 'RAD -After optavg_grt, sw band info'
          do ib = 1, NBDSW
           mb = ib + NSWSTR - 1
           print *,'RAD -wvnsw1,wvnsw2: ',ib,wvnsw1(mb),wvnsw2(mb)
           print *,'RAD -lamda1,lamda2: ',ib,10000./wvnsw1(mb),         &
     &                                   10000./wvnsw2(mb)
           print *,'RAD -extrhi_grt:', extrhi_grt(:,ib)
!          do i = 1, KRHLEV
           do i = 1, KRHLEV, 10
             print *, 'RAD -extrhd_grt:',i,rhlev_grt(i),                &
     &                                extrhd_grt(i,:,ib)
           enddo
          enddo
          print *, 'RAD -After optavg_grt, lw band info'
          do ib = 1, NLWBND
           ii = NBDSW + ib
           print *,'RAD -wvnlw1,wvnlw2: ',ib,wvnlw1(ib),wvnlw2(ib)
           print *,'RAD -lamda1,lamda2: ',ib,10000./wvnlw1(ib),         &
     &                                   10000./wvnlw2(ib)
           print *,'RAD -extrhi_grt:', extrhi_grt(:,ii)
!          do i = 1, KRHLEV
           do i = 1, KRHLEV, 10
             print *, 'RAD -extrhd_grt:',i,rhlev_grt(i),                &
     &                                extrhd_grt(i,:,ii)
           enddo
          enddo
        endif

!  --- ...  dealoocate input data arrays no longer needed
        deallocate ( iendwv_grt   )
        if ( allocated(rhidext0_grt) ) then
          deallocate ( rhidext0_grt )
          deallocate ( rhidssa0_grt )
          deallocate ( rhidasy0_grt )
        endif
        if ( allocated(rhdpext0_grt) ) then
          deallocate ( rhdpext0_grt )
          deallocate ( rhdpssa0_grt )
          deallocate ( rhdpasy0_grt )
        endif

      endif  ! end if_lgrtinit_block

! =================
      contains
! =================

!>\ingroup module_radiation_aerosols
!> This subroutine determines merging coefficients ctaer; setup aerosol
!!  specification. The current version only supports prognostic aerosols
!! (from GOCART in-line calculations) and climo aerosols (from GEOS-GOCART
!! runs).
!!\section set_aerspc_gen set_aerspc General Algorithm
!! place holder
!-----------------------------
      subroutine set_aerspc(raddt,fdaer)
!.............................
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

! ==================================================================== !
!                                                                      !
! subprogram: set_aerspc                                               !
!                                                                      !
! determine merging coefficients ctaer;                                !
! set up aerosol specification: num_gridcomp, gridcomp, dm_indx,       !
!                       dmfcs_indx, isoot, iwaso, isuso, issam, isscm  !
!                                                                      !
! Aerosol optical properties (ext, ssa, asy) are determined from       !
! NMGX (<=12) aerosol species                                          !
! ==> DU: dust1 (4 sub-micron bins), dust2, dust3, dust4, dust5        !
!     BC: soot_phobic, soot_philic                                     !
!     OC: waso_phobic, waso_philic                                     !
!     SU: suso (=so4)                                                  !
!     SS: ssam (accumulation mode), sscm (coarse mode)                 !
!                                                                      !
! The current version only supports prognostic aerosols (from GOCART   !
! in-line calculations) and climo aerosols (from GEOS-GOCART runs)    !
!                                                                      !
!  ==================================================================  !
!
      implicit none

!  ---  inputs:
      real (kind=kind_phys), intent(in) :: raddt, fdaer
!  ---  output:

!  ---  local:
!     real (kind=kind_phys)     :: raddt
      integer                   :: i, indxr
      character*2               :: tp, gridcomp_tmp(max_num_gridcomp)

!! ===> determine ctaer (user specified weight for fcst fields)
!     raddt = min(fhswr,fhlwr) / 24.
      if( fdaer >= 99999. ) ctaer = f_one
      if((fdaer>0.).and.(fdaer<99999.)) ctaer=exp(-raddt/fdaer)

      if(me==0 .and. lckprnt) then
        print *, 'RAD -raddt, fdaer,ctaer: ', raddt, fdaer, ctaer
        if (ctaer == f_one ) then
          print *, 'LU -aerosol fields determined from fcst'
        elseif (ctaer == f_zero) then
          print *, 'LU -aerosol fields determined from clim'
        else
          print *, 'LU -aerosol fields determined from fcst/clim'
        endif
      endif

!! ===> determine get_fcst and get_clim
!!    if fcst is chosen (ctaer == f_one ), set get_clim to F
!!    if clim is chosen (ctaer == f_zero), set get_fcst to F
      if ( ctaer == f_one  )  get_clim = .false.
      if ( ctaer == f_zero )  get_fcst = .false.

!! ===> determine aerosol species to be included in the calculations
!!      of aerosol optical properties (ext, ssa, asy)

!*  If climo option is chosen, the aerosol composition is hardwired
!*  to full package. If not, the composition is determined from
!*  tracer_config on-the-fly (full package or subset)
      lab_if_fcst : if ( get_fcst ) then

!!      use tracer_config to determine num_gridcomp and gridcomp
        if ( gfs_phy_tracer%doing_GOCART )  then
         if ( gfs_phy_tracer%doing_DU )  then
            num_gridcomp  =  num_gridcomp  + 1
            gridcomp_tmp(num_gridcomp) = 'DU'
         endif
         if ( gfs_phy_tracer%doing_SU ) then
            num_gridcomp  =  num_gridcomp  + 1
            gridcomp_tmp(num_gridcomp) = 'SU'
         endif
         if ( gfs_phy_tracer%doing_SS ) then
            num_gridcomp  =  num_gridcomp  + 1
            gridcomp_tmp(num_gridcomp) = 'SS'
         endif
         if ( gfs_phy_tracer%doing_OC ) then
            num_gridcomp  =  num_gridcomp  + 1
            gridcomp_tmp(num_gridcomp) = 'OC'
         endif
         if ( gfs_phy_tracer%doing_BC ) then
            num_gridcomp  =  num_gridcomp  + 1
            gridcomp_tmp(num_gridcomp) = 'BC'
         endif
!
         if ( num_gridcomp > 0 ) then
           allocate ( gridcomp(num_gridcomp) )
           gridcomp(1:num_gridcomp) = gridcomp_tmp(1:num_gridcomp)
         else
           print *,'ERROR: prognostic aerosols not found,abort',me
           stop 1000
         endif

        else      ! gfs_phy_tracer%doing_GOCART=F

         print *,'ERROR: prognostic aerosols option off, abort',me
        stop 1001

        endif     ! end_if_gfs_phy_tracer%doing_GOCART_if_

      else lab_if_fcst

!!      set to full package (max_num_gridcomp and max_gridcomp)
        num_gridcomp = max_num_gridcomp
        allocate ( gridcomp(num_gridcomp) )
        gridcomp(1:num_gridcomp) = max_gridcomp(1:num_gridcomp)

      endif lab_if_fcst

!!
!! Aerosol specification is determined as such:
!! A. For radiation-aerosol feedback, the specification is based on the aeropt
!!    routine from Mian Chin and Hongbin Yu (hydrophobic and hydrophilic for
!!    OC/BC; submicron and supermicron for SS, 8-bins (with 4 subgroups for the
!!    the submicron bin) for DU, and SO4 for SU)
!! B. For transport, the specification is determined from GOCART in-line module
!! C. For LUTS, (waso, soot, ssam, sscm, suso, dust) is used, based on the
!!    the OPAC climo aerosol scheme (implemented by Yu-Tai Hou)

!!=== <A> determine dm_indx and NMXG
      indxr = 0
      dm_indx%waso_phobic = -999         ! OC
      dm_indx%soot_phobic = -999         ! BC
      dm_indx%ssam = -999                ! SS
      dm_indx%suso = -999                ! SU
      dm_indx%dust1 = -999               ! DU
      do i = 1, num_gridcomp
         tp = gridcomp(i)
         select case ( tp )
         case ( 'OC')    ! consider hydrophobic and hydrophilic
           dm_indx%waso_phobic = indxr + 1
           dm_indx%waso_philic = indxr + 2
           indxr = indxr + 2
         case ( 'BC')    ! consider hydrophobic and hydrophilic
           dm_indx%soot_phobic = indxr + 1
           dm_indx%soot_philic = indxr + 2
           indxr = indxr + 2
         case ( 'SS')    ! consider submicron and supermicron
           dm_indx%ssam = indxr + 1
           dm_indx%sscm = indxr + 2
           indxr = indxr + 2
         case ( 'SU')    ! consider SO4 only
           dm_indx%suso = indxr + 1
           indxr = indxr + 1
         case ( 'DU')    ! consider all 5 bins
           dm_indx%dust1 = indxr + 1
           dm_indx%dust2 = indxr + 2
           dm_indx%dust3 = indxr + 3
           dm_indx%dust4 = indxr + 4
           dm_indx%dust5 = indxr + 5
           indxr = indxr + 5
         case default
           print *,'ERROR: aerosol species not supported, abort',me
           stop 1002
         end select
      enddo
!!
      NMXG       = indxr      ! num of gocart aer spec for opt cal
!!

!!=== <B> determine dmfcs_indx
!!    SS: 5-bins are considered for transport while only two groups
!!        (accumulation/coarse modes) are considered for radiation
!!    DU: 5-bins are considered for transport while 8 bins (with the
!!        submicorn bin exptended to 4 bins) are considered for radiation
!!    SU: DMS, SO2, and MSA are not considered for radiation

      if ( get_fcst ) then
         if ( gfs_phy_tracer%doing_OC )  then
            dmfcs_indx%ocphobic = trcindx ('ocphobic', gfs_phy_tracer)
            dmfcs_indx%ocphilic = trcindx ('ocphilic', gfs_phy_tracer)
         endif
         if ( gfs_phy_tracer%doing_BC )  then
            dmfcs_indx%bcphobic = trcindx ('bcphobic', gfs_phy_tracer)
            dmfcs_indx%bcphilic = trcindx ('bcphilic', gfs_phy_tracer)
         endif
         if ( gfs_phy_tracer%doing_SS )  then
            dmfcs_indx%ss001 = trcindx ('ss001', gfs_phy_tracer)
            dmfcs_indx%ss002 = trcindx ('ss002', gfs_phy_tracer)
            dmfcs_indx%ss003 = trcindx ('ss003', gfs_phy_tracer)
            dmfcs_indx%ss004 = trcindx ('ss004', gfs_phy_tracer)
            dmfcs_indx%ss005 = trcindx ('ss005', gfs_phy_tracer)
         endif
         if ( gfs_phy_tracer%doing_SU )  then
            dmfcs_indx%so4 = trcindx ('so4', gfs_phy_tracer)
         endif
         if ( gfs_phy_tracer%doing_DU )  then
            dmfcs_indx%du001 = trcindx ('du001', gfs_phy_tracer)
            dmfcs_indx%du002 = trcindx ('du002', gfs_phy_tracer)
            dmfcs_indx%du003 = trcindx ('du003', gfs_phy_tracer)
            dmfcs_indx%du004 = trcindx ('du004', gfs_phy_tracer)
            dmfcs_indx%du005 = trcindx ('du005', gfs_phy_tracer)
         endif
      endif

!!
!!=== <C> determin KCM, KCM1, KCM2
!!    DU: submicron bin (dust1) contains 4 sub-groups (e.g., hardwire
!!        8 bins for aerosol optical properties luts)
!!    OC/BC: while hydrophobic aerosols are rh-independent, the luts
!!        for hydrophilic aerosols are used (e.g., use the coeff
!!        corresponding to rh=0)
!!
      indxr = 1
      isoot = -999
      iwaso = -999
      isuso = -999
      issam = -999
      isscm = -999
      do i = 1, num_gridcomp
         tp = gridcomp(i)
         if ( tp /= 'DU' ) then  !<--- non-dust aerosols
           select case ( tp )
           case ( 'OC ')
             iwaso = indxr
           case ( 'BC ')
             isoot = indxr
           case ( 'SU ')
             isuso = indxr
           case ( 'SS ')
             issam = indxr
             isscm = indxr + 1
           end select
           if ( tp /= 'SS' ) then
             indxr = indxr + 1
           else
             indxr = indxr + 2
           endif
         else                   !<--- dust aerosols
           KCM1 =  8            ! num of rh independent aer species
         endif
      enddo
      KCM2 = indxr - 1          ! num of rh dependent aer species
      KCM  = KCM1 + KCM2        ! total num of aer species

!!
!! check print starts here
      if( me == 0 .and. lckprnt) then
       print *, 'RAD -num_gridcomp:', num_gridcomp
       print *, 'RAD -gridcomp    :', gridcomp(:)
       print *, 'RAD -NMXG:',  NMXG
       print *, 'RAD -dm_indx ===> '
       print *, 'RAD -aerspc: dust1=', dm_indx%dust1
       print *, 'RAD -aerspc: dust2=', dm_indx%dust2
       print *, 'RAD -aerspc: dust3=', dm_indx%dust3
       print *, 'RAD -aerspc: dust4=', dm_indx%dust4
       print *, 'RAD -aerspc: dust5=', dm_indx%dust5
       print *, 'RAD -aerspc: ssam=',  dm_indx%ssam
       print *, 'RAD -aerspc: sscm=',  dm_indx%sscm
       print *, 'RAD -aerspc: suso=',  dm_indx%suso
       print *, 'RAD -aerspc: waso_phobic=',dm_indx%waso_phobic
       print *, 'RAD -aerspc: waso_philic=',dm_indx%waso_philic
       print *, 'RAD -aerspc: soot_phobic=',dm_indx%soot_phobic
       print *, 'RAD -aerspc: soot_philic=',dm_indx%soot_philic

       print *, 'RAD -KCM1 =', KCM1
       print *, 'RAD -KCM2 =', KCM2
       print *, 'RAD -KCM  =', KCM
       if ( KCM2 > 0 ) then
         print *, 'RAD -aerspc: issam=', issam
         print *, 'RAD -aerspc: isscm=', isscm
         print *, 'RAD -aerspc: isuso=', isuso
         print *, 'RAD -aerspc: iwaso=', iwaso
         print *, 'RAD -aerspc: isoot=', isoot
       endif

       if ( get_fcst ) then
         print *, 'RAD -dmfcs_indx ===> '
         print *, 'RAD -trc_du001=',dmfcs_indx%du001
         print *, 'RAD -trc_du002=',dmfcs_indx%du002
         print *, 'RAD -trc_du003=',dmfcs_indx%du003
         print *, 'RAD -trc_du004=',dmfcs_indx%du004
         print *, 'RAD -trc_du005=',dmfcs_indx%du005
         print *, 'RAD -trc_so4  =',dmfcs_indx%so4
         print *, 'RAD -trc_ocphobic=',dmfcs_indx%ocphobic
         print *, 'RAD -trc_ocphilic=',dmfcs_indx%ocphilic
         print *, 'RAD -trc_bcphobic=',dmfcs_indx%bcphobic
         print *, 'RAD -trc_bcphilic=',dmfcs_indx%bcphilic
         print *, 'RAD -trc_ss001=',dmfcs_indx%ss001
         print *, 'RAD -trc_ss002=',dmfcs_indx%ss002
         print *, 'RAD -trc_ss003=',dmfcs_indx%ss003
         print *, 'RAD -trc_ss004=',dmfcs_indx%ss004
         print *, 'RAD -trc_ss005=',dmfcs_indx%ss005
       endif
      endif
!! check print ends here

      return
!                                                                      !
      end subroutine set_aerspc

!-----------------------------------
!>\ingroup module_radiation_aerosols
!> This subroutine reads input gocart aerosol optical data from Mie
!! code calculations.
!>\section rd_gocart_luts_gen rd_gocart_luts General Algorithm
!-----------------------------
      subroutine rd_gocart_luts
!.............................
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

! ==================================================================== !
! subprogram: rd_gocart_luts                                           !
!   read input gocart aerosol optical data from Mie code calculations  !
!                                                                      !
! Remarks (Quanhua (Mark) Liu, JCSDA, June 2008)                       !
!  The LUT is for NCEP selected 61 wave numbers and 6 aerosols         !
!  (dust, soot, suso, waso, ssam, and sscm) and 36 aerosol effective   !
!  size in microns.                                                    !
!                                                                      !
!  The LUT is computed using Mie code with a logorithm size            !
!  distribution for each of 36 effective sizes. The standard deviation !
!  sigma of the size, and min/max size follows Chin et al. 2000        !
!  For each effective size, it corresponds a relative humidity value.  !
!                                                                      !
!  The LUT contains the density, sigma, relative humidity, mean mode   !
!  radius, effective size, mass extinction coefficient, single         !
!  scattering albedo, asymmetry factor, and phase function             !
!                                                                      !
!  ==================================================================  !
!
      implicit none

!  ---  inputs:
!  ---  output:

!  ---  locals:
      INTEGER, PARAMETER :: NP = 100, NP2 = 2*NP, nWave=100,            &
     &                      nAero=6, n_p=36
      INTEGER :: NW, NS, nH, n_bin
      real (kind=kind_io8), Dimension( NP2 ) :: Angle, Cos_Angle,       &
     &                                          Cos_Weight
      real (kind=kind_io8), Dimension(n_p,nAero) :: RH, rm, reff
      real (kind=kind_io8), Dimension(nWave,n_p,nAero) ::               &
     &                      ext0, sca0, asy0
      real (kind=kind_io8), Dimension(NP2,n_p,nWave,nAero) :: ph0
      real (kind=kind_io8) :: wavelength(nWave), density(nAero),        &
     &                        sigma(nAero), wave,n_fac,PI,t1,s1,g1
      CHARACTER(len=80) :: AerosolName(nAero)
      INTEGER    :: i, j, k, l, ij

      character  :: aerosol_file*30
      logical    :: file_exist
      integer    :: indx_dust(8)          ! map 36 dust bins to gocart size bins

      data aerosol_file  /"NCEP_AEROSOL.bin"/
      data AerosolName/ ' Dust ', ' Soot ', ' SUSO ', ' WASO ',         &
     &                  ' SSAM ', ' SSCM '/

!! 8 dust bins
!!  1       2       3       4       5       6       7       8
!! .1-.18, .18-.3, .3-.6, 0.6-1.0, 1.0-1.8, 1.8-3, 3-6,  6-10  <-- def
!!  0.1399  0.2399  0.4499 0.8000 1.3994  2.3964 4.4964  7.9887 <-- reff
      data indx_dust/4, 8, 12, 18, 21, 24, 30, 36/

      PI = acos(-1.d0)

! -- allocate aerosol optical data
      if ( .not. allocated( iendwv_grt ) ) then
        allocate ( iendwv_grt (KAERBND) )
      endif
      if (.not. allocated(rhidext0_grt) .and. KCM1 > 0 ) then
        allocate ( rhidext0_grt(KAERBND,KCM1))
        allocate ( rhidssa0_grt(KAERBND,KCM1))
        allocate ( rhidasy0_grt(KAERBND,KCM1))
      endif
      if (.not. allocated(rhdpext0_grt) .and. KCM2 > 0 ) then
        allocate ( rhdpext0_grt(KAERBND,KRHLEV,KCM2))
        allocate ( rhdpssa0_grt(KAERBND,KRHLEV,KCM2))
        allocate ( rhdpasy0_grt(KAERBND,KRHLEV,KCM2))
      endif

! -- read luts
      inquire (file = aerosol_file, exist = file_exist)

      if ( file_exist ) then
        if(me==0 .and. lckprnt) print *,'RAD -open :',aerosol_file
        close (NIAERCM)
        open (unit=NIAERCM,file=aerosol_file,status='OLD',              &
     &        action='read',form='UNFORMATTED')
      else
        print *,'    Requested aerosol data file "',aerosol_file,       &
     &          '" not found!', me
        print *,'    *** Stopped in subroutine RD_GOCART_LUTS !!'
        stop 1003
      endif              ! end if_file_exist_block

      READ(NIAERCM) (Cos_Angle(i),i=1,NP)
      READ(NIAERCM) (Cos_Weight(i),i=1,NP)
      READ(NIAERCM)
      READ(NIAERCM)
      READ(NIAERCM) NW,NS
      READ(NIAERCM)
      READ(NIAERCM) (wavelength(i),i=1,NW)

! --- check nAero and NW
      if (NW /= KAERBND) then
        print *, "Incorrect spectral band, abort ", NW
        stop 1004
      endif

! --- convert wavelength to wavenumber
      do i = 1, KAERBND
       iendwv_grt(i) = 10000. / wavelength(i)
       if(me==0 .and. lckprnt) print *,'RAD -wn,lamda:',                &
     &           i,iendwv_grt(i),wavelength(i)
      enddo

      DO j = 1, nAero
        if(me==0 .and. lckprnt) print *,'RAD -read LUTs:',              &
     &            j,AerosolName(j)
        READ(NIAERCM)
        READ(NIAERCM)
        READ(NIAERCM) n_bin, density(j), sigma(j)
        READ(NIAERCM)
        READ(NIAERCM) (RH(i,j),i=1, n_bin)
        READ(NIAERCM)
        READ(NIAERCM) (rm(i,j),i=1, n_bin)
        READ(NIAERCM)
        READ(NIAERCM) (reff(i,j),i=1, n_bin)

! --- check n_bin
        if (n_bin /= KRHLEV ) then
          print *, "Incorrect rh levels, abort ", n_bin
          stop 1005
        endif

! --- read luts
        DO k = 1, NW
          READ(NIAERCM) wave,(ext0(k,L,j),L=1,n_bin)
          READ(NIAERCM) (sca0(k,L,j),L=1,n_bin)
          READ(NIAERCM) (asy0(k,L,j),L=1,n_bin)
          READ(NIAERCM) (ph0(1:NP2,L,k,j),L=1,n_bin)
        END DO

! --- map luts input to module variables
        if (AerosolName(j) == ' Dust ' ) then
         if ( KCM1 > 0) then    !<-- only if rh independent aerosols are needed
          do i = 1, KCM1
           rhidext0_grt(1:KAERBND,i)=ext0(1:KAERBND,indx_dust(i),j)
           rhidssa0_grt(1:KAERBND,i)=sca0(1:KAERBND,indx_dust(i),j)
           rhidasy0_grt(1:KAERBND,i)=asy0(1:KAERBND,indx_dust(i),j)
          enddo
         endif
        else
         if ( KCM2 > 0) then    !<-- only if rh dependent aerosols are needed
          if (AerosolName(j) == ' Soot ') ij = isoot
          if (AerosolName(j) == ' SUSO ') ij = isuso
          if (AerosolName(j) == ' WASO ') ij = iwaso
          if (AerosolName(j) == ' SSAM ') ij = issam
          if (AerosolName(j) == ' SSCM ') ij = isscm
          if ( ij .ne. -999 ) then
             rhdpext0_grt(1:KAERBND,1:KRHLEV,ij) =                      &
     &               ext0(1:KAERBND,1:KRHLEV,j)
             rhdpssa0_grt(1:KAERBND,1:KRHLEV,ij) =                      &
     &               sca0(1:KAERBND,1:KRHLEV,j)
             rhdpasy0_grt(1:KAERBND,1:KRHLEV,ij) =                      &
     &               asy0(1:KAERBND,1:KRHLEV,j)
          endif   ! if_ij
         endif    ! if_KCM2
        endif
      END DO

      return
!...................................
      end subroutine rd_gocart_luts
!-----------------------------------
!                                                                      !
!>\ingroup module_radiation_aerosols
!> This subroutine computes mean aerosols optical properties over each
!! SW/LW radiation spectral band for each of the species components.
!! This program follows GFDL's approach for thick cloud optical property
!! in SW radiation scheme (2000).
!>\section optavg_grt_gen optavg_grt General Algorithm
!! @{
!-----------------------------
      subroutine optavg_grt
!.............................
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

! ==================================================================== !
!                                                                      !
! subprogram: optavg_grt                                               !
!                                                                      !
!   compute mean aerosols optical properties over each sw/lw radiation !
!   spectral band for each of the species components.  This program    !
!   follows gfdl's approach for thick cloud opertical property in      !
!   sw radiation scheme (2000).                                        !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
! input arguments:                                                     !
!   nv1,nv2 (NBDSW)  - start/end spectral band indices of aerosol data !
!                      for each sw radiation spectral band             !
!   nr1,nr2 (NLWBND)  - start/end spectral band indices of aerosol data !
!                      for each ir radiation spectral band             !
!   solwaer (NBDSW,KAERBND)                                            !
!                    - solar flux weight over each sw radiation band   !
!                      vs each aerosol data spectral band              !
!   eirwaer (NLWBND,KAERBND)                                            !
!                    - ir flux weight over each lw radiation band      !
!                      vs each aerosol data spectral band              !
!   solbnd  (NBDSW)  - solar flux weight over each sw radiation band   !
!   eirbnd  (NLWBND) - ir flux weight over each lw radiation band      !
!   NBDSW            - total number of sw spectral bands               !
!   NLWBND           - total number of lw spectral bands               !
!   NSWLWBD          - total number of sw+lw spectral bands            !
!                                                                      !
! output arguments: (to module variables)                              !
!                                                                      !
!  ==================================================================  !
!
      implicit none

!  ---  inputs:
!  ---  output:

!  ---  locals:
      real (kind=kind_phys) :: sumk, sumok, sumokg, sumreft,            &
     &       sp, refb, reft, rsolbd, rirbd

      integer :: ib, nb, ni, nh, nc
!
!===> ...  begin here

!  --- ...  allocate aerosol optical data
      if (.not. allocated(extrhd_grt) .and. KCM2 > 0 ) then
        allocate ( extrhd_grt(KRHLEV,KCM2,NSWLWBD) )
        allocate ( ssarhd_grt(KRHLEV,KCM2,NSWLWBD) )
        allocate ( asyrhd_grt(KRHLEV,KCM2,NSWLWBD) )
      endif
      if (.not. allocated(extrhi_grt) .and. KCM1 > 0 ) then
        allocate ( extrhi_grt(KCM1,NSWLWBD) )
        allocate ( ssarhi_grt(KCM1,NSWLWBD) )
        allocate ( asyrhi_grt(KCM1,NSWLWBD) )
      endif
!
!  --- ...  loop for each sw radiation spectral band

      do nb = 1, NBDSW
        rsolbd = f_one / solbnd(nb)

!  ---  for rh independent aerosol species

        lab_rhi: if (KCM1 >  0 ) then
        do nc = 1, KCM1
          sumk    = f_zero
          sumok   = f_zero
          sumokg  = f_zero
          sumreft = f_zero

          do ni = nv1(nb), nv2(nb)
            sp   = sqrt( (f_one - rhidssa0_grt(ni,nc))                  &
     &           / (f_one - rhidssa0_grt(ni,nc)*rhidasy0_grt(ni,nc)) )
            reft = (f_one - sp) / (f_one + sp)
            sumreft = sumreft + reft*solwaer(nb,ni)

            sumk    = sumk    + rhidext0_grt(ni,nc)*solwaer(nb,ni)
            sumok   = sumok   + rhidssa0_grt(ni,nc)*solwaer(nb,ni)      &
     &              * rhidext0_grt(ni,nc)
            sumokg  = sumokg  + rhidssa0_grt(ni,nc)*solwaer(nb,ni)      &
     &              * rhidext0_grt(ni,nc)*rhidasy0_grt(ni,nc)
          enddo

          refb = sumreft * rsolbd

          extrhi_grt(nc,nb) = sumk   * rsolbd
          asyrhi_grt(nc,nb) = sumokg / (sumok + 1.0e-10)
          ssarhi_grt(nc,nb) = 4.0*refb                                  &
     &      / ( (f_one+refb)**2 - asyrhi_grt(nc,nb)*(f_one-refb)**2 )

        enddo   ! end do_nc_block for rh-ind aeros
        endif lab_rhi

!  ---  for rh dependent aerosols species

        lab_rhd: if (KCM2 > 0 ) then
        do nc = 1, KCM2
          do nh = 1, KRHLEV
            sumk    = f_zero
            sumok   = f_zero
            sumokg  = f_zero
            sumreft = f_zero

            do ni = nv1(nb), nv2(nb)
              sp   = sqrt( (f_one - rhdpssa0_grt(ni,nh,nc))             &
     &        /(f_one-rhdpssa0_grt(ni,nh,nc)*rhdpasy0_grt(ni,nh,nc)))
              reft = (f_one - sp) / (f_one + sp)
              sumreft = sumreft + reft*solwaer(nb,ni)

              sumk    = sumk   + rhdpext0_grt(ni,nh,nc)*solwaer(nb,ni)
              sumok   = sumok  + rhdpssa0_grt(ni,nh,nc)*solwaer(nb,ni)  &
     &                * rhdpext0_grt(ni,nh,nc)
              sumokg  = sumokg + rhdpssa0_grt(ni,nh,nc)*solwaer(nb,ni)  &
     &                * rhdpext0_grt(ni,nh,nc)*rhdpasy0_grt(ni,nh,nc)
            enddo

            refb = sumreft * rsolbd

            extrhd_grt(nh,nc,nb) = sumk   * rsolbd
            asyrhd_grt(nh,nc,nb) = sumokg / (sumok + 1.0e-10)
            ssarhd_grt(nh,nc,nb) = 4.0*refb                             &
     &      /((f_one+refb)**2 - asyrhd_grt(nh,nc,nb)*(f_one-refb)**2)
          enddo   ! end do_nh_block
        enddo   ! end do_nc_block for rh-dep aeros
        endif lab_rhd

      enddo   !  end do_nb_block for sw

!  --- ...  loop for each lw radiation spectral band

      do nb = 1, NLWBND

        ib = NBDSW + nb
        rirbd = f_one / eirbnd(nb)

!  ---  for rh independent aerosol species

        lab_rhi_lw: if (KCM1 > 0 ) then
        do nc = 1, KCM1
          sumk    = f_zero
          sumok   = f_zero
          sumokg  = f_zero
          sumreft = f_zero

          do ni = nr1(nb), nr2(nb)
            sp   = sqrt( (f_one - rhidssa0_grt(ni,nc))                  &
     &      / (f_one - rhidssa0_grt(ni,nc)*rhidasy0_grt(ni,nc)) )
            reft = (f_one - sp) / (f_one + sp)
            sumreft = sumreft + reft*eirwaer(nb,ni)

            sumk    = sumk    + rhidext0_grt(ni,nc)*eirwaer(nb,ni)
            sumok   = sumok   + rhidssa0_grt(ni,nc)*eirwaer(nb,ni)      &
     &              * rhidext0_grt(ni,nc)
            sumokg  = sumokg  + rhidssa0_grt(ni,nc)*eirwaer(nb,ni)      &
     &              * rhidext0_grt(ni,nc)*rhidasy0_grt(ni,nc)
          enddo

          refb = sumreft * rirbd

          extrhi_grt(nc,ib) = sumk   * rirbd
          asyrhi_grt(nc,ib) = sumokg / (sumok + 1.0e-10)
          ssarhi_grt(nc,ib) = 4.0*refb                                  &
     &    / ( (f_one+refb)**2 - asyrhi_grt(nc,ib)*(f_one-refb)**2 )
        enddo   ! end do_nc_block for rh-ind aeros
        endif lab_rhi_lw

!  ---  for rh dependent aerosols species

        lab_rhd_lw: if (KCM2 > 0 ) then
        do nc = 1, KCM2
          do nh = 1, KRHLEV
            sumk    = f_zero
            sumok   = f_zero
            sumokg  = f_zero
            sumreft = f_zero

            do ni = nr1(nb), nr2(nb)
              sp   = sqrt( (f_one - rhdpssa0_grt(ni,nh,nc))             &
     &        /(f_one - rhdpssa0_grt(ni,nh,nc)*rhdpasy0_grt(ni,nh,nc)) )
              reft = (f_one - sp) / (f_one + sp)
              sumreft = sumreft + reft*eirwaer(nb,ni)

              sumk    = sumk  + rhdpext0_grt(ni,nh,nc)*eirwaer(nb,ni)
              sumok   = sumok + rhdpssa0_grt(ni,nh,nc)*eirwaer(nb,ni)   &
     &                * rhdpext0_grt(ni,nh,nc)
              sumokg  = sumokg+ rhdpssa0_grt(ni,nh,nc)*eirwaer(nb,ni)   &
     &                * rhdpext0_grt(ni,nh,nc)*rhdpasy0_grt(ni,nh,nc)
            enddo

            refb = sumreft * rirbd

            extrhd_grt(nh,nc,ib) = sumk   * rirbd
            asyrhd_grt(nh,nc,ib) = sumokg / (sumok + 1.0e-10)
            ssarhd_grt(nh,nc,ib) = 4.0*refb                             &
     &      /((f_one+refb)**2 - asyrhd_grt(nh,nc,ib)*(f_one-refb)**2 )
          enddo   ! end do_nh_block
        enddo   ! end do_nc_block for rh-dep aeros
        endif lab_rhd_lw

      enddo   !  end do_nb_block for lw

!
      return
!................................
      end subroutine optavg_grt
!! @}
!--------------------------------
!
!>\ingroup module_radiation_aerosols
!> This subroutine:
!! - 1. read in aerosol dry mass and surface pressure from GEOS3-GOCART
!! C3.1 2000 monthly dataset or aerosol mixing ratio and surface
!! pressure from GEOS4-GOCART 2000-2007 averaged monthly data set.
!! - 2. compute goes lat/lon array (for horizontal mapping)
!>\section rd_gocart_clim_gen rd_gocart_clim General Algorithm
!! @{
!-----------------------------------
      subroutine rd_gocart_clim
!...................................
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

!  ==================================================================  !
!                                                                      !
! subprogram: rd_gocart_clim                                           !
!                                                                      !
!   1. read in aerosol dry mass and surface pressure from GEOS3-GOCART !
!      C3.1 2000 monthly data set                                      !
!      or aerosol mixing ratio and surface pressure from GEOS4-GOCART  !
!      2000-2007 averaged monthly data set                             !
!   2. compute goes lat/lon array (for horizontal mapping)             !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
! inputs arguments:                                                    !
!     imon    - month of the year                                      !
!     me      - print message control flag                             !
!                                                                      !
! outputs arguments: (to the module variables)                         !
!     psclmg   - pressure (sfc to toa)    cb   IMXG*JMXG*KMXG          !
!     dmclmg   - aerosol dry mass/mixing ratio     IMXG*JMXG*KMXG*NMXG !
!     geos_rlon - goes longitude          deg      IMXG                !
!     geos_rlat - goes latitude           deg      JMXG                !
!                                                                      !
!  usage:    call rd_gocart_clim                                       !
!                                                                      !
!  program history:                                                    !
!    05/18/2010  ---  Lu    Add the option to read GEOS4-GOCART climo  !
!  ==================================================================  !
!
      implicit none

!  ---  inputs:
!  ---  output:

!  ---  locals:
      integer, parameter :: MAXSPC = 5
      real (kind=kind_io4), parameter  :: PINT = 0.01
      real (kind=kind_io4), parameter  :: EPSQ = 0.0

      integer         :: i, j, k, numspci, ii
      integer         :: icmp, nrecl, nt1, nt2, nn(MAXSPC)
      character       :: ymd*6, yr*4, mn*2, tp*2,                       &
     &                   fname*30, fin*30, aerosol_file*40
      logical         :: file_exist

      real (kind=kind_io4), dimension(KMXG)             :: sig
      real (kind=kind_io4), dimension(IMXG,JMXG)        :: ps
      real (kind=kind_io4), dimension(IMXG,JMXG,KMXG)   :: temp
      real (kind=kind_io4), dimension(IMXG,JMXG,KMXG,MAXSPC):: buff
      real (kind=kind_phys)   :: pstmp

!     Add the following variables for GEOS4-GOCART
      real (kind=kind_io4), dimension(KMXG):: hyam, hybm
      real (kind=kind_io4)                 :: p0

      data yr /'2000'/           !!<=== use 2000 as the climo proxy

!* sigma_coordinate for GEOS3-GOCART
!* P(i,j,k) = PINT + SIG(k) * (PS(i,j) - PINT)
      data SIG  /                                                       &
     &     9.98547E-01,9.94147E-01,9.86350E-01,9.74300E-01,9.56950E-01, &
     &     9.33150E-01,9.01750E-01,8.61500E-01,8.11000E-01,7.50600E-01, &
     &     6.82900E-01,6.10850E-01,5.37050E-01,4.63900E-01,3.93650E-01, &
     &     3.28275E-01,2.69500E-01,2.18295E-01,1.74820E-01,1.38840E-01, &
     &     1.09790E-01,8.66900E-02,6.84150E-02,5.39800E-02,4.25750E-02, &
     &     3.35700E-02,2.39900E-02,1.36775E-02,5.01750E-03,5.30000E-04 /

!* hybrid_sigma_pressure_coordinate for GEOS4-GOCART
!* p(i,j,k) = a(k)*p0 + b(k)*ps(i,j)
      data hyam/                                                        &
     &   0, 0.0062694, 0.02377049, 0.05011813, 0.08278809, 0.1186361,   &
     &   0.1540329, 0.1836373, 0.2043698, 0.2167788, 0.221193,          &
     &   0.217729, 0.2062951, 0.1865887, 0.1615213, 0.1372958,          &
     &   0.1167039, 0.09920014, 0.08432171, 0.06656809, 0.04765031,     &
     &   0.03382346, 0.0237648, 0.01435208, 0.00659734, 0.002826232,    &
     &   0.001118959, 0.0004086494, 0.0001368611, 3.750308e-05/

      data hybm /                                                       &
     &   0.992555, 0.9642, 0.90556, 0.816375, 0.703815, 0.576585,       &
     &   0.44445, 0.324385, 0.226815, 0.149165, 0.089375,               &
     &   0.045865, 0.017485, 0.00348, 0, 0, 0, 0, 0,                    &
     &   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /

      data p0 /1013.25 /

!===>  ...  begin here

! --- allocate and initialize gocart climatological data
      if ( .not. allocated (dmclmg) ) then
        allocate ( dmclmg(IMXG,JMXG,KMXG,NMXG) )
        allocate ( psclmg(IMXG,JMXG,KMXG) )
        allocate ( molwgt(NMXG) )
      endif

      dmclmg(:,:,:,:)  = f_zero
      psclmg(:,:,:)    = f_zero
      molwgt(:)        = f_zero

! --- allocate and initialize geos lat and lon arrays
      if ( .not. allocated ( geos_rlon  )) then
          allocate (geos_rlon(IMXG))
          allocate (geos_rlat(JMXG))
      endif

      geos_rlon(:) = f_zero
      geos_rlat(:) = f_zero

! --- compute geos lat and lon arrays
      do i = 1, IMXG
         geos_rlon(i)     = -180. + (i-1)* dltx
      end do
      do j = 2, JMXG-1
         geos_rlat(j)     = -90. + (j-1)* dlty
      end do
      geos_rlat(1)      = -89.5
      geos_rlat(JMXG)   =  89.5

! --- determine whether GEOS3 or GEOS4 data set is provided
      if ( gocart_climo == 'xxxx' ) then
        gocart_climo='0000'
! check geos3-gocart climo
        aerosol_file = '200001.PS.avg'
        inquire (file = aerosol_file, exist = file_exist)
        if ( file_exist ) gocart_climo='ver3'
! check geos4-gocart climo
        aerosol_file = 'gocart_climo_2000x2007_ps_01.bin'
        inquire (file = aerosol_file, exist = file_exist)
        if ( file_exist ) gocart_climo='ver4'
      endif
!
!
! --- read ps (sfc pressure) and compute 3d pressure field (psclmg)
!
      write(mn,'(i2.2)') imon
      ymd = yr//mn
      aerosol_file = 'null'
      if ( gocart_climo == 'ver3' ) then
        aerosol_file = ymd//'.PS.avg'
      elseif ( gocart_climo == 'ver4' ) then
        aerosol_file = 'gocart_climo_2000x2007_ps_'//mn//'.bin'
      endif
!
      inquire (file = aerosol_file, exist = file_exist)
      lab_if_ps : if ( file_exist ) then

       close(NIAERCM)
       if ( gocart_climo == 'ver3' ) then
        nrecl = 4 * (IMXG * JMXG)
        open(NIAERCM, file=trim(aerosol_file),                          &
     &       action='read',access='direct',recl=nrecl)
        read(NIAERCM, rec=1) ps
        do j = 1, JMXG
          do i = 1, IMXG
            do k = 1, KMXG
              pstmp = pint + sig(k) * (ps(i,j) - pint)
              psclmg(i,j,k) = 0.1 * pstmp       ! convert mb to cb
            enddo
          enddo
        enddo

       elseif ( gocart_climo == 'ver4' ) then
         open(NIAERCM, file=trim(aerosol_file),                         &
     &        action='read',status='old', form='unformatted')
         read(NIAERCM) ps(:,:)
         do j = 1, JMXG
           do i = 1, IMXG
             do k = 1, KMXG
              pstmp = hyam(k)*p0 + hybm(k)*ps(i,j)
              psclmg(i,j,k) = 0.1 * pstmp       ! convert mb to cb
            enddo
          enddo
         enddo

       endif     !  ---- end if_gocart_climo

      else lab_if_ps

        print *,' *** Requested aerosol data file "',                    &
     &          trim(aerosol_file),  '" not found!'
        print *,' *** Stopped in RD_GOCART_CLIM ! ', me
        stop 1006
      endif   lab_if_ps
!
! --- read aerosol dry mass (g/m3) or mixing ratios (mol/mol,kg/kg)
!
      lab_do_icmp : do icmp = 1, num_gridcomp

         tp = gridcomp(icmp)

!        determine aerosol_file
         aerosol_file = 'null'
         if ( gocart_climo == 'ver3' ) then
           if(tp == 'DU')   fname='.DU.STD.tv20.g.avg'
           if(tp == 'SS')   fname='.SS.STD.tv17.g.avg'
           if(tp == 'SU')   fname='.SU.STD.tv15.g.avg'
           if(tp == 'OC')   fname='.CC.STD.tv15.g.avg'
           if(tp == 'BC')   fname='.CC.STD.tv15.g.avg'
           aerosol_file=ymd//trim(fname)
         elseif ( gocart_climo == 'ver4' ) then
           fin = 'gocart_climo_2000x2007_'
           if(tp == 'DU')   fname=trim(fin)//'du_'
           if(tp == 'SS')   fname=trim(fin)//'ss_'
           if(tp == 'SU')   fname=trim(fin)//'su_'
           if(tp == 'OC')   fname=trim(fin)//'cc_'
           if(tp == 'BC')   fname=trim(fin)//'cc_'
           aerosol_file=trim(fname)//mn//'.bin'
         endif

         numspci = 4
         if(tp == 'DU')   numspci = 5
         inquire (file=trim(aerosol_file), exist = file_exist)
         lab_if_aer: if ( file_exist ) then
!
          close(NIAERCM)
          if ( gocart_climo == 'ver3' ) then
          nrecl = 4 * numspci * (IMXG * JMXG * KMXG + 3)
          open (NIAERCM, file=trim(aerosol_file),                       &
     &         action='read',access='direct', recl=nrecl)
          read(NIAERCM,rec=1)(nt1,nt2,nn(i),buff(:,:,:,i),i=1,numspci)

          elseif ( gocart_climo == 'ver4' ) then
           open (NIAERCM, file=trim(aerosol_file),                      &
     &           action='read',status='old', form='unformatted')
           do i = 1, numspci
             do k = 1, KMXG
              read(NIAERCM) temp(:,:,k)
              buff(:,:,k,i) =  temp(:,:,k)
             enddo
           enddo
          endif

!!===> fill dmclmg with working array buff
          select case ( tp )

! fill in DU from DU: du1, du2, du3, du4, du5
          case ('DU' )
           if ( dm_indx%dust1 /= -999) then
            do ii = 1, 5
             dmclmg(:,:,:,dm_indx%dust1+ii-1) = buff(:,:,:,ii)
            enddo
           else
            print *, 'ERROR: invalid DU index, abort! ',me
            stop 1007
           endif

! fill in BC from CC: bc_phobic, oc_phobic, bc_philic, oc_philic
          case ('BC' )
           if ( dm_indx%soot_phobic /= -999) then
            dmclmg(:,:,:,dm_indx%soot_phobic)=buff(:,:,:,1)
            dmclmg(:,:,:,dm_indx%soot_philic)=buff(:,:,:,3)
            molwgt(dm_indx%soot_phobic) = 12.
            molwgt(dm_indx%soot_philic) = 12.
           else
            print *, 'ERROR: invalid BC index, abort! ',me
            stop 1008
           endif

! fill in SU from SU: dms, so2, so4, msa
          case ('SU' )
           if ( dm_indx%suso /= -999) then
            dmclmg(:,:,:,dm_indx%suso) = buff(:,:,:,3)
            molwgt(dm_indx%suso) = 96.
           else
            print *, 'ERROR: invalid SU index, abort! ',me
            stop 1009
           endif

! fill in OC from CC: bc_phobic, oc_phobic, bc_philic, oc_philic
          case ('OC' )
           if ( dm_indx%waso_phobic /= -999) then
            dmclmg(:,:,:,dm_indx%waso_phobic) = 1.4*buff(:,:,:,2)
            dmclmg(:,:,:,dm_indx%waso_philic) = 1.4*buff(:,:,:,4)
            molwgt(dm_indx%waso_phobic) = 12.
            molwgt(dm_indx%waso_philic) = 12.
           else
            print *, 'ERROR: invalid OC index, abort! ',me
            stop 1010
           endif

! fill in SS from SS: ss1, ss2, ss3, ss4
          case ('SS' )
           if ( dm_indx%ssam /= -999) then
            dmclmg(:,:,:,dm_indx%ssam) = buff(:,:,:,1)
            dmclmg(:,:,:,dm_indx%sscm) = buff(:,:,:,2) +                &
     &                     buff(:,:,:,3)+buff(:,:,:,4)
           else
            print *, 'ERROR: invalid SS index, abort! ',me
            stop  1011
           endif

          case default

            print *, 'ERROR: invalid aerosol species, abort ',tp
            stop  1012

          end select

         else   lab_if_aer
          print *,' *** Requested aerosol data file "',aerosol_file,    &
     &            '" not found!'
          print *,' *** Stopped in RD_GOCART_CLIM ! ', me
          stop 1013
         endif  lab_if_aer

       enddo lab_do_icmp

      return
!...................................
      end subroutine rd_gocart_clim
!! @}
!-----------------------------------
!
!...................................
      end subroutine gocart_init
!-----------------------------------
!! @}

!>\ingroup module_radiation_aerosols
!> This subroutine computes SW + LW aerosol optical properties for
!! gocart aerosol species (merged from fcst and clim fields).
!!
!>\param alon      IMAX, longitude of given points in degree
!!\param alat      IMAX, latitude of given points in degree
!!\param prslk     (IMAX,NLAY), pressure in cb
!!\param rhlay     (IMAX,NLAY), layer mean relative humidity
!!\param dz        (IMAX,NLAY), layer thickness in m
!!\param hz        (IMAX,NLP1), level high in m
!!\param NSWLWBD   total number of sw+ir bands for aeros opt prop
!!\param prsl      (IMAX,NLAY), layer mean pressure in mb
!!\param tvly      (IMAX,NLAY), layer mean virtual temperature in K
!!\param trcly     (IMAX,NLAY,NTRAC), layer mean specific tracer in g/g
!!\param IMAX      horizontal dimension of arrays
!!\param NLAY,NLP1 vertical dimensions of arrays
!!\param ivflip    control flag for direction of vertical index
!!\n               =0: index from toa to surface
!!\n               =1: index from surface to toa
!!\param lsswr,lslwr       logical flag for sw/lw radiation calls
!!\param aerosw    (IMAX,NLAY,NBDSW,NF_AESW), aeros opt properties for SW
!!\n               (:,:,:,1): optical depth
!!\n               (:,:,:,2): single scattering albedo
!!\n               (:,:,:,3): asymmetry parameter
!!\param aerolw    (IMAX,NLAY,NBDLW,NF_AELW), aeros opt properties for LW
!!\n               (:,:,:,1): optical depth
!!\n               (:,:,:,2): single scattering albedo
!!\n               (:,:,:,3): asymmetry parameter
!>\section gen_setgo setgocartaer General Algorithm
!!@{
!-----------------------------------
      subroutine setgocartaer                                           &
     &     ( alon,alat,prslk,rhlay,dz,hz,NSWLWBD,                       &       !  ---  inputs:
     &       prsl,tvly,trcly,                                           &
     &       IMAX,NLAY,NLP1, ivflip, lsswr,lslwr,                       &
     &       aerosw,aerolw                                              &       !  ---  outputs:
     &     )


!  ==================================================================  !
!                                                                      !
!  setgocartaer computes sw + lw aerosol optical properties for gocart !
!  aerosol species (merged from fcst and clim fields)                  !
!                                                                      !
!  inputs:                                                             !
!     alon, alat                                             IMAX      !
!             - longitude and latitude of given points in degree       !
!     prslk   - pressure                           cb   IMAX*NLAY      !
!     rhlay   - layer mean relative humidity            IMAX*NLAY      !
!     dz      - layer thickness                    m    IMAX*NLAY      !
!     hz      - level high                         m    IMAX*NLP1      !
!     NSWLWBD - total number of sw+ir bands for aeros opt prop  1      !
!     prsl    - layer mean pressure                mb   IMAX*NLAY      !
!     tvly    - layer mean virtual temperature     k   IMAX*NLAY      !
!     trcly   - layer mean specific tracer         g/g  IMAX*NLAY*NTRAC!
!     IMAX    - horizontal dimension of arrays                  1      !
!     NLAY,NLP1-vertical dimensions of arrays                   1      !
!     ivflip   - control flag for direction of vertical index   1      !
!               =0: index from toa to surface                          !
!               =1: index from surface to toa                          !
!     lsswr,lslwr                                                      !
!             - logical flag for sw/lw radiation calls          1      !
!                                                                      !
!  outputs:                                                            !
!     aerosw - aeros opt properties for sw      IMAX*NLAY*NBDSW*NF_AESW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!     aerolw - aeros opt properties for lw      IMAX*NLAY*NBDLW*NF_AELW!
!               (:,:,:,1): optical depth                               !
!               (:,:,:,2): single scattering albedo                    !
!               (:,:,:,3): asymmetry parameter                         !
!     tau_gocart - 550nm aeros opt depth     IMAX*NLAY*MAX_NUM_GRIDCOMP!
!                                                                      !
!  module parameters and constants:                                    !
!     NBDSW   - total number of sw bands for aeros opt prop     1      !
!     NLWBND  - total number of ir bands for aeros opt prop     1      !
!                                                                      !
!  module variable: (set by subroutine gocart_init)                    !
!     dmclmg  - aerosols dry mass/mixing ratios   IMXG*JMXG*KMXG*NMXG  !
!     psclmg  - pressure                cb        IMXG*JMXG*KMXG       !
!                                                                      !
!  usage:    call setgocartaer                                         !
!                                                                      !
!  subprograms called:  map_aermr, aeropt_grt                          !
!                                                                      !
!  ==================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: IMAX,NLAY,NLP1,ivflip,NSWLWBD
      logical, intent(in) :: lsswr, lslwr

      real (kind=kind_phys), dimension(:,:),   intent(in) :: prslk,     &
     &       prsl, rhlay, tvly, dz, hz
      real (kind=kind_phys), dimension(:),     intent(in) :: alon, alat
      real (kind=kind_phys), dimension(:,:,:), intent(in) :: trcly

!  ---  outputs:
      real (kind=kind_phys), dimension(:,:,:,:), intent(out) ::         &
     &       aerosw, aerolw

!  ---  locals:
      real (kind=kind_phys), dimension(NLAY) :: rh1, dz1
      real (kind=kind_phys), dimension(NLAY,NSWLWBD)::tauae,ssaae,asyae
      real (kind=kind_phys), dimension(NLAY,max_num_gridcomp) ::        &
     &                       tauae_gocart

      real (kind=kind_phys) :: tmp1, tmp2

      integer               :: i, i1, i2, j1, j2, k, m, m1, kp

!     prognostic aerosols on gfs grids
      real (kind=kind_phys), dimension(:,:,:),allocatable:: aermr,dmfcs

! aerosol (dry mass) on gfs grids/levels
      real (kind=kind_phys), dimension(:,:), allocatable ::             &
     &  dmanl,dmclm, dmclmx
      real (kind=kind_phys), dimension(KMXG)     :: pstmp, pkstr
      real (kind=kind_phys) :: ptop, psfc, tem, plv, tv, rho

!  ---  conversion constants
      real (kind=kind_phys), parameter :: hdltx = 0.5 * dltx
      real (kind=kind_phys), parameter :: hdlty = 0.5 * dlty

!===>  ...  begin here
!
      if ( .not. allocated(dmanl) ) then
        allocate ( dmclmx(KMXG,NMXG) )
        allocate ( dmanl(NLAY,NMXG) )
        allocate ( dmclm(NLAY,NMXG) )

        allocate ( aermr(IMAX,NLAY,NMXG) )
        allocate ( dmfcs(IMAX,NLAY,NMXG) )
      endif
!
!> -# Call map_aermr() to map input tracer array (trcly) to local
!!    tracer array (aermr).
      dmfcs(:,:,:) = f_zero
      lab_if_fcst : if ( get_fcst ) then

        call map_aermr
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

      endif       lab_if_fcst
!
!> -# Map geos-gocart climo (dmclmg) to gfs grids (dmclm).
      lab_do_IMAX : do i = 1, IMAX

        dmclm(:,:) = f_zero

        lab_if_clim : if ( get_clim ) then
!  ---  map grid in longitude direction
          i2 = 1
          j2 = 1
          tmp1 = alon(i)
          if (tmp1 > 180.) tmp1 = tmp1 - 360.0
          lab_do_IMXG : do i1 = 1, IMXG
            tmp2 = geos_rlon(i1)
            if (tmp2 > 180.) tmp2 = tmp2 - 360.0
            if (abs(tmp1-tmp2) <= hdltx) then
              i2 = i1
              exit lab_do_IMXG
            endif
          enddo  lab_do_IMXG

!  ---  map grid in latitude direction
          lab_do_JMXG : do j1 = 1, JMXG
            if (abs(alat(i)-geos_rlat(j1)) <= hdlty) then
              j2 = j1
              exit lab_do_JMXG
            endif
          enddo  lab_do_JMXG

!  ---  update local arrays pstmp and dmclmx
          pstmp(:)= psclmg(i2,j2,:)*1000.0      ! cb to Pa
          dmclmx(:,:) = dmclmg(i2,j2,:,:)

!  ---  map geos-gocart climo (dmclmx) to gfs level (dmclm)
          pkstr(:)=fpkap(pstmp(:))
          psfc = pkstr(1)                       ! pressure at sfc
          ptop = pkstr(KMXG)                    ! pressure at toa

!  ---  map grid in verical direction (follow how ozone is mapped
!       in radiation_gases routine)
          do k = 1, NLAY
           kp = k                              ! from sfc to toa
           if(ivflip==0) kp = NLAY - k + 1     ! from toa to sfc
           tmp1 = prslk(i,kp)

           do m1 = 1, KMXG - 1                 ! from sfc to toa
             if(tmp1 > pkstr(m1+1) .and. tmp1 <= pkstr(m1)) then
               tmp2 = f_one / (pkstr(m1)-pkstr(m1+1))
               tem = (pkstr(m1) - tmp1) * tmp2
               dmclm(kp,:) = tem * dmclmx(m1+1,:)+                      &
     &                   (f_one-tem) * dmclmx(m1,:)
             endif
           enddo

!*         if(tmp1 > psfc) dmclm(kp,:) = dmclmx(1,:)
!*         if(tmp1 < ptop) dmclm(kp,:) = dmclmx(KMXG,:)

          enddo
        endif    lab_if_clim
!
!  ---  compute fcst/clim merged aerosol loading (dmanl) and the
!       radiation optical properties (aerosw, aerolw)
!
        do k = 1, NLAY

!  ---  map global to local arrays (rh1 and dz1)
          rh1(k) = rhlay(i,k)
          dz1(k) = dz   (i,k)

!  ---  convert from mixing ratio to dry mass (g/m3)
          plv = 100. * prsl(i,k)       ! convert pressure from mb to Pa
          tv  = tvly(i,k)              ! virtual temp in K
          rho = plv / (con_rd * tv)    ! air density in kg/m3
          if ( get_fcst ) then
            do m = 1,  NMXG            ! mixing ratio (g/g)
            dmfcs(i,k,m) = max(1000.*(rho*aermr(i,k,m)),f_zero)
            enddo     ! m_do_loop
          endif
          if ( get_clim .and. (gocart_climo == 'ver4') ) then
            do m = 1,  NMXG
             dmclm(k,m)=1000.*dmclm(k,m)*rho !mixing ratio (g/g)
             if ( molwgt(m) /= 0. ) then     !mixing ratio (mol/mol)
              dmclm(k,m)=dmclm(k,m) * (molwgt(m)/con_amd)
             endif
            enddo     ! m_do_loop
          endif


!  ---  determine dmanl from dmclm and dmfcs
          do m = 1, NMXG
             dmanl(k,m)= ctaer*dmfcs(i,k,m) +                           &
     &                 ( f_one-ctaer)*dmclm(k,m)
          enddo
        enddo

!> -# Call aeropt_grt() to alculate sw/lw aerosol optical properties
!!    for the corresponding frequency bands.

        call aeropt_grt
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

        if ( lsswr ) then

          if ( laswflg ) then

            do m = 1, NBDSW
              do k = 1, NLAY
                aerosw(i,k,m,1) = tauae(k,m)
                aerosw(i,k,m,2) = ssaae(k,m)
                aerosw(i,k,m,3) = asyae(k,m)
              enddo
            enddo

          else

            aerosw(:,:,:,:) = f_zero

          endif

        endif     ! end if_lsswr_block

        if ( lslwr ) then

          if ( lalwflg ) then

            if ( NLWBND == 1 ) then
              m1 = NBDSW + 1
              do m = 1, NBDLW
                do k = 1, NLAY
                  aerolw(i,k,m,1) = tauae(k,m1)
                  aerolw(i,k,m,2) = ssaae(k,m1)
                  aerolw(i,k,m,3) = asyae(k,m1)
                enddo
              enddo
            else
              do m = 1, NBDLW
                m1 = NBDSW + m
                do k = 1, NLAY
                  aerolw(i,k,m,1) = tauae(k,m1)
                  aerolw(i,k,m,2) = ssaae(k,m1)
                  aerolw(i,k,m,3) = asyae(k,m1)
                enddo
              enddo
            endif

          else

            aerolw(:,:,:,:) = f_zero

          endif
        endif     ! end if_lslwr_block

      enddo  lab_do_IMAX

! =================
      contains
! =================

!>\ingroup module_radiation_aerosols
!> This subroutine maps input tracer fields (trcly) to local tracer
!! array (aermr).
!>\section map_aermr_gen map_aermr General Algorithm
!! @{
!-----------------------------
      subroutine map_aermr
!.............................
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

! ==================================================================== !
!                                                                      !
! subprogram: map_aermr                                                !
!                                                                      !
!   map input tracer fields (trcly) to local tracer array (aermr)      !
!                                                                      !
!  ====================  defination of variables  ===================  !
!                                                                      !
! input arguments:                                                     !
!     IMAX    - horizontal dimension of arrays                  1      !
!     NLAY    - vertical dimensions of arrays                   1      !
!     trcly   - layer tracer mass mixing ratio     g/g  IMAX*NLAY*NTRAC!
! output arguments: (to module variables)                              !
!     aermr   - layer aerosol mass mixing ratio    g/g  IMAX*NLAY*NMXG !
!                                                                      !
! note:                                                                !
!  NTRAC is the number of tracers excluding water vapor                !
!  NMXG is the number of prognostic aerosol species                    !
!  ==================================================================  !
!
      implicit none

!  ---  inputs:
!  ---  output:

!  ---  local:
      integer    :: i, indx, ii
      character  :: tp*2

! initialize
      aermr(:,:,:) = f_zero
      ii = 1        !! <---- trcly does not contain q

! ==>  DU: du1 (submicron bins), du2, du3, du4, du5
       if( gfs_phy_tracer%doing_DU ) then
         aermr(:,:,dm_indx%dust1) = trcly(:,:,dmfcs_indx%du001-ii)
         aermr(:,:,dm_indx%dust2) = trcly(:,:,dmfcs_indx%du002-ii)
         aermr(:,:,dm_indx%dust3) = trcly(:,:,dmfcs_indx%du003-ii)
         aermr(:,:,dm_indx%dust4) = trcly(:,:,dmfcs_indx%du004-ii)
         aermr(:,:,dm_indx%dust5) = trcly(:,:,dmfcs_indx%du005-ii)
       endif

! ==>  OC: oc_phobic, oc_philic
       if( gfs_phy_tracer%doing_OC ) then
         aermr(:,:,dm_indx%waso_phobic) =                               &
     &                     trcly(:,:,dmfcs_indx%ocphobic-ii)
         aermr(:,:,dm_indx%waso_philic) =                               &
     &                     trcly(:,:,dmfcs_indx%ocphilic-ii)
       endif

! ==>  BC: bc_phobic, bc_philic
       if( gfs_phy_tracer%doing_BC ) then
         aermr(:,:,dm_indx%soot_phobic) =                               &
     &                     trcly(:,:,dmfcs_indx%bcphobic-ii)
         aermr(:,:,dm_indx%soot_philic) =                               &
     &                     trcly(:,:,dmfcs_indx%bcphilic-ii)
       endif

! ==>  SS: ss1, ss2 (submicron bins), ss3, ss4, ss5
       if( gfs_phy_tracer%doing_SS ) then
          aermr(:,:,dm_indx%ssam) = trcly(:,:,dmfcs_indx%ss001-ii)      &
     &                            + trcly(:,:,dmfcs_indx%ss002-ii)
          aermr(:,:,dm_indx%sscm) = trcly(:,:,dmfcs_indx%ss003-ii)      &
     &                            + trcly(:,:,dmfcs_indx%ss004-ii)      &
     &                            + trcly(:,:,dmfcs_indx%ss005-ii)
       endif

! ==>  SU: so4
       if( gfs_phy_tracer%doing_SU ) then
          aermr(:,:,dm_indx%suso) = trcly(:,:,dmfcs_indx%so4-ii)
       endif

      return
!...................................
      end subroutine map_aermr
!! @}
!-----------------------------------


!>\ingroup module_radiation_aerosols
!! This subroutine computes aerosols optical properties in NSWLWBD
!! SW/LW bands. Aerosol distribution at each grid point is composed
!! from up to NMXG aerosol species (from NUM_GRIDCOMP components).
!>\section aeropt_grt_gen aeropt_grt General Algorithm
!! @{ 
!-----------------------------------
      subroutine aeropt_grt
!...................................
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

!  ==================================================================  !
!                                                                      !
!  subprogram: aeropt_grt                                              !
!                                                                      !
!  compute aerosols optical properties in NSWLWBD sw/lw bands.         !
!  Aerosol distribution at each grid point is composed from up to      !
!  NMXG aerosol species (from NUM_GRIDCOMP components).                !
!                                                                      !
!  input variables:                                                    !
!     dmanl  - aerosol dry mass                     g/m3   NLAY*NMXG   !
!     rh1    - relative humidity                     %     NLAY        !
!     dz1    - layer thickness                       km    NLAY        !
!     NLAY   - vertical dimensions                   -     1           !
!     ivflip - control flag for direction of vertical index            !
!               =0: index from toa to surface                          !
!               =1: index from surface to toa                          !
!                                                                      !
!  output variables:                                                   !
!     tauae  - aerosol optical depth                 -   NLAY*NSWLWBD  !
!     ssaae  - aerosol single scattering albedo      -   NLAY*NSWLWBD  !
!     asyae  - aerosol asymmetry parameter           -   NLAY*NSWLWBD  !
!                                                                      !
!  ==================================================================  !
!
      implicit none

!  ---  inputs:
!  ---  outputs:

!  ---  locals:
      real (kind=kind_phys) :: aerdm
      real (kind=kind_phys) :: ext1, ssa1, asy1, ex00, ss00, as00,      &
     &                         ex01, ss01, as01, exint
      real (kind=kind_phys) :: tau, ssa, asy,                           &
     &                         sum_tau, sum_ssa, sum_asy

!  ---  subgroups for sub-micron dust
!  ---  corresponds to 0.1-0.18, 0.18-0.3, 0.3-0.6, 0.6-1.0 micron

      real (kind=kind_phys) ::   fd(4)
      data fd  / 0.01053,0.08421,0.25263,0.65263 /

      character    :: tp*2
      integer      :: icmp, n, kk, ib, ih2, ih1, ii, ij, ijk
      real (kind=kind_phys) :: drh0, drh1, rdrh

      real (kind=kind_phys) :: qmin  !<--lower bound for opt calc
      data qmin  / 1.e-20 /

!===>  ...  begin here

! --- initialize (assume no aerosols)
      tauae = f_zero
      ssaae = f_one
      asyae = f_zero

      tauae_gocart = f_zero

!===> ... loop over vertical layers
!
      lab_do_layer : do kk = 1, NLAY

! --- linear interp coeffs for rh-dep species

        ih2 = 1
        do while ( rh1(kk) > rhlev_grt(ih2) )
          ih2 = ih2 + 1
          if ( ih2 > KRHLEV ) exit
        enddo
        ih1 = max( 1, ih2-1 )
        ih2 = min( KRHLEV, ih2 )

        drh0 = rhlev_grt(ih2) - rhlev_grt(ih1)
        drh1 = rh1(kk) - rhlev_grt(ih1)
        if ( ih1 == ih2 ) then
          rdrh = f_zero
        else
          rdrh = drh1 / drh0
        endif

! --- loop through sw/lw spectral bands

        lab_do_ib : do ib = 1, NSWLWBD
          sum_tau = f_zero
          sum_ssa = f_zero
          sum_asy = f_zero

! --- loop through aerosol grid components
          lab_do_icmp : do icmp = 1, NUM_GRIDCOMP
            ext1 = f_zero
            ssa1 = f_zero
            asy1 = f_zero

            tp = gridcomp(icmp)

            select case ( tp )

! -- dust aerosols: no humidification effect
            case ( 'DU')
              do n = 1, KCM1

                if (n <= 4) then
                  aerdm = dmanl(kk,dm_indx%dust1) * fd(n)
                else
                  aerdm = dmanl(kk,dm_indx%dust1+n-4 )
                endif

                if (aerdm < qmin) aerdm = f_zero
                ex00 = extrhi_grt(n,ib)*(1000.*dz1(kk))*aerdm
                ss00 = ssarhi_grt(n,ib)
                as00 = asyrhi_grt(n,ib)
                ext1 = ext1 + ex00
                ssa1 = ssa1 + ex00 * ss00
                asy1 = asy1 + ex00 * ss00 * as00

              enddo

! -- suso aerosols: with humidification effect
            case ( 'SU')
              ij = isuso
              exint = extrhd_grt(ih1,ij,ib)                             &
     &          + rdrh*(extrhd_grt(ih2,ij,ib) - extrhd_grt(ih1,ij,ib))
              ss00 = ssarhd_grt(ih1,ij,ib)                              &
     &          + rdrh*(ssarhd_grt(ih2,ij,ib) - ssarhd_grt(ih1,ij,ib))
              as00 = asyrhd_grt(ih1,ij,ib)                              &
     &          + rdrh*(asyrhd_grt(ih2,ij,ib) - asyrhd_grt(ih1,ij,ib))

              aerdm = dmanl(kk, dm_indx%suso)
              if (aerdm < qmin) aerdm = f_zero
              ex00 = exint*(1000.*dz1(kk))*aerdm
              ext1 = ex00
              ssa1 = ex00 * ss00
              asy1 = ex00 * ss00 * as00

! -- seasalt aerosols: with humidification effect
            case ( 'SS')
              do n = 1, 2                  !<---- ssam, sscm
                ij = issam + (n-1)
                exint = extrhd_grt(ih1,ij,ib)                           &
     &          + rdrh*(extrhd_grt(ih2,ij,ib) - extrhd_grt(ih1,ij,ib))
                ss00 = ssarhd_grt(ih1,ij,ib)                            &
     &          + rdrh*(ssarhd_grt(ih2,ij,ib) - ssarhd_grt(ih1,ij,ib))
                as00 = asyrhd_grt(ih1,ij,ib)                            &
     &          + rdrh*(asyrhd_grt(ih2,ij,ib) - asyrhd_grt(ih1,ij,ib))

                aerdm = dmanl(kk, dm_indx%ssam+n-1)
                if (aerdm < qmin) aerdm = f_zero
                ex00 = exint*(1000.*dz1(kk))*aerdm
                ext1 = ext1 + ex00
                ssa1 = ssa1 + ex00 * ss00
                asy1 = asy1 + ex00 * ss00 * as00

              enddo

! -- organic carbon/black carbon:
!    using 'waso' and 'soot' for hydrophilic OC and BC
!    using 'waso' and 'soot' at RH=0 for hydrophobic OC and BC
            case ( 'OC', 'BC')
              if(tp == 'OC') then
                 ii = dm_indx%waso_phobic
                 ij = iwaso
              else
                 ii = dm_indx%soot_phobic
                 ij = isoot
              endif

! --- hydrophobic
              aerdm = dmanl(kk, ii)
              if (aerdm < qmin) aerdm = f_zero
              ex00 = extrhd_grt(1,ij,ib)*(1000.*dz1(kk))*aerdm
              ss00 = ssarhd_grt(1,ij,ib)
              as00 = asyrhd_grt(1,ij,ib)
! --- hydrophilic
              aerdm = dmanl(kk, ii+1)
              if (aerdm < qmin) aerdm = f_zero
              exint = extrhd_grt(ih1,ij,ib)                             &
     &         + rdrh*(extrhd_grt(ih2,ij,ib) - extrhd_grt(ih1,ij,ib))
              ex01 = exint*(1000.*dz1(kk))*aerdm
              ss01 = ssarhd_grt(ih1,ij,ib)                              &
     &          + rdrh*(ssarhd_grt(ih2,ij,ib) - ssarhd_grt(ih1,ij,ib))
              as01 = asyrhd_grt(ih1,ij,ib)                              &
     &          + rdrh*(asyrhd_grt(ih2,ij,ib) - asyrhd_grt(ih1,ij,ib))

              ext1 = ex00 + ex01
              ssa1 = (ex00 * ss00) + (ex01 * ss01)
              asy1 = (ex00 * ss00 * as00) + (ex01 * ss01 * as01)

            end select

! --- determine tau, ssa, asy for each grid component
            tau = ext1
            if (ext1 > f_zero) ssa=min(f_one,ssa1/ext1)
            if (ssa1 > f_zero) asy=min(f_one,asy1/ssa1)

! --- save tau at 550 nm for each grid component
            if ( ib == nv_aod ) then
              do ijk = 1, max_num_gridcomp
                if ( tp == max_gridcomp(ijk) )  then
                   tauae_gocart(kk,ijk) = tau
                endif
              enddo
            endif

! --- update sum_tau, sum_ssa, sum_asy
            sum_tau = sum_tau + tau
            sum_ssa = sum_ssa + tau * ssa
            sum_asy = sum_asy + tau * ssa * asy

          enddo lab_do_icmp


! --- determine total tau, ssa, asy for aerosol mixture
          tauae(kk,ib) = sum_tau
          if (sum_tau > f_zero) ssaae(kk,ib) = sum_ssa / sum_tau
          if (sum_ssa > f_zero) asyae(kk,ib) = sum_asy / sum_ssa

        enddo   lab_do_ib

      enddo lab_do_layer

!
      return
!...................................
      end subroutine aeropt_grt
!! @}
!--------------------------------

!................................
      end subroutine setgocartaer
!--------------------------------
!! @}
!
! GOCART code modification end here (Sarah Lu)  ------------------------!
! =======================================================================

!..........................................!
      end module module_radiation_aerosols !
!==========================================!
