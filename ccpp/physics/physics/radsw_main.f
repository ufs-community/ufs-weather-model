!>  \file radsw_main.f
!!  This file contains NCEP's modifications of the rrtmg-sw radiation
!!  code from AER.

!  ==============================================================  !!!!!
!             sw-rrtm3 radiation package description              !!!!!
!  ==============================================================  !!!!!
!                                                                          !
!   this package includes ncep's modifications of the rrtm-sw radiation    !
!   code from aer inc.                                                     !
!                                                                          !
!   the sw-rrtm3 package includes these parts:                             !
!                                                                          !
!      'radsw_rrtm3_param.f'                                               !
!      'radsw_rrtm3_datatb.f'                                              !
!      'radsw_rrtm3_main.f'                                                !
!                                                                          !
!   the 'radsw_rrtm3_param.f' contains:                                    !
!                                                                          !
!      'module_radsw_parameters'  -- band parameters set up                !
!                                                                          !
!   the 'radsw_rrtm3_datatb.f' contains:                                   !
!                                                                          !
!      'module_radsw_ref'         -- reference temperature and pressure    !
!      'module_radsw_cldprtb'     -- cloud property coefficients table     !
!      'module_radsw_sflux'       -- spectral distribution of solar flux   !
!      'module_radsw_kgbnn'       -- absorption coeffients for 14          !
!                                    bands, where nn = 16-29               !
!                                                                          !
!   the 'radsw_rrtm3_main.f' contains:                                     !
!                                                                          !
!      'rrtmg_sw'                 -- main sw radiation transfer            !
!                                                                          !
!   in the main module 'rrtmg_sw' there are only two                       !
!   externally callable subroutines:                                       !
!                                                                          !
!      'swrad'      -- main sw radiation routine                           !
!         inputs:                                                          !
!           (plyr,plvl,tlyr,tlvl,qlyr,olyr,gasvmr,                         !
!            clouds,icseed,aerosols,sfcalb,                                !
!            cosz,solcon,NDAY,idxday,                                      !
!            npts, nlay, nlp1, lprnt,                                      !
!         outputs:                                                         !
!            hswc,topflx,sfcflx,                                           !
!!        optional outputs:                                                !
!            HSW0,HSWB,FLXPRF,FDNCMP)                                      !
!           )                                                              !
!                                                                          !
!      'rswinit'    -- initialization routine                              !
!         inputs:                                                          !
!           ( me )                                                         !
!         outputs:                                                         !
!           (none)                                                         !
!                                                                          !
!   all the sw radiation subprograms become contained subprograms          !
!   in module 'rrtmg_sw' and many of them are not directly                 !
!   accessable from places outside the module.                             !
!                                                                          !
!    derived data type constructs used:                                    !
!                                                                          !
!     1. radiation flux at toa: (from module 'module_radsw_parameters')    !
!          topfsw_type   -  derived data type for toa rad fluxes           !
!            upfxc              total sky upward flux at toa               !
!            dnfxc              total sky downward flux at toa             !
!            upfx0              clear sky upward flux at toa               !
!                                                                          !
!     2. radiation flux at sfc: (from module 'module_radsw_parameters')    !
!          sfcfsw_type   -  derived data type for sfc rad fluxes           !
!            upfxc              total sky upward flux at sfc               !
!            dnfxc              total sky downward flux at sfc             !
!            upfx0              clear sky upward flux at sfc               !
!            dnfx0              clear sky downward flux at sfc             !
!                                                                          !
!     3. radiation flux profiles(from module 'module_radsw_parameters')    !
!          profsw_type    -  derived data type for rad vertical prof       !
!            upfxc              level upward flux for total sky            !
!            dnfxc              level downward flux for total sky          !
!            upfx0              level upward flux for clear sky            !
!            dnfx0              level downward flux for clear sky          !
!                                                                          !
!     4. surface component fluxes(from module 'module_radsw_parameters'    !
!          cmpfsw_type    -  derived data type for component sfc flux      !
!            uvbfc              total sky downward uv-b flux at sfc        !
!            uvbf0              clear sky downward uv-b flux at sfc        !
!            nirbm              surface downward nir direct beam flux      !
!            nirdf              surface downward nir diffused flux         !
!            visbm              surface downward uv+vis direct beam flx    !
!            visdf              surface downward uv+vis diffused flux      !
!                                                                          !
!   external modules referenced:                                           !
!                                                                          !
!       'module physparam'                                                 !
!       'module physcons'                                                  !
!       'mersenne_twister'                                                 !
!                                                                          !
!   compilation sequence is:                                               !
!                                                                          !
!      'radsw_rrtm3_param.f'                                               !
!      'radsw_rrtm3_datatb.f'                                              !
!      'radsw_rrtm3_main.f'                                                !
!                                                                          !
!   and all should be put in front of routines that use sw modules         !
!                                                                          !
!==========================================================================!
!                                                                          !
!   the original program declarations:                                     !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
!  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  !
!  This software may be used, copied, or redistributed as long as it is    !
!  not sold and this copyright notice is reproduced on each copy made.     !
!  This model is provided as is without any express or implied warranties. !
!                       (http://www.rtweb.aer.com/)                        !
!                                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                          !
! ************************************************************************ !
!                                                                          !
!                              rrtmg_sw                                    !
!                                                                          !
!                                                                          !
!                   a rapid radiative transfer model                       !
!                    for the solar spectral region                         !
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
!        following people:  steven j. taubman, patrick d. brown,           !
!        ronald e. farren, luke chen, robert bergstrom.                    !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    references:                                                           !
!    (rrtm_sw/rrtmg_sw):                                                   !
!      clough, s.a., m.w. shephard, e.j. mlawer, j.s. delamere,            !
!      m.j. iacono, k. cady-pereira, s. boukabara, and p.d. brown:         !
!      atmospheric radiative transfer modeling: a summary of the aer       !
!      codes, j. quant. spectrosc. radiat. transfer, 91, 233-244, 2005.    !
!                                                                          !
!    (mcica):                                                              !
!      pincus, r., h. w. barker, and j.-j. morcrette: a fast, flexible,    !
!      approximation technique for computing radiative transfer in         !
!      inhomogeneous cloud fields, j. geophys. res., 108(d13), 4376,       !
!      doi:10.1029/2002jd003322, 2003.                                     !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!    aer's revision history:                                               !
!     this version of rrtmg_sw has been modified from rrtm_sw to use a     !
!     reduced set of g-point intervals and a two-stream model for          !
!     application to gcms.                                                 !
!                                                                          !
! --  original version (derived from rrtm_sw)                              !
!        2002: aer. inc.                                                   !
! --  conversion to f90 formatting; addition of 2-stream radiative transfer!
!        feb 2003: j.-j. morcrette, ecmwf                                  !
! --  additional modifications for gcm application                         !
!        aug 2003: m. j. iacono, aer inc.                                  !
! --  total number of g-points reduced from 224 to 112.  original          !
!     set of 224 can be restored by exchanging code in module parrrsw.f90  !
!     and in file rrtmg_sw_init.f90.                                       !
!        apr 2004: m. j. iacono, aer, inc.                                 !
! --  modifications to include output for direct and diffuse               !
!     downward fluxes.  there are output as "true" fluxes without          !
!     any delta scaling applied.  code can be commented to exclude         !
!     this calculation in source file rrtmg_sw_spcvrt.f90.                 !
!        jan 2005: e. j. mlawer, m. j. iacono, aer, inc.                   !
! --  revised to add mcica capability.                                     !
!        nov 2005: m. j. iacono, aer, inc.                                 !
! --  reformatted for consistency with rrtmg_lw.                           !
!        feb 2007: m. j. iacono, aer, inc.                                 !
! --  modifications to formatting to use assumed-shape arrays.             !
!        aug 2007: m. j. iacono, aer, inc.                                 !
!                                                                          !
! ************************************************************************ !
!                                                                          !
!   ncep modifications history log:                                        !
!                                                                          !
!       sep 2003,  yu-tai hou        -- received aer's rrtm-sw gcm version !
!                    code (v224)                                           !
!       nov 2003,  yu-tai hou        -- corrected errors in direct/diffuse !
!                    surface alabedo components.                           !
!       jan 2004,  yu-tai hou        -- modified code into standard modular!
!                    f9x code for ncep models. the original three cloud    !
!                    control flags are simplified into two: iflagliq and   !
!                    iflagice. combined the org subr sw_224 and setcoef    !
!                    into radsw (the main program); put all kgb##together  !
!                    and reformat into a separated data module; combine    !
!                    reftra and vrtqdr as swflux; optimized taumol and all !
!                    taubgs to form a contained subroutines.               !
!       jun 2004,  yu-tai hou        -- modified code based on aer's faster!
!                    version rrtmg_sw (v2.0) with 112 g-points.            !
!       mar 2005,  yu-tai hou        -- modified to aer v2.3, correct cloud!
!                    scaling error, total sky properties are delta scaled  !
!                    after combining clear and cloudy parts. the testing   !
!                    criterion of ssa is saved before scaling. added cloud !
!                    layer rain and snow contributions. all cloud water    !
!                    partical contents are treated the same way as other   !
!                    atmos particles.                                      !
!       apr 2005,  yu-tai hou        -- modified on module structures (this!
!                    version of code was given back to aer in jun 2006)    !
!       nov 2006,  yu-tai hou        -- modified code to include the       !
!                    generallized aerosol optical property scheme for gcms.!
!       apr 2007,  yu-tai hou        -- added spectral band heating as an  !
!                    optional output to support the 500km model's upper    !
!                    stratospheric radiation calculations. restructure     !
!                    optional outputs for easy access by different models. !
!       oct 2008,  yu-tai hou        -- modified to include new features   !
!                    from aer's newer release v3.5-v3.61, including mcica  !
!                    sub-grid cloud option and true direct/diffuse fluxes  !
!                    without delta scaling. added rain/snow opt properties !
!                    support to cloudy sky calculations. simplified and    !
!                    unified sw and lw sub-column cloud subroutines into   !
!                    one module by using optional parameters.              !
!       mar 2009,  yu-tai hou        -- replaced the original random number!
!                    generator coming with the original code with ncep w3  !
!                    library to simplify the program and moved sub-column  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       mar 2009,  yu-tai hou        -- replaced random number generator   !
!                    programs coming from the original code with the ncep  !
!                    w3 library to simplify the program and moved sub-col  !
!                    cloud subroutines inside the main module. added       !
!                    option of user provided permutation seeds that could  !
!                    be randomly generated from forecast time stamp.       !
!       nov 2009,  yu-tai hou        -- updated to aer v3.7-v3.8 version.  !
!                    notice the input cloud ice/liquid are assumed as      !
!                    in-cloud quantities, not grid average quantities.     !
!       aug 2010,  yu-tai hou        -- uptimized code to improve efficiency
!                    splited subroutine spcvrt into two subs, spcvrc and   !
!                    spcvrm, to handling non-mcica and mcica type of calls.!
!       apr 2012,  b. ferrier and y. hou -- added conversion factor to fu's!
!                    cloud-snow optical property scheme.                   !
!       jul 2012,  s. moorthi and Y. hou  -- eliminated the pointer array  !
!                     in subr 'spcvrt' for multi-threading issue running   !
!                     under intel's fortran compiler.                      !
!       nov 2012,  yu-tai hou        -- modified control parameters thru   !
!                     module 'physparam'.                                  !
!       jun 2013,  yu-tai hou        -- moving band 9 surface treatment    !
!                     back as in the rrtm2 version, spliting surface flux  !
!                     into two spectral regions (vis & nir), instead of    !
!                     designated it in nir region only.                    !
!       may 2016   yu-tai hou        --reverting swflux name back to vrtqdr!
!                                                                          !
!!!!!  ==============================================================  !!!!!
!!!!!                         end descriptions                         !!!!!
!!!!!  ==============================================================  !!!!!


!========================================!
      module rrtmg_sw                    !
!........................................!
!
      use physparam,        only : iswrate, iswrgas, iswcliq, iswcice,  &
     &                             isubcsw, icldflg, iovrsw,  ivflip,   &
     &                             iswmode, kind_phys
      use physcons,         only : con_g, con_cp, con_avgd, con_amd,    &
     &                             con_amw, con_amo3

      use module_radsw_parameters
      use mersenne_twister, only : random_setseed, random_number,       &
     &                             random_stat
      use module_radsw_ref, only : preflog, tref
      use module_radsw_sflux
!
      implicit none
!
      private
!
!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGSW='NCEP SW v5.1  Nov 2012 -RRTMG-SW v3.8   '
!    &   VTAGSW='NCEP SW v5.0  Aug 2012 -RRTMG-SW v3.8   '
!    &   VTAGSW='RRTMG-SW v3.8   Nov 2009'
!    &   VTAGSW='RRTMG-SW v3.7   Nov 2009'
!    &   VTAGSW='RRTMG-SW v3.61  Oct 2008'
!    &   VTAGSW='RRTMG-SW v3.5   Oct 2008'
!    &   VTAGSW='RRTM-SW 112v2.3 Apr 2007'
!    &   VTAGSW='RRTM-SW 112v2.3 Mar 2005'
!    &   VTAGSW='RRTM-SW 112v2.0 Jul 2004'

! \name constant values

      real (kind=kind_phys), parameter :: eps     = 1.0e-6
      real (kind=kind_phys), parameter :: oneminus= 1.0 - eps
! pade approx constant
      real (kind=kind_phys), parameter :: bpade   = 1.0/0.278
      real (kind=kind_phys), parameter :: stpfac  = 296.0/1013.0
      real (kind=kind_phys), parameter :: ftiny   = 1.0e-12
      real (kind=kind_phys), parameter :: flimit  = 1.0e-20
! internal solar constant
      real (kind=kind_phys), parameter :: s0      = 1368.22

      real (kind=kind_phys), parameter :: f_zero  = 0.0
      real (kind=kind_phys), parameter :: f_one   = 1.0

! \name atomic weights for conversion from mass to volume mixing ratios
      real (kind=kind_phys), parameter :: amdw    = con_amd/con_amw
      real (kind=kind_phys), parameter :: amdo3   = con_amd/con_amo3

! \name band indices
      integer, dimension(nblow:nbhgh) :: nspa, nspb
! band index for sfc flux
      integer, dimension(nblow:nbhgh) :: idxsfc
! band index for cld prop
      integer, dimension(nblow:nbhgh) :: idxebc

      data nspa(:) /  9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1 /
      data nspb(:) /  1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1 /

!     data idxsfc(:) / 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1 /  ! band index for sfc flux
      data idxsfc(:) / 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 2, 2, 2, 1 /  ! band index for sfc flux
      data idxebc(:) / 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1, 5 /  ! band index for cld prop

!  ---  band wavenumber intervals
!     real (kind=kind_phys), dimension(nblow:nbhgh):: wavenum1,wavenum2
!     data wavenum1(:)  /                                               &
!    &         2600.0, 3250.0, 4000.0, 4650.0, 5150.0, 6150.0, 7700.0,  &
!    &         8050.0,12850.0,16000.0,22650.0,29000.0,38000.0,  820.0 /
!     data wavenum2(:)  /                                               &
!              3250.0, 4000.0, 4650.0, 5150.0, 6150.0, 7700.0, 8050.0,  &
!    &        12850.0,16000.0,22650.0,29000.0,38000.0,50000.0, 2600.0 /
!     real (kind=kind_phys), dimension(nblow:nbhgh) :: delwave
!     data delwave(:)   /                                               &
!    &          650.0,  750.0,  650.0,  500.0, 1000.0, 1550.0,  350.0,  &
!    &         4800.0, 3150.0, 6650.0, 6350.0, 9000.0,12000.0, 1780.0 /

! uv-b band index
      integer, parameter :: nuvb = 27

!\name logical flags for optional output fields
      logical :: lhswb  = .false.
      logical :: lhsw0  = .false.
      logical :: lflxprf= .false.
      logical :: lfdncmp= .false.


! those data will be set up only once by "rswinit"
      real (kind=kind_phys) :: exp_tbl(0:NTBMX)


! the factor for heating rates (in k/day, or k/sec set by subroutine
!!  'rswinit')
      real (kind=kind_phys) :: heatfac


! initial permutation seed used for sub-column cloud scheme
      integer, parameter :: ipsdsw0 = 1

!  ---  public accessable subprograms

      public rrtmg_sw_init, rrtmg_sw_run, rrtmg_sw_finalize, rswinit 


! =================
      contains
! =================

      subroutine rrtmg_sw_init ()
      end subroutine rrtmg_sw_init

!> \ingroup RRTMG
!! \defgroup module_radsw_main GFS radsw Main
!! This module includes NCEP's modifications of the RRTMG-SW radiation
!! code from AER.
!!
!! The SW radiation model in the current NOAA Environmental Modeling
!! System (NEMS) was adapted from the RRTM radiation model developed by
!! AER Inc. (\cite clough_et_al_2005; \cite mlawer_et_al_1997). It contains 14
!! spectral bands spanning a spectral wavenumber range of
!! \f$50000-820 cm^{-1}\f$ (corresponding to a wavelength range
!! \f$0.2-12.2\mu m\f$), each spectral band focuses on a specific set of
!! atmospheric absorbing species as shown in Table 1. To achieve great
!! computation efficiency while at the same time to maintain a high
!! degree of accuracy, the RRTM radiation model employs a corrected-k
!! distribution method (i.e. mapping the highly spectral changing
!! absorption coefficient, k, into a monotonic and smooth varying
!! cumulative probability function, g). In the RRTM-SW, there are 16
!! unevenly distributed g points for each of the 14 bands for a total
!! of 224 g points. The GCM version of the code (RRTMG-SW) uses a reduced
!! number (various between 2 to 16) of g points for each of the bands
!! that totals to 112 instead of the full set of 224. To get high
!! quality for the scheme, many advanced techniques are used in RRTM
!! such as carefully selecting the band structure to handle various
!! major (key-species) and minor absorbers; deriving a binary parameter
!! for a paired key molecular species in the same domain; and using two
!! pressure regions (dividing level is at about 96mb) for optimal
!! treatment of various species, etc.
!!\tableofcontents
!! Table 1. RRTMG-SW spectral bands and the corresponding absorbing species
!! |Band #| Wavenumber Range | Lower Atm (Key)| Lower Atm (Minor)| Mid/Up Atm (Key)| Mid/Up Atm (Minor)|
!! |------|------------------|----------------|------------------|-----------------|-------------------|
!! |  16  |     2600-3250    |H2O,CH4         |                  |CH4              |                   |
!! |  17  |     3250-4000    |H2O,CO2         |                  |H2O,CO2          |                   |
!! |  18  |     4000-4650    |H2O,CH4         |                  |CH4              |                   |
!! |  19  |     4650-5150    |H2O,CO2         |                  |CO2              |                   |
!! |  20  |     5150-6150    |H2O             |CH4               |H2O              |CH4                |
!! |  21  |     6150-7700    |H2O,CO2         |                  |H2O,CO2          |                   |
!! |  22  |     7700-8050    |H2O,O2          |                  |O2               |                   |
!! |  23  |     8050-12850   |H2O             |                  |---              |                   |
!! |  24  |    12850-16000   |H2O,O2          |O3                |O2               |O3                 |
!! |  25  |    16000-22650   |H2O             |O3                |---              |O3                 |
!! |  26  |    22650-29000   |---             |                  |---              |                   |
!! |  27  |    29000-38000   |O3              |                  |O3               |                   |
!! |  28  |    38000-50000   |O3,O2           |                  |O3,O2            |                   |
!! |  29  |      820-2600    |H2O             |CO2               |CO2              |H2O                |
!!\tableofcontents
!!\n scattering due to clouds greatly complicate the SW radiative
!! transfer computations. To balance the trade-off between computation
!! and speed, RRTMG-SW uses a two-stream approximation method with a
!! delta-function adjustment. Several variations of the delta-two
!! method are included in the radiation transfer code; each holds its
!! own strength and shortcomings (\cite king_and_harshvardhan_1986;
!! \cite raisanen_2002; \cite barker_et_al_2015). The default (the same
!! in operation runs) selection (iswmode=2) activates the Practical
!! Improved Flux Method (PIFM) by 
!! \cite zdunkowski_et_al_1980 . In dealing with a column of cloudy
!! atmosphere, two approaches are included in the RRTMG-SW. One is the
!! commonly used treatment that sees each of the cloud contaminated
!! layers as independent, partially and uniformly filled slabs. Cloud
!! inhomogeneity within and the nature coherence among adjacent cloud
!! layers are largely ignored to reduce the overwhelm complexities
!! associated with scattering process. The effective layer reflectance
!! and transmittance are weighted mean according to cloud fraction. The
!! approach may overestimate cloud effect, especially for multi-layered
!! cloud system associated with deep convection. In NEMS radiation code,
!! to mitigate this shortcoming without increase computation cost, the
!! cloud contaminated column is divided into two parts based on the
!! column's total cloud coverage (a maximum-random overlapping is used
!! in the operational models) to form a cloud free part and an overcast
!! part. Layered clouds are then normalized by the total cloud amount
!! before going through radiative transfer calculations. Fluxes from the
!! cloud-free part and cloudy part are combined together to obtain the
!! final result.
!!\n On the other hand, the Monte-Carlo Independent Column Approximation
!! (McICA) (\cite pincus_et_al_2003;
!! \cite raisanen_and_barker_2004), provides a simple and effective way
!! to solve cloud overlapping issue without increasing computational
!! burden. The method is based on the concept of an ICA scheme that
!! divides each grid column into a large number of sub-columns, and
!! statistically redistributes layered clouds (under an assumed overlapping
!! condition, such as the maximum-random method) into the sub-columns
!! (i.e. at any layer it will be either clear or overcast). Thus the
!! grid domain averaged flux under ICA scheme can be expressed as:
!! \f[
!! \overline{F}=\frac{1}{N}\sum_{n=1}^N F_{n}
!! =\frac{1}{N}\sum_{n=1}^N\sum_{k=1}^K F_{n,k}
!! \f]
!! Where \f$N\f$ is the number of total sub-columns, and \f$K\f$ is the
!! number of spectral terms in integration.\f$F_{n}\f$ is flux obtained
!! in the \f$n^{th}\f$ sub-column, that is the summation of total of
!! \f$K\f$ spectral corresponding fluxes, \f$F_{n,k}\f$ . The double
!! integrations (summations) make ICA impractical for GCM applications.
!! The McICA method is to divide a model grid into \f$K\f$ sub-columns
!! and randomly to pair a sub-column's cloud profile with one of the
!! radiative spectral intervals (e.g. the g-point in RRTM). The double
!! summations will then be reduced to only one:
!! \f[
!! \overline{F}=\frac{1}{N}\sum_{n=1}^N\sum_{k=1}^K F_{n,k}
!! \approx\overline{F}=\sum_{k=1}^K F_{S_{k},k}
!! \f]
!!
!! The RRTM-SW package includes three files:
!! - radsw_param.f, which contains:
!!  - module_radsw_parameters: specifies major parameters of the spectral
!!    bands and defines the construct structures of derived-type variables
!!    for holding the output results.
!! - radsw_datatb.f, which contains:
!!  - module_radsw_ref: reference temperature and pressure
!!  - module_radsw_cldprtb: cloud property coefficients table
!!  - module_radsw_sflux: indexes and coefficients for spectral
!!                        distribution of solar flux
!!  - module_radsw_kgbnn: absorption coefficents for 14 bands, where
!!    nn = 16-29
!! - radsw_main.f, which contains:
!!  - module_radsw_main: the main SW radiation computation source codes
!!
!!\author   Eli J. Mlawer, emlawer@aer.com
!!\author   Jennifer S. Delamere, jdelamer@aer.com
!!\author   Michael J. Iacono, miacono@aer.com
!!\author   Shepard A. Clough
!!\version NCEP SW v5.1  Nov 2012 -RRTMG-SW v3.8
!!
!! The authors wish to acknowledge the contributions of the
!! following people:  Steven J. Taubman, Karen Cady-Pereira,
!! Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert Bergstrom.
!!
!!\copyright  2002-2007, Atmospheric & Environmental Research, Inc. (AER).
!!  This software may be used, copied, or redistributed as long as it is
!!  not sold and this copyright notice is reproduced on each copy made.
!!  This model is provided as is without any express or implied warranties.
!!  (http://www.rtweb.aer.com/)
!! \section arg_table_rrtmg_sw_run Argument Table
!! | local_name      | standard_name                                                                                  | long_name                                                                | units   | rank | type        |    kind   | intent | optional |
!! |-----------------|------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|---------|------|-------------|-----------|--------|----------|
!! | plyr            | air_pressure_at_layer_for_radiation_in_hPa                                                     | air pressure layer                                                       | hPa     |    2 | real        | kind_phys | in     | F        |
!! | plvl            | air_pressure_at_interface_for_radiation_in_hPa                                                 | air pressure level                                                       | hPa     |    2 | real        | kind_phys | in     | F        |
!! | tlyr            | air_temperature_at_layer_for_radiation                                                         | air temperature layer                                                    | K       |    2 | real        | kind_phys | in     | F        |
!! | tlvl            | air_temperature_at_interface_for_radiation                                                     | air temperature level                                                    | K       |    2 | real        | kind_phys | in     | F        |
!! | qlyr            | water_vapor_specific_humidity_at_layer_for_radiation                                           | specific humidity layer                                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | olyr            | ozone_concentration_at_layer_for_radiation                                                     | ozone concentration layer                                                | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_co2      | volume_mixing_ratio_co2                                                                        | volume mixing ratio co2                                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_n2o      | volume_mixing_ratio_n2o                                                                        | volume mixing ratio no2                                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_ch4      | volume_mixing_ratio_ch4                                                                        | volume mixing ratio ch4                                                  | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_o2       | volume_mixing_ratio_o2                                                                         | volume mixing ratio o2                                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_co       | volume_mixing_ratio_co                                                                         | volume mixing ratio co                                                   | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_cfc11    | volume_mixing_ratio_cfc11                                                                      | volume mixing ratio cfc11                                                | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_cfc12    | volume_mixing_ratio_cfc12                                                                      | volume mixing ratio cfc12                                                | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_cfc22    | volume_mixing_ratio_cfc22                                                                      | volume mixing ratio cfc22                                                | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | gasvmr_ccl4     | volume_mixing_ratio_ccl4                                                                       | volume mixing ratio ccl4                                                 | kg kg-1 |    2 | real        | kind_phys | in     | F        |
!! | icseed          | seed_random_numbers_sw                                                                         | seed for random number generation for shortwave radiation                | none    |    1 | integer     |           | in     | F        |
!! | aeraod          | aerosol_optical_depth_for_shortwave_bands_01-16                                                | aerosol optical depth for shortwave bands 01-16                          | none    |    3 | real        | kind_phys | in     | F        |
!! | aerssa          | aerosol_single_scattering_albedo_for_shortwave_bands_01-16                                     | aerosol single scattering albedo for shortwave bands 01-16               | frac    |    3 | real        | kind_phys | in     | F        |
!! | aerasy          | aerosol_asymmetry_parameter_for_shortwave_bands_01-16                                          | aerosol asymmetry paramter for shortwave bands 01-16                     | none    |    3 | real        | kind_phys | in     | F        |
!! | sfcalb_nir_dir  | surface_albedo_due_to_near_IR_direct                                                           | surface albedo due to near IR direct beam                                | frac    |    1 | real        | kind_phys | in     | F        |
!! | sfcalb_nir_dif  | surface_albedo_due_to_near_IR_diffused                                                         | surface albedo due to near IR diffused beam                              | frac    |    1 | real        | kind_phys | in     | F        |
!! | sfcalb_uvis_dir | surface_albedo_due_to_UV_and_VIS_direct                                                        | surface albedo due to UV+VIS direct beam                                 | frac    |    1 | real        | kind_phys | in     | F        |
!! | sfcalb_uvis_dif | surface_albedo_due_to_UV_and_VIS_diffused                                                      | surface albedo due to UV+VIS diffused beam                               | frac    |    1 | real        | kind_phys | in     | F        |
!! | cosz            | cosine_of_zenith_angle                                                                         | cosine of the solar zenit angle                                          | none    |    1 | real        | kind_phys | in     | F        |
!! | solcon          | solar_constant                                                                                 | solar constant                                                           | W m-2   |    0 | real        | kind_phys | in     | F        |
!! | nday            | daytime_points_dimension                                                                       | daytime points dimension                                                 | count   |    0 | integer     |           | in     | F        |
!! | idxday          | daytime_points                                                                                 | daytime points                                                           | index   |    1 | integer     |           | in     | F        |
!! | npts            | horizontal_loop_extent                                                                         | horizontal dimension                                                     | count   |    0 | integer     |           | in     | F        |
!! | nlay            | adjusted_vertical_layer_dimension_for_radiation                                                | number of vertical layers for radiation                                  | count   |    0 | integer     |           | in     | F        |
!! | nlp1            | adjusted_vertical_level_dimension_for_radiation                                                | number of vertical levels for radiation                                  | count   |    0 | integer     |           | in     | F        |
!! | lprnt           | flag_print                                                                                     | flag to print                                                            | flag    |    0 | logical     |           | in     | F        |
!! | cld_cf          | total_cloud_fraction                                                                           | total cloud fraction                                                     | frac    |    2 | real        | kind_phys | in     | F        |
!! | lsswr           | flag_to_calc_sw                                                                                | flag to calculate SW irradiances                                         | flag    |    0 | logical     |           | in     | F        |
!! | hswc            | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | shortwave total sky heating rate                                         | K s-1   |    2 | real        | kind_phys | inout  | F        |
!! | topflx          | sw_fluxes_top_atmosphere                                                                       | shortwave total sky fluxes at the top of the atm                         | W m-2   |    1 | topfsw_type |           | inout  | F        |
!! | sfcflx          | sw_fluxes_sfc                                                                                  | shortwave total sky fluxes at the Earth surface                          | W m-2   |    1 | sfcfsw_type |           | inout  | F        |
!! | hsw0            | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | shortwave clear sky heating rate                                         | K s-1   |    2 | real        | kind_phys | inout  | T        |
!! | hswb            | sw_heating_rate_spectral                                                                       | shortwave total sky heating rate (spectral)                              | K s-1   |    3 | real        | kind_phys | inout  | T        |
!! | flxprf          | sw_fluxes                                                                                      | sw fluxes total sky / csk and up / down at levels                        | W m-2   |    2 | profsw_type |           | inout  | T        |
!! | fdncmp          | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes | W m-2   |    1 | cmpfsw_type |           | inout  | T        |
!! | cld_lwp         | cloud_liquid_water_path                                                                        | cloud liquid water path                                                  | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_liq     | mean_effective_radius_for_liquid_cloud                                                         | mean effective radius for liquid cloud                                   | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_iwp         | cloud_ice_water_path                                                                           | cloud ice water path                                                     | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_ice     | mean_effective_radius_for_ice_cloud                                                            | mean effective radius for ice cloud                                      | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_rwp         | cloud_rain_water_path                                                                          | cloud rain water path                                                    | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_rain    | mean_effective_radius_for_rain_drop                                                            | mean effective radius for rain drop                                      | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_swp         | cloud_snow_water_path                                                                          | cloud snow water path                                                    | g m-2   |    2 | real        | kind_phys | in     | T        |
!! | cld_ref_snow    | mean_effective_radius_for_snow_flake                                                           | mean effective radius for snow flake                                     | micron  |    2 | real        | kind_phys | in     | T        |
!! | cld_od          | cloud_optical_depth                                                                            | cloud optical depth                                                      | none    |    2 | real        | kind_phys | in     | T        |
!! | cld_ssa         | cloud_single_scattering_albedo                                                                 | cloud single scattering albedo                                           | frac    |    2 | real        | kind_phys | in     | T        |
!! | cld_asy         | cloud_asymmetry_parameter                                                                      | cloud asymmetry parameter                                                | none    |    2 | real        | kind_phys | in     | T        |
!! | errmsg          | error_message                                                                                  | error message for error handling in CCPP                                 | none    |    0 | character   | len=*     | out    | F        |
!! | errflg          | error_flag                                                                                     | error flag for error handling in CCPP                                    | flag    |    0 | integer     |           | out    | F        |
!!
!> \section gen_swrad RRTMG Shortwave Radiation Scheme General Algorithm
!> @{
!-----------------------------------
      subroutine rrtmg_sw_run                                           &
     &     ( plyr,plvl,tlyr,tlvl,qlyr,olyr,                             &
     &       gasvmr_co2,gasvmr_n2o,gasvmr_ch4,gasvmr_o2,gasvmr_co,      &
     &       gasvmr_cfc11,gasvmr_cfc12,gasvmr_cfc22,gasvmr_ccl4,        &   !  ---  inputs
     &       icseed, aeraod, aerssa, aerasy,                            &
     &       sfcalb_nir_dir, sfcalb_nir_dif,                            &
     &       sfcalb_uvis_dir, sfcalb_uvis_dif,                          &
     &       cosz,solcon,NDAY,idxday,                                   &
     &       npts, nlay, nlp1, lprnt,                                   &
     &       cld_cf, lsswr,                                             &
     &       hswc,topflx,sfcflx,                                        &   !  ---  outputs
     &       HSW0,HSWB,FLXPRF,FDNCMP,                                   &   ! ---  optional
     &       cld_lwp, cld_ref_liq, cld_iwp, cld_ref_ice,                &
     &       cld_rwp,cld_ref_rain, cld_swp, cld_ref_snow,               &
     &       cld_od, cld_ssa, cld_asy, errmsg, errflg
     &     )

!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!   plyr (npts,nlay) : model layer mean pressure in mb                  !
!   plvl (npts,nlp1) : model level pressure in mb                       !
!   tlyr (npts,nlay) : model layer mean temperature in k                !
!   tlvl (npts,nlp1) : model level temperature in k    (not in use)     !
!   qlyr (npts,nlay) : layer specific humidity in gm/gm   *see inside   !
!   olyr (npts,nlay) : layer ozone concentration in gm/gm               !
!   gasvmr(npts,nlay,:): atmospheric constent gases:                    !
!                      (check module_radiation_gases for definition)    !
!      gasvmr(:,:,1)  - co2 volume mixing ratio                         !
!      gasvmr(:,:,2)  - n2o volume mixing ratio                         !
!      gasvmr(:,:,3)  - ch4 volume mixing ratio                         !
!      gasvmr(:,:,4)  - o2  volume mixing ratio                         !
!      gasvmr(:,:,5)  - co  volume mixing ratio        (not used)       !
!      gasvmr(:,:,6)  - cfc11 volume mixing ratio      (not used)       !
!      gasvmr(:,:,7)  - cfc12 volume mixing ratio      (not used)       !
!      gasvmr(:,:,8)  - cfc22 volume mixing ratio      (not used)       !
!      gasvmr(:,:,9)  - ccl4  volume mixing ratio      (not used)       !
!   clouds(npts,nlay,:): cloud profile                                  !
!                      (check module_radiation_clouds for definition)   !
!                ---  for  iswcliq > 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer in-cloud liq water path   (g/m**2)     !
!       clouds(:,:,3)  -   mean eff radius for liq cloud   (micron)     !
!       clouds(:,:,4)  -   layer in-cloud ice water path   (g/m**2)     !
!       clouds(:,:,5)  -   mean eff radius for ice cloud   (micron)     !
!       clouds(:,:,6)  -   layer rain drop water path      (g/m**2)     !
!       clouds(:,:,7)  -   mean eff radius for rain drop   (micron)     !
!       clouds(:,:,8)  -   layer snow flake water path     (g/m**2)     !
!       clouds(:,:,9)  -   mean eff radius for snow flake  (micron)     !
!                ---  for  iswcliq = 0  ---                             !
!       clouds(:,:,1)  -   layer total cloud fraction                   !
!       clouds(:,:,2)  -   layer cloud optical depth                    !
!       clouds(:,:,3)  -   layer cloud single scattering albedo         !
!       clouds(:,:,4)  -   layer cloud asymmetry factor                 !
!     icseed(npts)   : auxiliary special cloud related array            !
!                      when module variable isubcsw=2, it provides      !
!                      permutation seed for each column profile that    !
!                      are used for generating random numbers.          !
!                      when isubcsw /=2, it will not be used.           !
!   aerosols(npts,nlay,nbdsw,:) : aerosol optical properties            !
!                      (check module_radiation_aerosols for definition) !
!         (:,:,:,1)   - optical depth                                   !
!         (:,:,:,2)   - single scattering albedo                        !
!         (:,:,:,3)   - asymmetry parameter                             !
!   sfcalb(npts, : ) : surface albedo in fraction                       !
!                      (check module_radiation_surface for definition)  !
!         ( :, 1 )    - near ir direct beam albedo                      !
!         ( :, 2 )    - near ir diffused albedo                         !
!         ( :, 3 )    - uv+vis direct beam albedo                       !
!         ( :, 4 )    - uv+vis diffused albedo                          !
!   cosz  (npts)     : cosine of solar zenith angle                     !
!   solcon           : solar constant                      (w/m**2)     !
!   NDAY             : num of daytime points                            !
!   idxday(npts)     : index array for daytime points                   !
!   npts             : number of horizontal points                      !
!   nlay,nlp1        : vertical layer/lavel numbers                     !
!   lprnt            : logical check print flag                         !
!                                                                       !
!  output variables:                                                    !
!   hswc  (npts,nlay): total sky heating rates (k/sec or k/day)         !
!   topflx(npts)     : radiation fluxes at toa (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at toa                   !
!     dnflx            - total sky downward flux at toa                 !
!     upfx0            - clear sky upward flux at toa                   !
!   sfcflx(npts)     : radiation fluxes at sfc (w/m**2), components:    !
!                      (check module_radsw_parameters for definition)   !
!     upfxc            - total sky upward flux at sfc                   !
!     dnfxc            - total sky downward flux at sfc                 !
!     upfx0            - clear sky upward flux at sfc                   !
!     dnfx0            - clear sky downward flux at sfc                 !
!                                                                       !
!!optional outputs variables:                                           !
!   hswb(npts,nlay,nbdsw): spectral band total sky heating rates        !
!   hsw0  (npts,nlay): clear sky heating rates (k/sec or k/day)         !
!   flxprf(npts,nlp1): level radiation fluxes (w/m**2), components:     !
!                      (check module_radsw_parameters for definition)   !
!     dnfxc            - total sky downward flux at interface           !
!     upfxc            - total sky upward flux at interface             !
!     dnfx0            - clear sky downward flux at interface           !
!     upfx0            - clear sky upward flux at interface             !
!   fdncmp(npts)     : component surface downward fluxes (w/m**2):      !
!                      (check module_radsw_parameters for definition)   !
!     uvbfc            - total sky downward uv-b flux at sfc            !
!     uvbf0            - clear sky downward uv-b flux at sfc            !
!     nirbm            - downward surface nir direct beam flux          !
!     nirdf            - downward surface nir diffused flux             !
!     visbm            - downward surface uv+vis direct beam flux       !
!     visdf            - downward surface uv+vis diffused flux          !
!                                                                       !
!  external module variables:  (in physparam)                           !
!   iswrgas - control flag for rare gases (ch4,n2o,o2, etc.)            !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   iswcliq - control flag for liq-cloud optical properties             !
!           =0: input cloud optical depth, fixed ssa, asy               !
!           =1: use hu and stamnes(1993) method for liq cld             !
!           =2: not used                                                !
!   iswcice - control flag for ice-cloud optical properties             !
!           *** if iswcliq==0, iswcice is ignored                       !
!           =1: use ebert and curry (1992) scheme for ice clouds        !
!           =2: use streamer v3.0 (2001) method for ice clouds          !
!           =3: use fu's method (1996) for ice clouds                   !
!   iswmode - control flag for 2-stream transfer scheme                 !
!           =1; delta-eddington    (joseph et al., 1976)                !
!           =2: pifm               (zdunkowski et al., 1980)            !
!           =3: discrete ordinates (liou, 1973)                         !
!   isubcsw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   iovrsw  - cloud overlapping control flag                            !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud                                   !
!   ivflip  - control flg for direction of vertical index               !
!           =0: index from toa to surface                               !
!           =1: index from surface to toa                               !
!                                                                       !
!  module parameters, control variables:                                !
!     nblow,nbhgh      - lower and upper limits of spectral bands       !
!     maxgas           - maximum number of absorbing gaseous            !
!     ngptsw           - total number of g-point subintervals           !
!     ng##             - number of g-points in band (##=16-29)          !
!     ngb(ngptsw)      - band indices for each g-point                  !
!     bpade            - pade approximation constant (1/0.278)          !
!     nspa,nspb(nblow:nbhgh)                                            !
!                      - number of lower/upper ref atm's per band       !
!     ipsdsw0          - permutation seed for mcica sub-col clds        !
!                                                                       !
!  major local variables:                                               !
!     pavel  (nlay)         - layer pressures (mb)                      !
!     delp   (nlay)         - layer pressure thickness (mb)             !
!     tavel  (nlay)         - layer temperatures (k)                    !
!     coldry (nlay)         - dry air column amount                     !
!                                   (1.e-20*molecules/cm**2)            !
!     cldfrc (nlay)         - layer cloud fraction (norm by tot cld)    !
!     cldfmc (nlay,ngptsw)  - layer cloud fraction for g-point          !
!     taucw  (nlay,nbdsw)   - cloud optical depth                       !
!     ssacw  (nlay,nbdsw)   - cloud single scattering albedo (weighted) !
!     asycw  (nlay,nbdsw)   - cloud asymmetry factor         (weighted) !
!     tauaer (nlay,nbdsw)   - aerosol optical depths                    !
!     ssaaer (nlay,nbdsw)   - aerosol single scattering albedo          !
!     asyaer (nlay,nbdsw)   - aerosol asymmetry factor                  !
!     colamt (nlay,maxgas)  - column amounts of absorbing gases         !
!                             1 to maxgas are for h2o, co2, o3, n2o,    !
!                             ch4, o2, co, respectively (mol/cm**2)     !
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
!     laytrop               - layer at which switch is made from one    !
!                             combination of key species to another     !
!     jp(nlay),jt(nlay),jt1(nlay)                                       !
!                           - lookup table indexes                      !
!     flxucb(nlp1,nbdsw)    - spectral bnd total-sky upward flx (w/m2)  !
!     flxdcb(nlp1,nbdsw)    - spectral bnd total-sky downward flx (w/m2)!
!     flxu0b(nlp1,nbdsw)    - spectral bnd clear-sky upward flx (w/m2)  !
!     flxd0b(nlp1,nbdsw)    - spectral b d clear-sky downward flx (w/m2)!
!                                                                       !
!                                                                       !
!  =====================    end of definitions    ====================  !

!  ---  inputs:
      integer, intent(in) :: npts, nlay, nlp1, NDAY

      integer, dimension(:), intent(in) :: idxday, icseed

      logical, intent(in) :: lprnt, lsswr

      real (kind=kind_phys), dimension(npts,nlp1), intent(in) ::        &
     &       plvl, tlvl
      real (kind=kind_phys), dimension(npts,nlay), intent(in) ::        &
     &       plyr, tlyr, qlyr, olyr

      real (kind=kind_phys),dimension(npts),intent(in):: sfcalb_nir_dir &
      real (kind=kind_phys),dimension(npts),intent(in):: sfcalb_nir_dif &
      real (kind=kind_phys),dimension(npts),intent(in):: sfcalb_uvis_dir&
      real (kind=kind_phys),dimension(npts),intent(in):: sfcalb_uvis_dif&

      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_co2
      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_n2o
      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_ch4
      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_o2
      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_co
      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_cfc11
      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_cfc12
      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_cfc22
      real(kind=kind_phys),dimension(npts,nlay),intent(in)::gasvmr_ccl4

      real (kind=kind_phys), dimension(npts,nlay),intent(in):: cld_cf
      real (kind=kind_phys), dimension(npts,nlay),intent(in),optional:: &
     &       cld_lwp, cld_ref_liq,  cld_iwp, cld_ref_ice,               &
     &       cld_rwp, cld_ref_rain, cld_swp, cld_ref_snow,              &
     &       cld_od, cld_ssa, cld_asy

      real(kind=kind_phys),dimension(npts,nlay,nbdsw),intent(in)::aeraod
      real(kind=kind_phys),dimension(npts,nlay,nbdsw),intent(in)::aerssa
      real(kind=kind_phys),dimension(npts,nlay,nbdsw),intent(in)::aerasy

      real (kind=kind_phys), intent(in) :: cosz(npts), solcon

!  ---  outputs:
      real (kind=kind_phys), dimension(npts,nlay), intent(inout) :: hswc

      type (topfsw_type),    dimension(npts), intent(inout) :: topflx
      type (sfcfsw_type),    dimension(npts), intent(inout) :: sfcflx

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!! ---  optional outputs:
      real (kind=kind_phys), dimension(npts,nlay,nbdsw), optional,      &
     &       intent(inout) :: hswb

      real (kind=kind_phys), dimension(npts,nlay),       optional,      &
     &       intent(inout) :: hsw0
      type (profsw_type),    dimension(npts,nlp1),       optional,      &
     &       intent(inout) :: flxprf
      type (cmpfsw_type),    dimension(npts),            optional,      &
     &       intent(inout) :: fdncmp

!  ---  locals:
      real (kind=kind_phys), dimension(nlay,ngptsw) :: cldfmc,          &
     &       taug, taur
      real (kind=kind_phys), dimension(nlp1,nbdsw):: fxupc, fxdnc,      &
     &       fxup0, fxdn0

      real (kind=kind_phys), dimension(nlay,nbdsw)  ::                  &
     &       tauae, ssaae, asyae, taucw, ssacw, asycw

      real (kind=kind_phys), dimension(ngptsw) :: sfluxzen

      real (kind=kind_phys), dimension(nlay)   :: cldfrc,     delp,     &
     &       pavel, tavel, coldry, colmol, h2ovmr, o3vmr, temcol,       &
     &       cliqp, reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4,    &
     &       cfrac, fac00, fac01, fac10, fac11, forfac, forfrac,        &
     &       selffac, selffrac, rfdelp

      real (kind=kind_phys), dimension(nlp1) :: fnet, flxdc, flxuc,     &
     &       flxd0, flxu0

      real (kind=kind_phys), dimension(2) :: albbm, albdf, sfbmc,       &
     &       sfbm0, sfdfc, sfdf0

      real (kind=kind_phys) :: cosz1, sntz1, tem0, tem1, tem2, s0fac,   &
     &       ssolar, zcf0, zcf1, ftoau0, ftoauc, ftoadc,                &
     &       fsfcu0, fsfcuc, fsfcd0, fsfcdc, suvbfc, suvbf0

!  ---  column amount of absorbing gases:
!       (:,m) m = 1-h2o, 2-co2, 3-o3, 4-n2o, 5-ch4, 6-o2, 7-co
      real (kind=kind_phys) ::  colamt(nlay,maxgas)

      integer, dimension(npts) :: ipseed
      integer, dimension(nlay) :: indfor, indself, jp, jt, jt1

      integer :: i, ib, ipt, j1, k, kk, laytrop, mb

!
!===> ... begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
      if (.not. lsswr) return
      if (nday <= 0) return

      lhswb  = present ( hswb )
      lhsw0  = present ( hsw0 )
      lflxprf= present ( flxprf )
      lfdncmp= present ( fdncmp )

!> -# Compute solar constant adjustment factor (s0fac) according to solcon.
!      ***  s0, the solar constant at toa in w/m**2, is hard-coded with
!           each spectra band, the total flux is about 1368.22 w/m**2.

      s0fac = solcon / s0

!> -# Initial output arrays (and optional) as zero.

      hswc(:,:) = f_zero
      topflx = topfsw_type ( f_zero, f_zero, f_zero )
      sfcflx = sfcfsw_type ( f_zero, f_zero, f_zero, f_zero )

!! --- ...  initial optional outputs
      if ( lflxprf ) then
        flxprf = profsw_type ( f_zero, f_zero, f_zero, f_zero )
      endif

      if ( lfdncmp ) then
        fdncmp = cmpfsw_type (f_zero,f_zero,f_zero,f_zero,f_zero,f_zero)
      endif

      if ( lhsw0 ) then
        hsw0(:,:) = f_zero
      endif

      if ( lhswb ) then
        hswb(:,:,:) = f_zero
      endif

!! --- check for optional input arguments, depending on cloud method
      if (iswcliq > 0) then    ! use prognostic cloud method
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
        if ( .not.present(cld_od) .or. .not.present(cld_ssa) .or.       &
     &                                 .not.present(cld_asy)) then
          write(errmsg,'(*(a))')                                        &
     &               'Logic error: iswcliq<=0 requires the following',  &
     &               ' optional arguments to be present:',              &
     &               ' cld_od, cld_ssa, cld_asy'
          errflg = 1
          return
        end if
      endif                    ! end if_iswcliq

!> -# Change random number seed value for each radiation invocation
!!    (isubcsw =1 or 2).

      if     ( isubcsw == 1 ) then     ! advance prescribed permutation seed
        do i = 1, npts
          ipseed(i) = ipsdsw0 + i
        enddo
      elseif ( isubcsw == 2 ) then     ! use input array of permutaion seeds
        do i = 1, npts
          ipseed(i) = icseed(i)
        enddo
      endif

      if ( lprnt ) then
        write(0,*)'  In radsw, isubcsw, ipsdsw0,ipseed =',              &
     &           isubcsw, ipsdsw0, ipseed
      endif

!  --- ...  loop over each daytime grid point

      lab_do_ipt : do ipt = 1, NDAY

        j1 = idxday(ipt)

        cosz1  = cosz(j1)
        sntz1  = f_one / cosz(j1)
        ssolar = s0fac * cosz(j1)

!> -# Prepare surface albedo: bm,df - dir,dif; 1,2 - nir,uvv.
        albbm(1) = sfcalb_nir_dir(j1)
        albdf(1) = sfcalb_nir_dif(j1)
        albbm(2) = sfcalb_uvis_dir(j1)
        albdf(2) = sfcalb_uvis_dif(j1)

!> -# Prepare atmospheric profile for use in rrtm.
!           the vertical index of internal array is from surface to top

        if (ivflip == 0) then       ! input from toa to sfc

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd

          do k = 1, nlay
            kk = nlp1 - k
            pavel(k) = plyr(j1,kk)
            tavel(k) = tlyr(j1,kk)
            delp (k) = plvl(j1,kk+1) - plvl(j1,kk)
!> -# Set absorber and gas column amount, convert from volume mixing
!!    ratio to molec/cm2 based on coldry (scaled to 1.0e-20)
!!    - colamt(nlay,maxgas):column amounts of absorbing gases 1 to
!!      maxgas are for h2o,co2,o3,n2o,ch4,o2,co, respectively
!!      (\f$ mol/cm^2 \f$)

!test use
!           h2ovmr(k)= max(f_zero,qlyr(j1,kk)*amdw)                     ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qlyr(j1,kk))                          ! input vol mixing ratio
!           o3vmr (k)= max(f_zero,olyr(j1,kk))                          ! input vol mixing ratio
!ncep model use
            h2ovmr(k)= max(f_zero,qlyr(j1,kk)*amdw/(f_one-qlyr(j1,kk))) ! input specific humidity
            o3vmr (k)= max(f_zero,olyr(j1,kk)*amdo3)                    ! input mass mixing ratio

            tem0 = (f_one - h2ovmr(k))*con_amd + h2ovmr(k)*con_amw
            coldry(k) = tem2 * delp(k) / (tem1*tem0*(f_one + h2ovmr(k)))
            temcol(k) = 1.0e-12 * coldry(k)

            colamt(k,1) = max(f_zero,    coldry(k)*h2ovmr(k))         ! h2o
            !colamt(k,2) = max(temcol(k), coldry(k)*gasvmr(j1,kk,1))   ! co2
            colamt(k,2) = max(temcol(k), coldry(k)*gasvmr_co2(j1,kk))   ! co2
            colamt(k,3) = max(f_zero,    coldry(k)*o3vmr(k))          ! o3
            colmol(k)   = coldry(k) + colamt(k,1)
          enddo

!  --- ...  set up gas column amount, convert from volume mixing ratio
!           to molec/cm2 based on coldry (scaled to 1.0e-20)

          if (iswrgas > 0) then
            do k = 1, nlay
              kk = nlp1 - k
              !colamt(k,4) = max(temcol(k), coldry(k)*gasvmr(j1,kk,2))  ! n2o
             colamt(k,4) = max(temcol(k), coldry(k)*gasvmr_n2o(j1,kk))  ! n2o
             ! colamt(k,5) = max(temcol(k), coldry(k)*gasvmr(j1,kk,3))  ! ch4
             colamt(k,5) = max(temcol(k), coldry(k)*gasvmr_ch4(j1,kk))  ! ch4
             ! colamt(k,6) = max(temcol(k), coldry(k)*gasvmr(j1,kk,4))  ! o2
             colamt(k,6) = max(temcol(k), coldry(k)*gasvmr_o2(j1,kk))  ! o2
!             colamt(k,7) = max(temcol(k), coldry(k)*gasvmr(j1,kk,5))  ! co - notused
            enddo
          else
            do k = 1, nlay
              colamt(k,4) = temcol(k)                                  ! n2o
              colamt(k,5) = temcol(k)                                  ! ch4
              colamt(k,6) = temcol(k)                                  ! o2
!             colamt(k,7) = temcol(k)                                  ! co - notused
            enddo
          endif

!> -# Read aerosol optical properties from 'aerosols'.

          do k = 1, nlay
            kk = nlp1 - k
            do ib = 1, nbdsw
              !tauae(k,ib) = aerosols(j1,kk,ib,1)
              !ssaae(k,ib) = aerosols(j1,kk,ib,2)
              !asyae(k,ib) = aerosols(j1,kk,ib,3)
              tauae(k,ib) = aeraod(j1,kk,ib)
              ssaae(k,ib) = aerssa(j1,kk,ib)
              asyae(k,ib) = aerasy(j1,kk,ib)
            enddo
          enddo

!> -# Read cloud optical properties from 'clouds'.
          if (iswcliq > 0) then    ! use prognostic cloud method
            do k = 1, nlay
              kk = nlp1 - k
              cfrac(k) = cld_cf(j1,kk)        ! cloud fraction
              cliqp(k) = cld_lwp(j1,kk)       ! cloud liq path
              reliq(k) = cld_ref_liq(j1,kk)   ! liq partical effctive radius
              cicep(k) = cld_iwp(j1,kk)       ! cloud ice path
              reice(k) = cld_ref_ice(j1,kk)   ! ice partical effctive radius
              cdat1(k) = cld_rwp(j1,kk)       ! cloud rain drop path
              cdat2(k) = cld_ref_rain(j1,kk)  ! rain partical effctive radius
              cdat3(k) = cld_swp(j1,kk)       ! cloud snow path
              cdat4(k) = cld_ref_snow(j1,kk)  ! snow partical effctive radius
            enddo
          else                     ! use diagnostic cloud method
            do k = 1, nlay
              kk = nlp1 - k
              cfrac(k) = cld_cf(j1,kk)      ! cloud fraction
              cdat1(k) = cld_od(j1,kk)      ! cloud optical depth
              cdat2(k) = cld_ssa(j1,kk)      ! cloud single scattering albedo
              cdat3(k) = cld_asy(j1,kk)      ! cloud asymmetry factor
            enddo
          endif                    ! end if_iswcliq

        else                        ! input from sfc to toa

          tem1 = 100.0 * con_g
          tem2 = 1.0e-20 * 1.0e3 * con_avgd

          do k = 1, nlay
            pavel(k) = plyr(j1,k)
            tavel(k) = tlyr(j1,k)
            delp (k) = plvl(j1,k) - plvl(j1,k+1)

!  --- ...  set absorber amount
!test use
!           h2ovmr(k)= max(f_zero,qlyr(j1,k)*amdw)                    ! input mass mixing ratio
!           h2ovmr(k)= max(f_zero,qlyr(j1,k))                         ! input vol mixing ratio
!           o3vmr (k)= max(f_zero,olyr(j1,k))                         ! input vol mixing ratio
!ncep model use
            h2ovmr(k)= max(f_zero,qlyr(j1,k)*amdw/(f_one-qlyr(j1,k))) ! input specific humidity
            o3vmr (k)= max(f_zero,olyr(j1,k)*amdo3)                   ! input mass mixing ratio

            tem0 = (f_one - h2ovmr(k))*con_amd + h2ovmr(k)*con_amw
            coldry(k) = tem2 * delp(k) / (tem1*tem0*(f_one + h2ovmr(k)))
            temcol(k) = 1.0e-12 * coldry(k)

            colamt(k,1) = max(f_zero,    coldry(k)*h2ovmr(k))         ! h2o
            !colamt(k,2) = max(temcol(k), coldry(k)*gasvmr(j1,k,1))    ! co2
            colamt(k,2) = max(temcol(k), coldry(k)*gasvmr_co2(j1,k))  ! co2
            colamt(k,3) = max(f_zero,    coldry(k)*o3vmr(k))          ! o3
            colmol(k)   = coldry(k) + colamt(k,1)
          enddo


          if (lprnt) then
            if (ipt == 1) then
              write(0,*)' pavel=',pavel
              write(0,*)' tavel=',tavel
              write(0,*)' delp=',delp
              write(0,*)' h2ovmr=',h2ovmr*1000
              write(0,*)' o3vmr=',o3vmr*1000000
            endif
          endif

!  --- ...  set up gas column amount, convert from volume mixing ratio
!           to molec/cm2 based on coldry (scaled to 1.0e-20)

          if (iswrgas > 0) then
            do k = 1, nlay
              !colamt(k,4) = max(temcol(k), coldry(k)*gasvmr(j1,k,2))  ! n2o
            colamt(k,4) = max(temcol(k), coldry(k)*gasvmr_n2o(j1,k))  ! n2o
              !colamt(k,5) = max(temcol(k), coldry(k)*gasvmr(j1,k,3))  ! ch4
            colamt(k,5) = max(temcol(k), coldry(k)*gasvmr_ch4(j1,k))  ! ch4
              !colamt(k,6) = max(temcol(k), coldry(k)*gasvmr(j1,k,4))  ! o2
            colamt(k,6) = max(temcol(k), coldry(k)*gasvmr_o2(j1,k))  ! o2
!             colamt(k,7) = max(temcol(k), coldry(k)*gasvmr(j1,k,5))  ! co - notused
            enddo
          else
            do k = 1, nlay
              colamt(k,4) = temcol(k)                                 ! n2o
              colamt(k,5) = temcol(k)                                 ! ch4
              colamt(k,6) = temcol(k)                                 ! o2
!             colamt(k,7) = temcol(k)                                 ! co - notused
            enddo
          endif

!  --- ...  set aerosol optical properties

          do ib = 1, nbdsw
            do k = 1, nlay
              !tauae(k,ib) = aerosols(j1,k,ib,1)
              !ssaae(k,ib) = aerosols(j1,k,ib,2)
              !asyae(k,ib) = aerosols(j1,k,ib,3)
              tauae(k,ib) = aeraod(j1,k,ib)
              ssaae(k,ib) = aerssa(j1,k,ib)
              asyae(k,ib) = aerasy(j1,k,ib)
            enddo
          enddo

          if (iswcliq > 0) then    ! use prognostic cloud method
            do k = 1, nlay
              cfrac(k) = cld_cf(j1,k)       ! cloud fraction
              cliqp(k) = cld_lwp(j1,k)       ! cloud liq path
              reliq(k) = cld_ref_liq(j1,k)   ! liq partical effctive radius
              cicep(k) = cld_iwp(j1,k)       ! cloud ice path
              reice(k) = cld_ref_ice(j1,k)   ! ice partical effctive radius
              cdat1(k) = cld_rwp(j1,k)       ! cloud rain drop path
              cdat2(k) = cld_ref_rain(j1,k)  ! rain partical effctive radius
              cdat3(k) = cld_swp(j1,k)       ! cloud snow path
              cdat4(k) = cld_ref_snow(j1,k)  ! snow partical effctive radius
            enddo
          else                     ! use diagnostic cloud method
            do k = 1, nlay
              cfrac(k) = cld_cf(j1,k)     ! cloud fraction
              cdat1(k) = cld_od(j1,k)     ! cloud optical depth
              cdat2(k) = cld_ssa(j1,k)    ! cloud single scattering albedo
              cdat3(k) = cld_asy(j1,k)    ! cloud asymmetry factor
            enddo
          endif                    ! end if_iswcliq

        endif                       ! if_ivflip

!> -# Compute fractions of clear sky view:
!!    - random overlapping
!!    - max/ran overlapping
!!    - maximum overlapping

        zcf0   = f_one
        zcf1   = f_one
        if (iovrsw == 0) then                    ! random overlapping
          do k = 1, nlay
            zcf0 = zcf0 * (f_one - cfrac(k))
          enddo
        else if (iovrsw == 1) then               ! max/ran overlapping
          do k = 1, nlay
            if (cfrac(k) > ftiny) then                ! cloudy layer
              zcf1 = min ( zcf1, f_one-cfrac(k) )
            elseif (zcf1 < f_one) then                ! clear layer
              zcf0 = zcf0 * zcf1
              zcf1 = f_one
            endif
          enddo
          zcf0 = zcf0 * zcf1
        else if (iovrsw == 2) then               ! maximum overlapping
          do k = 1, nlay
            zcf0 = min ( zcf0, f_one-cfrac(k) )
          enddo
        endif

        if (zcf0 <= ftiny) zcf0 = f_zero
        if (zcf0 > oneminus) zcf0 = f_one
        zcf1 = f_one - zcf0

!> -# For cloudy sky column, call cldprop() to compute the cloud
!!    optical properties for each cloudy layer.

        if (zcf1 > f_zero) then     ! cloudy sky column

          call cldprop                                                  &
!  ---  inputs:
     &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     &
     &       zcf1, nlay, ipseed(j1),                                    &
!  ---  outputs:
     &       taucw, ssacw, asycw, cldfrc, cldfmc                        &
     &     )

        else                        ! clear sky column
          cldfrc(:)  = f_zero
          cldfmc(:,:)= f_zero
          do i = 1, nbdsw
            do k = 1, nlay
              taucw(k,i) = f_zero
              ssacw(k,i) = f_zero
              asycw(k,i) = f_zero
            enddo
          enddo
        endif   ! end if_zcf1_block

!> -# Call setcoef() to compute various coefficients needed in
!!    radiative transfer calculations.
        call setcoef                                                    &
!  ---  inputs:
     &     ( pavel,tavel,h2ovmr, nlay,nlp1,                             &
!  ---  outputs:
     &       laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,                 &
     &       selffac,selffrac,indself,forfac,forfrac,indfor             &
     &     )

!> -# Call taumol() to calculate optical depths for gaseous absorption
!!    and rayleigh scattering
        call taumol                                                     &
!  ---  inputs:
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   &
     &       forfac,forfrac,indfor,selffac,selffrac,indself, NLAY,      &
!  ---  outputs:
     &       sfluxzen, taug, taur                                       &
     &     )

!> -# Call the 2-stream radiation transfer model:
!!    - if physparam::isubcsw .le.0, using standard cloud scheme,
!!      call spcvrtc().
!!    - if physparam::isubcsw .gt.0, using mcica cloud scheme,
!!      call spcvrtm().

        if ( isubcsw <= 0 ) then     ! use standard cloud scheme

          call spcvrtc                                                  &
!  ---  inputs:
     &     ( ssolar,cosz1,sntz1,albbm,albdf,sfluxzen,cldfrc,            &
     &       zcf1,zcf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,   &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       fxupc,fxdnc,fxup0,fxdn0,                                   &
     &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
     &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
     &     )

        else                         ! use mcica cloud scheme

          call spcvrtm                                                  &
!  ---  inputs:
     &     ( ssolar,cosz1,sntz1,albbm,albdf,sfluxzen,cldfmc,            &
     &       zcf1,zcf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,   &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       fxupc,fxdnc,fxup0,fxdn0,                                   &
     &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
     &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
     &     )

        endif

!> -# Save outputs.
!  --- ...  sum up total spectral fluxes for total-sky

        do k = 1, nlp1
          flxuc(k) = f_zero
          flxdc(k) = f_zero

          do ib = 1, nbdsw
            flxuc(k) = flxuc(k) + fxupc(k,ib)
            flxdc(k) = flxdc(k) + fxdnc(k,ib)
          enddo
        enddo

!! --- ...  optional clear sky fluxes

        if ( lhsw0 .or. lflxprf ) then
          do k = 1, nlp1
            flxu0(k) = f_zero
            flxd0(k) = f_zero

            do ib = 1, nbdsw
              flxu0(k) = flxu0(k) + fxup0(k,ib)
              flxd0(k) = flxd0(k) + fxdn0(k,ib)
            enddo
          enddo
        endif

!  --- ...  prepare for final outputs

        do k = 1, nlay
          rfdelp(k) = heatfac / delp(k)
        enddo

        if ( lfdncmp ) then
!! --- ...  optional uv-b surface downward flux
          fdncmp(j1)%uvbf0 = suvbf0
          fdncmp(j1)%uvbfc = suvbfc

!! --- ...  optional beam and diffuse sfc fluxes
          fdncmp(j1)%nirbm = sfbmc(1)
          fdncmp(j1)%nirdf = sfdfc(1)
          fdncmp(j1)%visbm = sfbmc(2)
          fdncmp(j1)%visdf = sfdfc(2)
        endif    ! end if_lfdncmp

!  --- ...  toa and sfc fluxes

        topflx(j1)%upfxc = ftoauc
        topflx(j1)%dnfxc = ftoadc
        topflx(j1)%upfx0 = ftoau0

        sfcflx(j1)%upfxc = fsfcuc
        sfcflx(j1)%dnfxc = fsfcdc
        sfcflx(j1)%upfx0 = fsfcu0
        sfcflx(j1)%dnfx0 = fsfcd0

        if (ivflip == 0) then       ! output from toa to sfc

!  --- ...  compute heating rates

          fnet(1) = flxdc(1) - flxuc(1)

          do k = 2, nlp1
            kk = nlp1 - k + 1
            fnet(k) = flxdc(k) - flxuc(k)
            hswc(j1,kk) = (fnet(k)-fnet(k-1)) * rfdelp(k-1)
          enddo

!! --- ...  optional flux profiles

          if ( lflxprf ) then
            do k = 1, nlp1
              kk = nlp1 - k + 1
              flxprf(j1,kk)%upfxc = flxuc(k)
              flxprf(j1,kk)%dnfxc = flxdc(k)
              flxprf(j1,kk)%upfx0 = flxu0(k)
              flxprf(j1,kk)%dnfx0 = flxd0(k)
            enddo
          endif

!! --- ...  optional clear sky heating rates

          if ( lhsw0 ) then
            fnet(1) = flxd0(1) - flxu0(1)

            do k = 2, nlp1
              kk = nlp1 - k + 1
              fnet(k) = flxd0(k) - flxu0(k)
              hsw0(j1,kk) = (fnet(k)-fnet(k-1)) * rfdelp(k-1)
            enddo
          endif

!! --- ...  optional spectral band heating rates

          if ( lhswb ) then
            do mb = 1, nbdsw
              fnet(1) = fxdnc(1,mb) - fxupc(1,mb)

              do k = 2, nlp1
                kk = nlp1 - k + 1
                fnet(k) = fxdnc(k,mb) - fxupc(k,mb)
                hswb(j1,kk,mb) = (fnet(k) - fnet(k-1)) * rfdelp(k-1)
              enddo
            enddo
          endif

        else                        ! output from sfc to toa

!  --- ...  compute heating rates

          fnet(1) = flxdc(1) - flxuc(1)

          do k = 2, nlp1
            fnet(k) = flxdc(k) - flxuc(k)
            hswc(j1,k-1) = (fnet(k)-fnet(k-1)) * rfdelp(k-1)
          enddo

!! --- ...  optional flux profiles

          if ( lflxprf ) then
            do k = 1, nlp1
              flxprf(j1,k)%upfxc = flxuc(k)
              flxprf(j1,k)%dnfxc = flxdc(k)
              flxprf(j1,k)%upfx0 = flxu0(k)
              flxprf(j1,k)%dnfx0 = flxd0(k)
            enddo
          endif

!! --- ...  optional clear sky heating rates

          if ( lhsw0 ) then
            fnet(1) = flxd0(1) - flxu0(1)

            do k = 2, nlp1
              fnet(k) = flxd0(k) - flxu0(k)
              hsw0(j1,k-1) = (fnet(k)-fnet(k-1)) * rfdelp(k-1)
            enddo
          endif

!! --- ...  optional spectral band heating rates

          if ( lhswb ) then
            do mb = 1, nbdsw
              fnet(1) = fxdnc(1,mb) - fxupc(1,mb)

              do k = 1, nlay
                fnet(k+1) = fxdnc(k+1,mb) - fxupc(k+1,mb)
                hswb(j1,k,mb) = (fnet(k+1) - fnet(k)) * rfdelp(k)
              enddo
            enddo
          endif

        endif                       ! if_ivflip

      enddo   lab_do_ipt

      return
!...................................
      end subroutine rrtmg_sw_run
!-----------------------------------
!> @}

      subroutine rrtmg_sw_finalize ()
      end subroutine rrtmg_sw_finalize


!>\ingroup module_radsw_main
!> This subroutine initializes non-varying module variables, conversion
!! factors, and look-up tables.
!!\param me             print control for parallel process
!>\section rswinit_gen rswinit General Algorithm
!! @{
!-----------------------------------
      subroutine rswinit                                                &
     &     ( me ) !  ---  inputs:
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
!  external module variables:  (in physparam)                           !
!   iswrate - heating rate unit selections                              !
!           =1: output in k/day                                         !
!           =2: output in k/second                                      !
!   iswrgas - control flag for rare gases (ch4,n2o,o2, etc.)            !
!           =0: do not include rare gases                               !
!           >0: include all rare gases                                  !
!   iswcliq - liquid cloud optical properties contrl flag               !
!           =0: input cloud opt depth from diagnostic scheme            !
!           >0: input cwp,rew, and other cloud content parameters       !
!   isubcsw - sub-column cloud approximation control flag               !
!           =0: no sub-col cld treatment, use grid-mean cld quantities  !
!           =1: mcica sub-col, prescribed seeds to get random numbers   !
!           =2: mcica sub-col, providing array icseed for random numbers!
!   icldflg - cloud scheme control flag                                 !
!           =0: diagnostic scheme gives cloud tau, omiga, and g.        !
!           =1: prognostic scheme gives cloud liq/ice path, etc.        !
!   iovrsw  - clouds vertical overlapping control flag                  !
!           =0: random overlapping clouds                               !
!           =1: maximum/random overlapping clouds                       !
!           =2: maximum overlap cloud                                   !
!   iswmode - control flag for 2-stream transfer scheme                 !
!           =1; delta-eddington    (joseph et al., 1976)                !
!           =2: pifm               (zdunkowski et al., 1980)            !
!           =3: discrete ordinates (liou, 1973)                         !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
! definitions:                                                          !
!     arrays for 10000-point look-up tables:                            !
!     tau_tbl  clear-sky optical depth                                  !
!     exp_tbl  exponential lookup table for transmittance               !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: me

!  ---  outputs: none

!  ---  locals:
      real (kind=kind_phys), parameter :: expeps = 1.e-20

      integer :: i

      real (kind=kind_phys) :: tfn, tau

!
!===> ... begin here
!
      if ( iovrsw<0 .or. iovrsw>2 ) then
        print *,'  *** Error in specification of cloud overlap flag',   &
     &          ' IOVRSW=',iovrsw,' in RSWINIT !!'
        stop
      endif

      if (me == 0) then
        print *,' - Using AER Shortwave Radiation, Version: ',VTAGSW

        if (iswmode == 1) then
          print *,'   --- Delta-eddington 2-stream transfer scheme'
        else if (iswmode == 2) then
          print *,'   --- PIFM 2-stream transfer scheme'
        else if (iswmode == 3) then
          print *,'   --- Discrete ordinates 2-stream transfer scheme'
        endif

        if (iswrgas <= 0) then
          print *,'   --- Rare gases absorption is NOT included in SW'
        else
          print *,'   --- Include rare gases N2O, CH4, O2, absorptions',&
     &            ' in SW'
        endif

        if ( isubcsw == 0 ) then
          print *,'   --- Using standard grid average clouds, no ',     &
     &            'sub-column clouds approximation applied'
        elseif ( isubcsw == 1 ) then
          print *,'   --- Using MCICA sub-colum clouds approximation ', &
     &            'with a prescribed sequence of permutation seeds'
        elseif ( isubcsw == 2 ) then
          print *,'   --- Using MCICA sub-colum clouds approximation ', &
     &            'with provided input array of permutation seeds'
        else
          print *,'  *** Error in specification of sub-column cloud ',  &
     &            ' control flag isubcsw =',isubcsw,' !!'
          stop
        endif
      endif

!> -# Check cloud flags for consistency.

      if ((icldflg == 0 .and. iswcliq /= 0) .or.                        &
     &    (icldflg == 1 .and. iswcliq == 0)) then
        print *,'  *** Model cloud scheme inconsistent with SW',        &
     &          ' radiation cloud radiative property setup !!'
        stop
      endif

!> -# Setup constant factors for heating rate
!! the 1.0e-2 is to convert pressure from mb to \f$N/m^2\f$ .

      if (iswrate == 1) then
!       heatfac = 8.4391
!       heatfac = con_g * 86400. * 1.0e-2 / con_cp  !   (in k/day)
        heatfac = con_g * 864.0 / con_cp            !   (in k/day)
      else
        heatfac = con_g * 1.0e-2 / con_cp           !   (in k/second)
      endif

!> -# Define exponential lookup tables for transmittance. 
!          tau is  computed as a function of the \a tau transition function, and
!           transmittance is calculated as a function of tau.  all tables
!           are computed at intervals of 0.0001.  the inverse of the
!           constant used in the Pade approximation to the tau transition
!           function is set to bpade.

      exp_tbl(0) = 1.0
      exp_tbl(NTBMX) = expeps

      do i = 1, NTBMX-1
        tfn = float(i) / float(NTBMX-i)
        tau = bpade * tfn
        exp_tbl(i) = exp( -tau )
      enddo

      return
!...................................
      end subroutine rswinit
!! @}
!-----------------------------------

!>\ingroup module_radsw_main
!> This subroutine computes the cloud optical properties for each
!! cloudy layer and g-point interval.
!!\param cfrac          layer cloud fraction
!!\n for  physparam::iswcliq > 0 (prognostic cloud scheme)  - - -
!!\param cliqp          layer in-cloud liq water path (\f$g/m^2\f$)
!!\param reliq          mean eff radius for liq cloud (micron)
!!\param cicep          layer in-cloud ice water path (\f$g/m^2\f$)
!!\param reice          mean eff radius for ice cloud (micron)
!!\param cdat1          layer rain drop water path (\f$g/m^2\f$)
!!\param cdat2          effective radius for rain drop (micron)
!!\param cdat3          layer snow flake water path(\f$g/m^2\f$)
!!\param cdat4          mean eff radius for snow flake(micron)
!!\n for physparam::iswcliq = 0  (diagnostic cloud scheme)  - - -
!!\param cliqp          not used
!!\param cicep          not used
!!\param reliq          not used
!!\param reice          not used
!!\param cdat1          layer cloud optical depth
!!\param cdat2          layer cloud single scattering albedo
!!\param cdat3          layer cloud asymmetry factor
!!\param cdat4          optional use
!!\param cf1            effective total cloud cover at surface
!!\param nlay           vertical layer number
!!\param ipseed         permutation seed for generating random numbers
!!                      (isubcsw>0)
!!\param taucw          cloud optical depth, w/o delta scaled
!!\param ssacw          weighted cloud single scattering albedo
!!                      (ssa = ssacw / taucw)
!!\param asycw          weighted cloud asymmetry factor
!!                      (asy = asycw / ssacw)
!!\param cldfrc         cloud fraction of grid mean value
!!\param cldfmc         cloud fraction for each sub-column
!!\section General_cldprop cldprop General Algorithm
!> @{
!-----------------------------------
      subroutine cldprop                                                &
     &     ( cfrac,cliqp,reliq,cicep,reice,cdat1,cdat2,cdat3,cdat4,     &   !  ---  inputs
     &       cf1, nlay, ipseed,                                         &
     &       taucw, ssacw, asycw, cldfrc, cldfmc                        &   !  ---  output
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
! Purpose: Compute the cloud optical properties for each cloudy layer   !
! and g-point interval.                                                 !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!    cfrac - real, layer cloud fraction                            nlay !
!        .....  for  iswcliq > 0 (prognostic cloud sckeme)  - - -       !
!    cliqp - real, layer in-cloud liq water path (g/m**2)          nlay !
!    reliq - real, mean eff radius for liq cloud (micron)          nlay !
!    cicep - real, layer in-cloud ice water path (g/m**2)          nlay !
!    reice - real, mean eff radius for ice cloud (micron)          nlay !
!    cdat1 - real, layer rain drop water path (g/m**2)             nlay !
!    cdat2 - real, effective radius for rain drop (micron)         nlay !
!    cdat3 - real, layer snow flake water path(g/m**2)             nlay !
!    cdat4 - real, mean eff radius for snow flake(micron)          nlay !
!        .....  for iswcliq = 0  (diagnostic cloud sckeme)  - - -       !
!    cdat1 - real, layer cloud optical depth                       nlay !
!    cdat2 - real, layer cloud single scattering albedo            nlay !
!    cdat3 - real, layer cloud asymmetry factor                    nlay !
!    cdat4 - real, optional use                                    nlay !
!    cliqp - real, not used                                        nlay !
!    cicep - real, not used                                        nlay !
!    reliq - real, not used                                        nlay !
!    reice - real, not used                                        nlay !
!                                                                       !
!    cf1   - real, effective total cloud cover at surface           1   !
!    nlay  - integer, vertical layer number                         1   !
!    ipseed- permutation seed for generating random numbers (isubcsw>0) !
!                                                                       !
!  outputs:                                                             !
!    taucw  - real, cloud optical depth, w/o delta scaled    nlay*nbdsw !
!    ssacw  - real, weighted cloud single scattering albedo  nlay*nbdsw !
!                             (ssa = ssacw / taucw)                     !
!    asycw  - real, weighted cloud asymmetry factor          nlay*nbdsw !
!                             (asy = asycw / ssacw)                     !
!    cldfrc - real, cloud fraction of grid mean value              nlay !
!    cldfmc - real, cloud fraction for each sub-column       nlay*ngptsw!
!                                                                       !
!                                                                       !
!  explanation of the method for each value of iswcliq, and iswcice.    !
!  set up in module "physparam"                                         !
!                                                                       !
!     iswcliq=0  : input cloud optical property (tau, ssa, asy).        !
!                  (used for diagnostic cloud method)                   !
!     iswcliq>0  : input cloud liq/ice path and effective radius, also  !
!                  require the user of 'iswcice' to specify the method  !
!                  used to compute aborption due to water/ice parts.    !
!  ...................................................................  !
!                                                                       !
!     iswcliq=1  : liquid water cloud optical properties are computed   !
!                  as in hu and stamnes (1993), j. clim., 6, 728-742.   !
!                                                                       !
!     iswcice used only when iswcliq > 0                                !
!                  the cloud ice path (g/m2) and ice effective radius   !
!                  (microns) are inputs.                                !
!     iswcice=1  : ice cloud optical properties are computed as in      !
!                  ebert and curry (1992), jgr, 97, 3831-3836.          !
!     iswcice=2  : ice cloud optical properties are computed as in      !
!                  streamer v3.0 (2001), key, streamer user's guide,    !
!                  cooperative institude for meteorological studies,95pp!
!     iswcice=3  : ice cloud optical properties are computed as in      !
!                  fu (1996), j. clim., 9.                              !
!                                                                       !
!  other cloud control module variables:                                !
!     isubcsw =0: standard cloud scheme, no sub-col cloud approximation !
!             >0: mcica sub-col cloud scheme using ipseed as permutation!
!                 seed for generating rundom numbers                    !
!                                                                       !
!  ======================  end of description block  =================  !
!
      use module_radsw_cldprtb

!  ---  inputs:
      integer, intent(in) :: nlay, ipseed
      real (kind=kind_phys), intent(in) :: cf1

      real (kind=kind_phys), dimension(nlay), intent(in) :: cliqp,      &
     &       reliq, cicep, reice, cdat1, cdat2, cdat3, cdat4, cfrac

!  ---  outputs:
      real (kind=kind_phys), dimension(nlay,ngptsw), intent(out) ::     &
     &       cldfmc
      real (kind=kind_phys), dimension(nlay,nbdsw),  intent(out) ::     &
     &       taucw, ssacw, asycw
      real (kind=kind_phys), dimension(nlay), intent(out) :: cldfrc

!  ---  locals:
      real (kind=kind_phys), dimension(nblow:nbhgh) :: tauliq, tauice,  &
     &       ssaliq, ssaice, ssaran, ssasnw, asyliq, asyice,            &
     &       asyran, asysnw
      real (kind=kind_phys), dimension(nlay)       :: cldf

      real (kind=kind_phys) :: dgeice, factor, fint, tauran, tausnw,    &
     &       cldliq, refliq, cldice, refice, cldran, cldsnw, refsnw,    &
     &       extcoliq, ssacoliq, asycoliq, extcoice, ssacoice, asycoice,&
     &       dgesnw

      logical :: lcloudy(nlay,ngptsw)
      integer :: ia, ib, ig, jb, k, index

!
!===> ...  begin here
!
      do ib = 1, nbdsw
        do k = 1, nlay
          taucw (k,ib) = f_zero
          ssacw (k,ib) = f_one
          asycw (k,ib) = f_zero
        enddo
      enddo

!> -# Compute cloud radiative properties for a cloudy column.

      lab_if_iswcliq : if (iswcliq > 0) then

        lab_do_k : do k = 1, nlay
          lab_if_cld : if (cfrac(k) > ftiny) then

!>  - Compute optical properties for rain and snow.
!!\n    For rain: tauran/ssaran/asyran
!!\n    For snow: tausnw/ssasnw/asysnw
!>  - Calculation of absorption coefficients due to water clouds
!!\n    For water clouds: tauliq/ssaliq/asyliq
!>  - Calculation of absorption coefficients due to ice clouds
!!\n    For ice clouds: tauice/ssaice/asyice
!>  - For Prognostic cloud scheme: sum up the cloud optical property:
!!\n     \f$ taucw=tauliq+tauice+tauran+tausnw \f$
!!\n     \f$ ssacw=ssaliq+ssaice+ssaran+ssasnw \f$
!!\n     \f$ asycw=asyliq+asyice+asyran+asysnw \f$

            cldran = cdat1(k)
!           refran = cdat2(k)
            cldsnw = cdat3(k)
            refsnw = cdat4(k)
            dgesnw = 1.0315 * refsnw        ! for fu's snow formula

            tauran = cldran * a0r

!>  - If use fu's formula it needs to be normalized by snow/ice density.
!! not use snow density = 0.1 g/cm**3 = 0.1 g/(mu * m**2)
!!\n use ice density = 0.9167 g/cm**3 = 0.9167 g/(mu * m**2)
!!\n       1/0.9167 = 1.09087
!!\n       factor 1.5396=8/(3*sqrt(3)) converts reff to generalized ice particle size
!!       use newer factor value 1.0315
            if (cldsnw>f_zero .and. refsnw>10.0_kind_phys) then
!             tausnw = cldsnw * (a0s + a1s/refsnw)
              tausnw = cldsnw*1.09087*(a0s + a1s/dgesnw)     ! fu's formula
            else
              tausnw = f_zero
            endif

            do ib = nblow, nbhgh
              ssaran(ib) = tauran * (f_one - b0r(ib))
              ssasnw(ib) = tausnw * (f_one - (b0s(ib)+b1s(ib)*dgesnw))
              asyran(ib) = ssaran(ib) * c0r(ib)
              asysnw(ib) = ssasnw(ib) * c0s(ib)
            enddo

            cldliq = cliqp(k)
            cldice = cicep(k)
            refliq = reliq(k)
            refice = reice(k)

!>  - Calculation of absorption coefficients due to water clouds.

            if ( cldliq <= f_zero ) then
              do ib = nblow, nbhgh
                tauliq(ib) = f_zero
                ssaliq(ib) = f_zero
                asyliq(ib) = f_zero
              enddo
            else
              if ( iswcliq == 1 ) then
                factor = refliq - 1.5
                index  = max( 1, min( 57, int( factor ) ))
                fint   = factor - float(index)

                do ib = nblow, nbhgh
                  extcoliq = max(f_zero,            extliq1(index,ib)   &
     &              + fint*(extliq1(index+1,ib)-extliq1(index,ib)) )
                  ssacoliq = max(f_zero, min(f_one, ssaliq1(index,ib)   &
     &              + fint*(ssaliq1(index+1,ib)-ssaliq1(index,ib)) ))

                  asycoliq = max(f_zero, min(f_one, asyliq1(index,ib)   &
     &              + fint*(asyliq1(index+1,ib)-asyliq1(index,ib)) ))
!                 forcoliq = asycoliq * asycoliq

                  tauliq(ib) = cldliq     * extcoliq
                  ssaliq(ib) = tauliq(ib) * ssacoliq
                  asyliq(ib) = ssaliq(ib) * asycoliq
                enddo
              endif   ! end if_iswcliq_block
            endif   ! end if_cldliq_block

!>  - Calculation of absorption coefficients due to ice clouds.

            if ( cldice <= f_zero ) then
              do ib = nblow, nbhgh
                tauice(ib) = f_zero
                ssaice(ib) = f_zero
                asyice(ib) = f_zero
              enddo
            else

!>   - ebert and curry approach for all particle sizes though somewhat
!! unjustified for large ice particles.

              if ( iswcice == 1 ) then
                refice = min(130.0_kind_phys,max(13.0_kind_phys,refice))

                do ib = nblow, nbhgh
                  ia = idxebc(ib)           ! eb_&_c band index for ice cloud coeff

                  extcoice = max(f_zero, abari(ia)+bbari(ia)/refice )
                  ssacoice = max(f_zero, min(f_one,                     &
     &                             f_one-cbari(ia)-dbari(ia)*refice ))
                  asycoice = max(f_zero, min(f_one,                     &
     &                                   ebari(ia)+fbari(ia)*refice ))
!                 forcoice = asycoice * asycoice

                  tauice(ib) = cldice     * extcoice
                  ssaice(ib) = tauice(ib) * ssacoice
                  asyice(ib) = ssaice(ib) * asycoice
                enddo

!>   - streamer approach for ice effective radius between 5.0 and 131.0 microns.

              elseif ( iswcice == 2 ) then
                refice = min(131.0_kind_phys,max(5.0_kind_phys,refice))

                factor = (refice - 2.0) / 3.0
                index  = max( 1, min( 42, int( factor ) ))
                fint   = factor - float(index)

                do ib = nblow, nbhgh
                  extcoice = max(f_zero,            extice2(index,ib)   &
     &                + fint*(extice2(index+1,ib)-extice2(index,ib)) )
                  ssacoice = max(f_zero, min(f_one, ssaice2(index,ib)   &
     &                + fint*(ssaice2(index+1,ib)-ssaice2(index,ib)) ))
                  asycoice = max(f_zero, min(f_one, asyice2(index,ib)   &
     &                + fint*(asyice2(index+1,ib)-asyice2(index,ib)) ))
!                 forcoice = asycoice * asycoice

                  tauice(ib) = cldice     * extcoice
                  ssaice(ib) = tauice(ib) * ssacoice
                  asyice(ib) = ssaice(ib) * asycoice
                enddo

!>   - fu's approach for ice effective radius between 4.8 and 135 microns
!! (generalized effective size from 5 to 140 microns).

              elseif ( iswcice == 3 ) then
                dgeice = max( 5.0, min( 140.0, 1.0315*refice ))

                factor = (dgeice - 2.0) / 3.0
                index  = max( 1, min( 45, int( factor ) ))
                fint   = factor - float(index)

                do ib = nblow, nbhgh
                  extcoice = max(f_zero,            extice3(index,ib)   &
     &                + fint*(extice3(index+1,ib)-extice3(index,ib)) )
                  ssacoice = max(f_zero, min(f_one, ssaice3(index,ib)   &
     &                + fint*(ssaice3(index+1,ib)-ssaice3(index,ib)) ))
                  asycoice = max(f_zero, min(f_one, asyice3(index,ib)   &
     &                + fint*(asyice3(index+1,ib)-asyice3(index,ib)) ))
!                 fdelta   = max(f_zero, min(f_one, fdlice3(index,ib)   &
!    &                + fint*(fdlice3(index+1,ib)-fdlice3(index,ib)) ))
!                 forcoice = min( asycoice, fdelta+0.5/ssacoice )           ! see fu 1996 p. 2067

                  tauice(ib) = cldice     * extcoice
                  ssaice(ib) = tauice(ib) * ssacoice
                  asyice(ib) = ssaice(ib) * asycoice
                enddo

              endif   ! end if_iswcice_block
            endif   ! end if_cldice_block

            do ib = 1, nbdsw
              jb = nblow + ib - 1
              taucw(k,ib) = tauliq(jb)+tauice(jb)+tauran+tausnw
              ssacw(k,ib) = ssaliq(jb)+ssaice(jb)+ssaran(jb)+ssasnw(jb)
              asycw(k,ib) = asyliq(jb)+asyice(jb)+asyran(jb)+asysnw(jb)
            enddo

          endif  lab_if_cld
        enddo  lab_do_k

      else  lab_if_iswcliq

        do k = 1, nlay
          if (cfrac(k) > ftiny) then
            do ib = 1, nbdsw
              taucw(k,ib) = cdat1(k)
              ssacw(k,ib) = cdat1(k)    * cdat2(k)
              asycw(k,ib) = ssacw(k,ib) * cdat3(k)
            enddo
          endif
        enddo

      endif  lab_if_iswcliq

!> -# if physparam::isubcsw > 0, call mcica_subcol() to distribute
!!    cloud properties to each g-point.

      if ( isubcsw > 0 ) then      ! mcica sub-col clouds approx

        cldf(:) = cfrac(:)
        where (cldf(:) < ftiny)
          cldf(:) = f_zero
        end where

!  --- ...  call sub-column cloud generator

        call mcica_subcol                                               &
!  ---  inputs:
     &     ( cldf, nlay, ipseed,                                        &
!  ---  outputs:
     &       lcloudy                                                    &
     &     )

        do ig = 1, ngptsw
          do k = 1, nlay
            if ( lcloudy(k,ig) ) then
              cldfmc(k,ig) = f_one
            else
              cldfmc(k,ig) = f_zero
            endif
          enddo
        enddo

      else                         ! non-mcica, normalize cloud

        do k = 1, nlay
          cldfrc(k) = cfrac(k) / cf1
        enddo
      endif   ! end if_isubcsw_block

      return
!...................................
      end subroutine cldprop
!-----------------------------------
!> @}

!>\ingroup module_radsw_main
!> This subroutine computes the sub-colum cloud profile flag array.
!!\param cldf        layer cloud fraction
!!\param nlay        number of model vertical layers
!!\param ipseed      permute seed for random num generator
!!\param lcloudy     sub-colum cloud profile flag array
!!\section mcica_subcol_gen mcica_subcol General Algorithm
!! @{
! ----------------------------------
      subroutine mcica_subcol                                           &
     &    ( cldf, nlay, ipseed,                                         &       !  ---  inputs
     &      lcloudy                                                     &       !  ---  outputs
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
!   lcloudy - logical, sub-colum cloud profile flag array    nlay*ngptsw!
!                                                                       !
!  other control flags from module variables:                           !
!     iovrsw    : control flag for cloud overlapping method             !
!                 =0:random; =1:maximum/random; =2:maximum              !
!                                                                       !
!                                                                       !
!  =====================    end of definitions    ====================  !

      implicit none

!  ---  inputs:
      integer, intent(in) :: nlay, ipseed

      real (kind=kind_phys), dimension(nlay), intent(in) :: cldf

!  ---  outputs:
      logical, dimension(nlay,ngptsw), intent(out):: lcloudy

!  ---  locals:
      real (kind=kind_phys) :: cdfunc(nlay,ngptsw), tem1,               &
     &       rand2d(nlay*ngptsw), rand1d(ngptsw)

      type (random_stat) :: stat          ! for thread safe random generator

      integer :: k, n, k1
!
!===> ...  begin here
!
!> -# Advance randum number generator by ipseed values.

      call random_setseed                                               &
!  ---  inputs:
     &    ( ipseed,                                                     &
!  ---  outputs:
     &      stat                                                        &
     &    )

!> -# Sub-column set up according to overlapping assumption.

      select case ( iovrsw )

        case( 0 )        ! random overlap, pick a random value at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptsw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(k,n) = rand2d(k1)
            enddo
          enddo

        case( 1 )        ! max-ran overlap

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand2d, stat )

          k1 = 0
          do n = 1, ngptsw
            do k = 1, nlay
              k1 = k1 + 1
              cdfunc(k,n) = rand2d(k1)
            enddo
          enddo

!  ---  first pick a random number for bottom/top layer.
!       then walk up the column: (aer's code)
!       if layer below is cloudy, use the same rand num in the layer below
!       if layer below is clear,  use a new random number

!  ---  from bottom up
          do k = 2, nlay
            k1 = k - 1
            tem1 = f_one - cldf(k1)

            do n = 1, ngptsw
              if ( cdfunc(k1,n) > tem1 ) then
                cdfunc(k,n) = cdfunc(k1,n)
              else
                cdfunc(k,n) = cdfunc(k,n) * tem1
              endif
            enddo
          enddo

!  ---  then walk down the column: (if use original author's method)
!       if layer above is cloudy, use the same rand num in the layer above
!       if layer above is clear,  use a new random number

!  ---  from top down
!         do k = nlay-1, 1, -1
!           k1 = k + 1
!           tem1 = f_one - cldf(k1)

!           do n = 1, ngptsw
!             if ( cdfunc(k1,n) > tem1 ) then
!               cdfunc(k,n) = cdfunc(k1,n)
!             else
!               cdfunc(k,n) = cdfunc(k,n) * tem1
!             endif
!           enddo
!         enddo

        case( 2 )        ! maximum overlap, pick same random numebr at every level

          call random_number                                            &
!  ---  inputs: ( none )
!  ---  outputs:
     &     ( rand1d, stat )

          do n = 1, ngptsw
            tem1 = rand1d(n)

            do k = 1, nlay
              cdfunc(k,n) = tem1
            enddo
          enddo

      end select

!> -# Generate subcolumns for homogeneous clouds.

      do k = 1, nlay
        tem1 = f_one - cldf(k)

        do n = 1, ngptsw
          lcloudy(k,n) = cdfunc(k,n) >= tem1
        enddo
      enddo

      return
! ..................................
      end subroutine mcica_subcol
!! @}
! ----------------------------------

!>\ingroup module_radsw_main
!> This subroutine computes various coefficients needed in radiative
!! transfer calculation.
!!\param pavel           layer pressure (mb)
!!\param tavel           layer temperature (k)
!!\param h2ovmr          layer w.v. volumn mixing ratio (kg/kg)
!!\param nlay            total number of vertical layers
!!\param nlp1            total number of vertical levels
!!\param laytrop         tropopause layer index (unitless)
!!\param jp              indices of lower reference pressure
!!\param jt,jt1          indices of lower reference temperatures at
!!                       levels of jp and jp+1
!!\param facij           factors mltiply the reference ks,i,j=0/1 for
!!                       lower/higher of the 2 appropriate temperature
!!                       and altitudes.
!!\param selffac         scale factor for w. v. self-continuum equals
!!                       (w.v. density)/(atmospheric density at 296k
!!                       and 1013 mb)
!!\param seffrac         factor for temperature interpolation of
!!                       reference w.v. self-continuum data
!!\param indself         index of lower ref temp for selffac
!!\param forfac          scale factor for w. v. foreign-continuum
!!\param forfrac         factor for temperature interpolation of
!!                       reference w.v. foreign-continuum data
!!\param indfor          index of lower ref temp for forfac
!>\section setcoef_gen_rw setcoef General Algorithm
!! @{
! ----------------------------------
      subroutine setcoef                                                &
     &     ( pavel,tavel,h2ovmr, nlay,nlp1,                             &    !  ---  inputs
     &       laytrop,jp,jt,jt1,fac00,fac01,fac10,fac11,                 &    !  ---  outputs
     &       selffac,selffrac,indself,forfac,forfrac,indfor             &
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
!   h2ovmr    - real, layer w.v. volum mixing ratio (kg/kg)        nlay !
!   nlay/nlp1 - integer, total number of vertical layers, levels    1   !
!                                                                       !
!  outputs:                                                             !
!   laytrop   - integer, tropopause layer index (unitless)          1   !
!   jp        - real, indices of lower reference pressure          nlay !
!   jt, jt1   - real, indices of lower reference temperatures      nlay !
!                 at levels of jp and jp+1                              !
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
!                                                                       !
!  ======================    end of definitions    ===================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(:), intent(in) :: pavel, tavel,  &
     &       h2ovmr

!  ---  outputs:
      integer, dimension(nlay), intent(out) :: indself, indfor,         &
     &       jp, jt, jt1
      integer, intent(out) :: laytrop

      real (kind=kind_phys), dimension(nlay), intent(out) :: fac00,     &
     &       fac01, fac10, fac11, selffac, selffrac, forfac, forfrac

!  ---  locals:
      real (kind=kind_phys) :: plog, fp, fp1, ft, ft1, tem1, tem2

      integer :: i, k, jp1
!
!===> ... begin here
!
      laytrop= nlay

      do k = 1, nlay

        forfac(k) = pavel(k)*stpfac / (tavel(k)*(f_one + h2ovmr(k)))

!> -# Find the two reference pressures on either side of the
!! layer pressure.  store them in jp and jp1.  store in fp the
!! fraction of the difference (in ln(pressure)) between these
!! two values that the layer pressure lies.

        plog  = log(pavel(k))
        jp(k) = max(1, min(58, int(36.0 - 5.0*(plog+0.04)) ))
        jp1   = jp(k) + 1
        fp    = 5.0 * (preflog(jp(k)) - plog)

!> -# Determine, for each reference pressure (jp and jp1), which
!! reference temperature (these are different for each reference
!! pressure) is nearest the layer temperature but does not exceed it.
!! store these indices in jt and jt1, resp. store in ft (resp. ft1)
!! the fraction of the way between jt (jt1) and the next highest
!! reference temperature that the layer temperature falls.

        tem1 = (tavel(k) - tref(jp(k))) / 15.0
        tem2 = (tavel(k) - tref(jp1  )) / 15.0
        jt (k) = max(1, min(4, int(3.0 + tem1) ))
        jt1(k) = max(1, min(4, int(3.0 + tem2) ))
        ft  = tem1 - float(jt (k) - 3)
        ft1 = tem2 - float(jt1(k) - 3)

!> -# We have now isolated the layer ln pressure and temperature,
!! between two reference pressures and two reference temperatures
!! (for each reference pressure).  we multiply the pressure
!! fraction fp with the appropriate temperature fractions to get
!! the factors that will be needed for the interpolation that yields
!! the optical depths (performed in routines taugbn for band n).

        fp1 = f_one - fp
        fac10(k) = fp1 * ft
        fac00(k) = fp1 * (f_one - ft)
        fac11(k) = fp  * ft1
        fac01(k) = fp  * (f_one - ft1)

!> -# If the pressure is less than ~100mb, perform a different
!! set of species interpolations.

        if ( plog > 4.56 ) then

          laytrop =  k

!> -# Set up factors needed to separately include the water vapor
!! foreign-continuum in the calculation of absorption coefficient.

          tem1 = (332.0 - tavel(k)) / 36.0
          indfor (k) = min(2, max(1, int(tem1)))
          forfrac(k) = tem1 - float(indfor(k))

!> -# Set up factors needed to separately include the water vapor
!! self-continuum in the calculation of absorption coefficient.

          tem2 = (tavel(k) - 188.0) / 7.2
          indself (k) = min(9, max(1, int(tem2)-7))
          selffrac(k) = tem2 - float(indself(k) + 7)
          selffac (k) = h2ovmr(k) * forfac(k)

        else

!  --- ...  set up factors needed to separately include the water vapor
!           foreign-continuum in the calculation of absorption coefficient.

          tem1 = (tavel(k) - 188.0) / 36.0
          indfor (k) = 3
          forfrac(k) = tem1 - f_one

          indself (k) = 0
          selffrac(k) = f_zero
          selffac (k) = f_zero

        endif

      enddo    ! end_do_k_loop

      return
! ..................................
      end subroutine setcoef
!! @}
! ----------------------------------

!>\ingroup module_radsw_main
!> This subroutine computes the shortwave radiative fluxes using
!! two-stream method.
!!\param ssolar           incoming solar flux at top
!!\param cosz             cosine solar zenith angle
!!\param sntz             secant solar zenith angle
!!\param albbm            surface albedo for direct beam radiation
!!\param albdf            surface albedo for diffused radiation
!!\param sfluxzen         spectral distribution of incoming solar flux
!!\param cldfrc           layer cloud fraction
!!\param cf1              >0: cloudy sky, otherwise: clear sky
!!\param cf0              =1-cf1
!!\param taug             spectral optical depth for gases
!!\param taur             optical depth for rayleigh scattering
!!\param tauae            aerosols optical depth
!!\param ssaae            aerosols single scattering albedo
!!\param asyae            aerosols asymmetry factor
!!\param taucw            weighted cloud optical depth
!!\param ssacw            weighted cloud single scat albedo
!!\param asycw            weighted cloud asymmetry factor
!!\param nlay,nlp1        number of layers/levels
!!\param fxupc            tot sky upward flux
!!\param fxdnc            tot sky downward flux
!!\param fxup0            clr sky upward flux
!!\param fxdn0            clr sky downward flux
!!\param ftoauc           tot sky toa upwd flux
!!\param ftoau0           clr sky toa upwd flux
!!\param ftoadc           toa downward (incoming) solar flux
!!\param fsfcuc           tot sky sfc upwd flux
!!\param fsfcu0           clr sky sfc upwd flux
!!\param fsfcdc           tot sky sfc dnwd flux
!!\param fsfcd0           clr sky sfc dnwd flux
!!\param sfbmc            tot sky sfc dnwd beam flux (nir/uv+vis)
!!\param sfdfc            tot sky sfc dnwd diff flux (nir/uv+vis)
!!\param sfbm0            clr sky sfc dnwd beam flux (nir/uv+vis)
!!\param sfdf0            clr sky sfc dnwd diff flux (nir/uv+vis)
!!\param suvbfc           tot sky sfc dnwd uv-b flux
!!\param suvbf0           clr sky sfc dnwd uv-b flux
!>\section General_spcvrtc spcvrtc General Algorithm
!! @{
!-----------------------------------
      subroutine spcvrtc                                                &
     &     ( ssolar,cosz,sntz,albbm,albdf,sfluxzen,cldfrc,              &  !  ---  inputs
     &       cf1,cf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,     &
     &       nlay, nlp1,                                                &
     &       fxupc,fxdnc,fxup0,fxdn0,                                   &  !  ---  outputs
     &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
     &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
!   purpose:  computes the shortwave radiative fluxes using two-stream  !
!             method                                                    !
!                                                                       !
!   subprograms called:  vrtqdr                                         !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!    ssolar  - real, incoming solar flux at top                    1    !
!    cosz    - real, cosine solar zenith angle                     1    !
!    sntz    - real, secant solar zenith angle                     1    !
!    albbm   - real, surface albedo for direct beam radiation      2    !
!    albdf   - real, surface albedo for diffused radiation         2    !
!    sfluxzen- real, spectral distribution of incoming solar flux ngptsw!
!    cldfrc  - real, layer cloud fraction                         nlay  !
!    cf1     - real, >0: cloudy sky, otherwise: clear sky          1    !
!    cf0     - real, =1-cf1                                        1    !
!    taug    - real, spectral optical depth for gases        nlay*ngptsw!
!    taur    - real, optical depth for rayleigh scattering   nlay*ngptsw!
!    tauae   - real, aerosols optical depth                  nlay*nbdsw !
!    ssaae   - real, aerosols single scattering albedo       nlay*nbdsw !
!    asyae   - real, aerosols asymmetry factor               nlay*nbdsw !
!    taucw   - real, weighted cloud optical depth            nlay*nbdsw !
!    ssacw   - real, weighted cloud single scat albedo       nlay*nbdsw !
!    asycw   - real, weighted cloud asymmetry factor         nlay*nbdsw !
!    nlay,nlp1 - integer,  number of layers/levels                 1    !
!                                                                       !
!  output variables:                                                    !
!    fxupc   - real, tot sky upward flux                     nlp1*nbdsw !
!    fxdnc   - real, tot sky downward flux                   nlp1*nbdsw !
!    fxup0   - real, clr sky upward flux                     nlp1*nbdsw !
!    fxdn0   - real, clr sky downward flux                   nlp1*nbdsw !
!    ftoauc  - real, tot sky toa upwd flux                         1    !
!    ftoau0  - real, clr sky toa upwd flux                         1    !
!    ftoadc  - real, toa downward (incoming) solar flux            1    !
!    fsfcuc  - real, tot sky sfc upwd flux                         1    !
!    fsfcu0  - real, clr sky sfc upwd flux                         1    !
!    fsfcdc  - real, tot sky sfc dnwd flux                         1    !
!    fsfcd0  - real, clr sky sfc dnwd flux                         1    !
!    sfbmc   - real, tot sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdfc   - real, tot sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    sfbm0   - real, clr sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdf0   - real, clr sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    suvbfc  - real, tot sky sfc dnwd uv-b flux                    1    !
!    suvbf0  - real, clr sky sfc dnwd uv-b flux                    1    !
!                                                                       !
!  internal variables:                                                  !
!    zrefb   - real, direct beam reflectivity for clear/cloudy    nlp1  !
!    zrefd   - real, diffuse reflectivity for clear/cloudy        nlp1  !
!    ztrab   - real, direct beam transmissivity for clear/cloudy  nlp1  !
!    ztrad   - real, diffuse transmissivity for clear/cloudy      nlp1  !
!    zldbt   - real, layer beam transmittance for clear/cloudy    nlp1  !
!    ztdbt   - real, lev total beam transmittance for clr/cld     nlp1  !
!                                                                       !
!  control parameters in module "physparam"                             !
!    iswmode - control flag for 2-stream transfer schemes               !
!              = 1 delta-eddington    (joseph et al., 1976)             !
!              = 2 pifm               (zdunkowski et al., 1980)         !
!              = 3 discrete ordinates (liou, 1973)                      !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  method:                                                              !
!  -------                                                              !
!     standard delta-eddington, p.i.f.m., or d.o.m. layer calculations. !
!     kmodts  = 1 eddington (joseph et al., 1976)                       !
!             = 2 pifm (zdunkowski et al., 1980)                        !
!             = 3 discrete ordinates (liou, 1973)                       !
!                                                                       !
!  modifications:                                                       !
!  --------------                                                       !
!   original: h. barker                                                 !
!   revision: merge with rrtmg_sw: j.-j.morcrette, ecmwf, feb 2003      !
!   revision: add adjustment for earth/sun distance:mjiacono,aer,oct2003!
!   revision: bug fix for use of palbp and palbd: mjiacono, aer, nov2003!
!   revision: bug fix to apply delta scaling to clear sky: aer, dec2004 !
!   revision: code modified so that delta scaling is not done in cloudy !
!             profiles if routine cldprop is used; delta scaling can be !
!             applied by swithcing code below if cldprop is not used to !
!             get cloud properties. aer, jan 2005                       !
!   revision: uniform formatting for rrtmg: mjiacono, aer, jul 2006     !
!   revision: use exponential lookup table for transmittance: mjiacono, !
!             aer, aug 2007                                             !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: zcrit = 0.9999995 ! thresold for conservative scattering
      real (kind=kind_phys), parameter :: zsr3  = sqrt(3.0)
      real (kind=kind_phys), parameter :: od_lo = 0.06
      real (kind=kind_phys), parameter :: eps1  = 1.0e-8

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nlay,ngptsw), intent(in) ::      &
     &       taug, taur
      real (kind=kind_phys), dimension(nlay,nbdsw),  intent(in) ::      &
     &       taucw, ssacw, asycw, tauae, ssaae, asyae

      real (kind=kind_phys), dimension(ngptsw), intent(in) :: sfluxzen
      real (kind=kind_phys), dimension(nlay),   intent(in) :: cldfrc

      real (kind=kind_phys), dimension(2),  intent(in) :: albbm, albdf

      real (kind=kind_phys), intent(in) :: cosz, sntz, cf1, cf0, ssolar

!  ---  outputs:
      real (kind=kind_phys), dimension(nlp1,nbdsw), intent(out) ::      &
     &       fxupc, fxdnc, fxup0, fxdn0

      real (kind=kind_phys), dimension(2), intent(out) :: sfbmc, sfdfc, &
     &       sfbm0, sfdf0

      real (kind=kind_phys), intent(out) :: suvbfc, suvbf0, ftoadc,     &
     &       ftoauc, ftoau0, fsfcuc, fsfcu0, fsfcdc, fsfcd0

!  ---  locals:
      real (kind=kind_phys), dimension(nlay) :: ztaus, zssas, zasys,    &
     &       zldbt0

      real (kind=kind_phys), dimension(nlp1) :: zrefb, zrefd, ztrab,    &
     &       ztrad, ztdbt, zldbt, zfu, zfd

      real (kind=kind_phys) :: ztau1, zssa1, zasy1, ztau0, zssa0,       &
     &       zasy0, zasy3, zssaw, zasyw, zgam1, zgam2, zgam3, zgam4,    &
     &       zc0, zc1, za1, za2, zb1, zb2, zrk, zrk2, zrp, zrp1, zrm1,  &
     &       zrpp, zrkg1, zrkg3, zrkg4, zexp1, zexm1, zexp2, zexm2,     &
     &       zexp3, zexp4, zden1, ze1r45, ftind, zsolar, zrefb1,        &
     &       zrefd1, ztrab1, ztrad1, ztdbt0, zr1, zr2, zr3, zr4, zr5,   &
     &       zt1, zt2, zt3, zf1, zf2, zrpp1

      integer :: ib, ibd, jb, jg, k, kp, itind
!
!===> ...  begin here

!> -# Initialize output fluxes.
      do ib = 1, nbdsw
        do k = 1, nlp1
          fxdnc(k,ib) = f_zero
          fxupc(k,ib) = f_zero
          fxdn0(k,ib) = f_zero
          fxup0(k,ib) = f_zero
        enddo
      enddo

      ftoadc = f_zero
      ftoauc = f_zero
      ftoau0 = f_zero
      fsfcuc = f_zero
      fsfcu0 = f_zero
      fsfcdc = f_zero
      fsfcd0 = f_zero

!! --- ...  uv-b surface downward fluxes
      suvbfc  = f_zero
      suvbf0  = f_zero

!! --- ...  output surface flux components
      sfbmc(1) = f_zero
      sfbmc(2) = f_zero
      sfdfc(1) = f_zero
      sfdfc(2) = f_zero
      sfbm0(1) = f_zero
      sfbm0(2) = f_zero
      sfdf0(1) = f_zero
      sfdf0(2) = f_zero

!> -# Loop over all g-points in each band.

      lab_do_jg : do jg = 1, ngptsw

        jb = NGB(jg)
        ib = jb + 1 - nblow
        ibd = idxsfc(jb)

        zsolar = ssolar * sfluxzen(jg)

!> -# Set up toa direct beam and surface values (beam and diff).

        ztdbt(nlp1) = f_one
        ztdbt0   = f_one

        zldbt(1) = f_zero
        if (ibd /= 0) then
          zrefb(1) = albbm(ibd)
          zrefd(1) = albdf(ibd)
        else
          zrefb(1) = 0.5 * (albbm(1) + albbm(2))
          zrefd(1) = 0.5 * (albdf(1) + albdf(2))
        endif
        ztrab(1) = f_zero
        ztrad(1) = f_zero

!> -# Compute clear-sky optical parameters, layer reflectance and
!!    transmittance.
!    - Set up toa direct beam and surface values (beam and diff).
!    - Delta scaling for clear-sky condition.
!    - General two-stream expressions for physparam::iswmode .
!    - Compute homogeneous reflectance and transmittance for both
!      conservative and non-conservative scattering.
!    - Pre-delta-scaling clear and cloudy direct beam transmittance.
!    - Call swflux() to compute the upward and downward radiation
!      fluxes.

        do k = nlay, 1, -1
          kp = k + 1

          ztau0 = max( ftiny, taur(k,jg)+taug(k,jg)+tauae(k,ib) )
          zssa0 = taur(k,jg) + tauae(k,ib)*ssaae(k,ib)
          zasy0 = asyae(k,ib)*ssaae(k,ib)*tauae(k,ib)
          zssaw = min( oneminus, zssa0 / ztau0 )
          zasyw = zasy0 / max( ftiny, zssa0 )

!>  - Saving clear-sky quantities for later total-sky usage.
          ztaus(k) = ztau0
          zssas(k) = zssa0
          zasys(k) = zasy0

!>  - Delta scaling for clear-sky condition.
          za1 = zasyw * zasyw
          za2 = zssaw * za1

          ztau1 = (f_one - za2) * ztau0
          zssa1 = (zssaw - za2) / (f_one - za2)
!org      zasy1 = (zasyw - za1) / (f_one - za1)   ! this line is replaced by the next
          zasy1 = zasyw / (f_one + zasyw)         ! to reduce truncation error
          zasy3 = 0.75 * zasy1

!>  - Perform general two-stream expressions:
!!\n  control parameters in module "physparam"                             
!!\n    iswmode - control flag for 2-stream transfer schemes               
!!\n              = 1 delta-eddington    (joseph et al., 1976)             
!!\n              = 2 pifm               (zdunkowski et al., 1980)         
!!\n              = 3 discrete ordinates (liou, 1973)                      
          if ( iswmode == 1 ) then
            zgam1 = 1.75 - zssa1 * (f_one + zasy3)
            zgam2 =-0.25 + zssa1 * (f_one - zasy3)
            zgam3 = 0.5  - zasy3 * cosz
          elseif ( iswmode == 2 ) then               ! pifm
            zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
            zgam2 = 0.75* zssa1 * (f_one- zasy1)
            zgam3 = 0.5 - zasy3 * cosz
          elseif ( iswmode == 3 ) then               ! discrete ordinates
            zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
            zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
            zgam3 = (1.0 - zsr3 * zasy1 * cosz) * 0.5
          endif
          zgam4 = f_one - zgam3

!>  - Compute homogeneous reflectance and transmittance for both conservative
!! scattering and non-conservative scattering.

          if ( zssaw >= zcrit ) then    ! for conservative scattering
            za1 = zgam1 * cosz - zgam3
            za2 = zgam1 * ztau1

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

            zb1 = min ( ztau1*sntz , 500.0 )
            if ( zb1 <= od_lo ) then
              zb2 = f_one - zb1 + 0.5*zb1*zb1
            else
              ftind = zb1 / (bpade + zb1)
              itind = ftind*NTBMX + 0.5
              zb2 = exp_tbl(itind)
            endif

!      ...  collimated beam
            zrefb(kp) = max(f_zero, min(f_one,                          &
     &                  (za2 - za1*(f_one - zb2))/(f_one + za2) ))
            ztrab(kp) = max(f_zero, min(f_one, f_one-zrefb(kp) ))

!      ...  isotropic incidence
            zrefd(kp) = max(f_zero, min(f_one, za2/(f_one + za2) ))
            ztrad(kp) = max(f_zero, min(f_one, f_one-zrefd(kp) ))

          else                          ! for non-conservative scattering
            za1 = zgam1*zgam4 + zgam2*zgam3
            za2 = zgam1*zgam3 + zgam2*zgam4
            zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
            zrk2= 2.0 * zrk

            zrp  = zrk * cosz
            zrp1 = f_one + zrp
            zrm1 = f_one - zrp
            zrpp1= f_one - zrp*zrp
            zrpp = sign( max(flimit, abs(zrpp1)), zrpp1 )    ! avoid numerical singularity
            zrkg1= zrk + zgam1
            zrkg3= zrk * zgam3
            zrkg4= zrk * zgam4

            zr1  = zrm1 * (za2 + zrkg3)
            zr2  = zrp1 * (za2 - zrkg3)
            zr3  = zrk2 * (zgam3 - za2*cosz)
            zr4  = zrpp * zrkg1
            zr5  = zrpp * (zrk - zgam1)

            zt1  = zrp1 * (za1 + zrkg4)
            zt2  = zrm1 * (za1 - zrkg4)
            zt3  = zrk2 * (zgam4 + za1*cosz)

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

            zb1 = min ( zrk*ztau1, 500.0 )
            if ( zb1 <= od_lo ) then
              zexm1 = f_one - zb1 + 0.5*zb1*zb1
            else
              ftind = zb1 / (bpade + zb1)
              itind = ftind*NTBMX + 0.5
              zexm1 = exp_tbl(itind)
            endif
            zexp1 = f_one / zexm1

            zb2 = min ( sntz*ztau1, 500.0 )
            if ( zb2 <= od_lo ) then
              zexm2 = f_one - zb2 + 0.5*zb2*zb2
            else
              ftind = zb2 / (bpade + zb2)
              itind = ftind*NTBMX + 0.5
              zexm2 = exp_tbl(itind)
            endif
            zexp2 = f_one / zexm2
            ze1r45 = zr4*zexp1 + zr5*zexm1

!      ...  collimated beam
            if (ze1r45>=-eps1 .and. ze1r45<=eps1) then
              zrefb(kp) = eps1
              ztrab(kp) = zexm2
            else
              zden1 = zssa1 / ze1r45
              zrefb(kp) = max(f_zero, min(f_one,                        &
     &                    (zr1*zexp1 - zr2*zexm1 - zr3*zexm2)*zden1 ))
              ztrab(kp) = max(f_zero, min(f_one, zexm2*(f_one           &
     &                  - (zt1*zexp1 - zt2*zexm1 - zt3*zexp2)*zden1) ))
            endif

!      ...  diffuse beam
            zden1 = zr4 / (ze1r45 * zrkg1)
            zrefd(kp) = max(f_zero, min(f_one,                          &
     &                  zgam2*(zexp1 - zexm1)*zden1 ))
            ztrad(kp) = max(f_zero, min(f_one, zrk2*zden1 ))
          endif    ! end if_zssaw_block

!>  -  Calculate direct beam transmittance. use exponential lookup table
!! for transmittance, or expansion of exponential for low optical depth.

          zr1 = ztau1 * sntz
          if ( zr1 <= od_lo ) then
            zexp3 = f_one - zr1 + 0.5*zr1*zr1
          else
            ftind = zr1 / (bpade + zr1)
            itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
            zexp3 = exp_tbl(itind)
          endif

          ztdbt(k)  = zexp3 * ztdbt(kp)
          zldbt(kp) = zexp3

!>  - Calculate pre-delta-scaling clear and cloudy direct beam transmittance.
!           (must use 'orig', unscaled cloud optical depth)

          zr1 = ztau0 * sntz
          if ( zr1 <= od_lo ) then
            zexp4 = f_one - zr1 + 0.5*zr1*zr1
          else
            ftind = zr1 / (bpade + zr1)
            itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
            zexp4 = exp_tbl(itind)
          endif

          zldbt0(k) = zexp4
          ztdbt0 = zexp4 * ztdbt0
        enddo    ! end do_k_loop

!> -# Call vrtqdr(), to compute the upward and downward radiation fluxes.
        call vrtqdr                                                     &
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!> -# Compute upward and downward fluxes at levels.
        do k = 1, nlp1
          fxup0(k,ib) = fxup0(k,ib) + zsolar*zfu(k)
          fxdn0(k,ib) = fxdn0(k,ib) + zsolar*zfd(k)
        enddo

!> -# Compute surface downward beam/diffused flux components.
        zb1 = zsolar*ztdbt0
        zb2 = zsolar*(zfd(1) - ztdbt0)

        if (ibd /= 0) then
          sfbm0(ibd) = sfbm0(ibd) + zb1
          sfdf0(ibd) = sfdf0(ibd) + zb2
        else
          zf1 = 0.5 * zb1
          zf2 = 0.5 * zb2
          sfbm0(1) = sfbm0(1) + zf1
          sfdf0(1) = sfdf0(1) + zf2
          sfbm0(2) = sfbm0(2) + zf1
          sfdf0(2) = sfdf0(2) + zf2
        endif
!       sfbm0(ibd) = sfbm0(ibd) + zsolar*ztdbt0
!       sfdf0(ibd) = sfdf0(ibd) + zsolar*(zfd(1) - ztdbt0)

!> -# Compute total sky optical parameters, layer reflectance and
!!    transmittance.
!    - Set up toa direct beam and surface values (beam and diff)
!    - Delta scaling for total-sky condition
!    - General two-stream expressions for physparam::iswmode
!    - Compute homogeneous reflectance and transmittance for
!      conservative scattering and non-conservative scattering
!    - Pre-delta-scaling clear and cloudy direct beam transmittance
!    - Call swflux() to compute the upward and downward radiation fluxes

        if ( cf1 > eps ) then

!>  - Set up toa direct beam and surface values (beam and diff).
          ztdbt0 = f_one
          zldbt(1) = f_zero

          do k = nlay, 1, -1
            kp = k + 1
            zc0 = f_one - cldfrc(k)
            zc1 = cldfrc(k)
            if ( zc1 > ftiny ) then          ! it is a cloudy-layer

              ztau0 = ztaus(k) + taucw(k,ib)
              zssa0 = zssas(k) + ssacw(k,ib)
              zasy0 = zasys(k) + asycw(k,ib)
              zssaw = min(oneminus, zssa0 / ztau0)
              zasyw = zasy0 / max(ftiny, zssa0)

!>  - Perform delta scaling for total-sky condition.
              za1 = zasyw * zasyw
              za2 = zssaw * za1

              ztau1 = (f_one - za2) * ztau0
              zssa1 = (zssaw - za2) / (f_one - za2)
!org          zasy1 = (zasyw - za1) / (f_one - za1)
              zasy1 = zasyw / (f_one + zasyw)
              zasy3 = 0.75 * zasy1

!>  - Perform general two-stream expressions:
!!\n  control parameters in module "physparam"
!!\n    iswmode - control flag for 2-stream transfer schemes
!!\n              = 1 delta-eddington    (joseph et al., 1976)
!!\n              = 2 pifm               (zdunkowski et al., 1980)
!!\n              = 3 discrete ordinates (liou, 1973)

              if ( iswmode == 1 ) then
                zgam1 = 1.75 - zssa1 * (f_one + zasy3)
                zgam2 =-0.25 + zssa1 * (f_one - zasy3)
                zgam3 = 0.5  - zasy3 * cosz
              elseif ( iswmode == 2 ) then               ! pifm
                zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
                zgam2 = 0.75* zssa1 * (f_one- zasy1)
                zgam3 = 0.5 - zasy3 * cosz
              elseif ( iswmode == 3 ) then               ! discrete ordinates
                zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
                zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
                zgam3 = (1.0 - zsr3 * zasy1 * cosz) * 0.5
              endif
              zgam4 = f_one - zgam3

              zrefb1 = zrefb(kp)
              zrefd1 = zrefd(kp)
              ztrab1 = ztrab(kp)
              ztrad1 = ztrad(kp)

!>  - Compute homogeneous reflectance and transmittance for both conservative
!! and non-conservative scattering.

              if ( zssaw >= zcrit ) then    ! for conservative scattering
                za1 = zgam1 * cosz - zgam3
                za2 = zgam1 * ztau1

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

                zb1 = min ( ztau1*sntz , 500.0 )
                if ( zb1 <= od_lo ) then
                  zb2 = f_one - zb1 + 0.5*zb1*zb1
                else
                  ftind = zb1 / (bpade + zb1)
                  itind = ftind*NTBMX + 0.5
                  zb2 = exp_tbl(itind)
                endif

!      ...  collimated beam
                zrefb(kp) = max(f_zero, min(f_one,                      &
     &                      (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                ztrab(kp) = max(f_zero, min(f_one, f_one-zrefb(kp)))

!      ...  isotropic incidence
                zrefd(kp) = max(f_zero, min(f_one, za2 / (f_one+za2) ))
                ztrad(kp) = max(f_zero, min(f_one, f_one - zrefd(kp) ))

              else                          ! for non-conservative scattering
                za1 = zgam1*zgam4 + zgam2*zgam3
                za2 = zgam1*zgam3 + zgam2*zgam4
                zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
                zrk2= 2.0 * zrk

                zrp  = zrk * cosz
                zrp1 = f_one + zrp
                zrm1 = f_one - zrp
                zrpp1= f_one - zrp*zrp
                zrpp = sign( max(flimit, abs(zrpp1)), zrpp1 )    ! avoid numerical singularity
                zrkg1= zrk + zgam1
                zrkg3= zrk * zgam3
                zrkg4= zrk * zgam4

                zr1  = zrm1 * (za2 + zrkg3)
                zr2  = zrp1 * (za2 - zrkg3)
                zr3  = zrk2 * (zgam3 - za2*cosz)
                zr4  = zrpp * zrkg1
                zr5  = zrpp * (zrk - zgam1)

                zt1  = zrp1 * (za1 + zrkg4)
                zt2  = zrm1 * (za1 - zrkg4)
                zt3  = zrk2 * (zgam4 + za1*cosz)

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

                zb1 = min ( zrk*ztau1, 500.0 )
                if ( zb1 <= od_lo ) then
                  zexm1 = f_one - zb1 + 0.5*zb1*zb1
                else
                  ftind = zb1 / (bpade + zb1)
                  itind = ftind*NTBMX + 0.5
                  zexm1 = exp_tbl(itind)
                endif
                zexp1 = f_one / zexm1

                zb2 = min ( ztau1*sntz, 500.0 )
                if ( zb2 <= od_lo ) then
                  zexm2 = f_one - zb2 + 0.5*zb2*zb2
                else
                  ftind = zb2 / (bpade + zb2)
                  itind = ftind*NTBMX + 0.5
                  zexm2 = exp_tbl(itind)
                endif
                zexp2 = f_one / zexm2
                ze1r45 = zr4*zexp1 + zr5*zexm1

!      ...  collimated beam
                if ( ze1r45>=-eps1 .and. ze1r45<=eps1 ) then
                  zrefb(kp) = eps1
                  ztrab(kp) = zexm2
                else
                  zden1 = zssa1 / ze1r45
                  zrefb(kp) = max(f_zero, min(f_one,                    &
     &                        (zr1*zexp1-zr2*zexm1-zr3*zexm2)*zden1 ))
                  ztrab(kp) = max(f_zero, min(f_one, zexm2*(f_one -     &
     &                        (zt1*zexp1-zt2*zexm1-zt3*zexp2)*zden1) ))
                endif

!      ...  diffuse beam
                zden1 = zr4 / (ze1r45 * zrkg1)
                zrefd(kp) = max(f_zero, min(f_one,                      &
     &                      zgam2*(zexp1 - zexm1)*zden1 ))
                ztrad(kp) = max(f_zero, min(f_one, zrk2*zden1 ))
              endif    ! end if_zssaw_block

!  --- ...  combine clear and cloudy contributions for total sky
!           and calculate direct beam transmittances

              zrefb(kp) = zc0*zrefb1 + zc1*zrefb(kp)
              zrefd(kp) = zc0*zrefd1 + zc1*zrefd(kp)
              ztrab(kp) = zc0*ztrab1 + zc1*ztrab(kp)
              ztrad(kp) = zc0*ztrad1 + zc1*ztrad(kp)

!  --- ...  direct beam transmittance. use exponential lookup table
!           for transmittance, or expansion of exponential for low
!           optical depth

              zr1 = ztau1 * sntz
              if ( zr1 <= od_lo ) then
                zexp3 = f_one - zr1 + 0.5*zr1*zr1
              else
                ftind = zr1 / (bpade + zr1)
                itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
                zexp3 = exp_tbl(itind)
              endif

              zldbt(kp) = zc0*zldbt(kp) + zc1*zexp3
              ztdbt(k) = zldbt(kp) * ztdbt(kp)

!>  - Calculate pre-delta-scaling clear and cloudy direct beam transmittance.
!           (must use 'orig', unscaled cloud optical depth)

              zr1 = ztau0 * sntz
              if ( zr1 <= od_lo ) then
                zexp4 = f_one - zr1 + 0.5*zr1*zr1
              else
                ftind = zr1 / (bpade + zr1)
                itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
                zexp4 = exp_tbl(itind)
              endif

              ztdbt0 = (zc0*zldbt0(k) + zc1*zexp4) * ztdbt0

            else     ! if_zc1_block  ---  it is a clear layer

!  --- ...  direct beam transmittance
              ztdbt(k) = zldbt(kp) * ztdbt(kp)

!  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
              ztdbt0 = zldbt0(k) * ztdbt0

            endif    ! end if_zc1_block
          enddo   ! end do_k_loop

!> -# Call vrtqdr(), to compute the upward and downward radiation fluxes.

          call vrtqdr                                                   &
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!> -# Compute upward and downward fluxes at levels.
          do k = 1, nlp1
            fxupc(k,ib) = fxupc(k,ib) + zsolar*zfu(k)
            fxdnc(k,ib) = fxdnc(k,ib) + zsolar*zfd(k)
          enddo

!> -# Process and save outputs.
!!  - surface downward beam/diffused flux components
          zb1 = zsolar*ztdbt0
          zb2 = zsolar*(zfd(1) - ztdbt0)

          if (ibd /= 0) then
            sfbmc(ibd) = sfbmc(ibd) + zb1
            sfdfc(ibd) = sfdfc(ibd) + zb2
          else
            zf1 = 0.5 * zb1
            zf2 = 0.5 * zb2
            sfbmc(1) = sfbmc(1) + zf1
            sfdfc(1) = sfdfc(1) + zf2
            sfbmc(2) = sfbmc(2) + zf1
            sfdfc(2) = sfdfc(2) + zf2
          endif
!         sfbmc(ibd) = sfbmc(ibd) + zsolar*ztdbt0
!         sfdfc(ibd) = sfdfc(ibd) + zsolar*(zfd(1) - ztdbt0)

        endif      ! end if_cf1_block

      enddo  lab_do_jg

!  --- ...  end of g-point loop

      do ib = 1, nbdsw
        ftoadc = ftoadc + fxdn0(nlp1,ib)
        ftoau0 = ftoau0 + fxup0(nlp1,ib)
        fsfcu0 = fsfcu0 + fxup0(1,ib)
        fsfcd0 = fsfcd0 + fxdn0(1,ib)
      enddo

!>  - uv-b surface downward flux
      ibd = nuvb - nblow + 1
      suvbf0 = fxdn0(1,ibd)

      if ( cf1 <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
        do ib = 1, nbdsw
          do k = 1, nlp1
            fxupc(k,ib) = fxup0(k,ib)
            fxdnc(k,ib) = fxdn0(k,ib)
          enddo
        enddo

        ftoauc = ftoau0
        fsfcuc = fsfcu0
        fsfcdc = fsfcd0

!>  - surface downward beam/diffused flux components
        sfbmc(1) = sfbm0(1)
        sfdfc(1) = sfdf0(1)
        sfbmc(2) = sfbm0(2)
        sfdfc(2) = sfdf0(2)

!>  - uv-b surface downward flux
        suvbfc = suvbf0
      else                        ! cloudy column, compute total-sky fluxes
        do ib = 1, nbdsw
          do k = 1, nlp1
            fxupc(k,ib) = cf1*fxupc(k,ib) + cf0*fxup0(k,ib)
            fxdnc(k,ib) = cf1*fxdnc(k,ib) + cf0*fxdn0(k,ib)
          enddo
        enddo

        do ib = 1, nbdsw
          ftoauc = ftoauc + fxupc(nlp1,ib)
          fsfcuc = fsfcuc + fxupc(1,ib)
          fsfcdc = fsfcdc + fxdnc(1,ib)
        enddo

!>  - uv-b surface downward flux
        suvbfc = fxdnc(1,ibd)

!>  - surface downward beam/diffused flux components
        sfbmc(1) = cf1*sfbmc(1) + cf0*sfbm0(1)
        sfbmc(2) = cf1*sfbmc(2) + cf0*sfbm0(2)
        sfdfc(1) = cf1*sfdfc(1) + cf0*sfdf0(1)
        sfdfc(2) = cf1*sfdfc(2) + cf0*sfdf0(2)
      endif    ! end if_cf1_block

      return
!...................................
      end subroutine spcvrtc
!-----------------------------------
!> @}

!>\ingroup module_radsw_main
!> This subroutine computes the shortwave radiative fluxes using
!! two-stream method of h. barder and mcica,the monte-carlo independent
!! column approximation, for the representation of sub-grid cloud
!! variability (i.e. cloud overlap).
!!\param ssolar        incoming solar flux at top
!!\param cosz          cosine solar zenith angle
!!\param sntz          secant solar zenith angle
!!\param albbm         surface albedo for direct beam radiation
!!\param albdf         surface albedo for diffused radiation
!!\param sfluxzen      spectral distribution of incoming solar flux
!!\param cldfmc        layer cloud fraction for g-point
!!\param cf1            >0: cloudy sky, otherwise: clear sky
!!\param cf0           =1-cf1
!!\param taug          spectral optical depth for gases
!!\param taur          optical depth for rayleigh scattering
!!\param tauae         aerosols optical depth
!!\param ssaae         aerosols single scattering albedo
!!\param asyae         aerosols asymmetry factor
!!\param taucw         weighted cloud optical depth
!!\param ssacw         weighted cloud single scat albedo
!!\param asycw         weighted cloud asymmetry factor
!!\param nlay,nlp1     number of layers/levels
!!\param fxupc         tot sky upward flux
!!\param fxdnc         tot sky downward flux
!!\param fxup0         clr sky upward flux
!!\param fxdn0         clr sky downward flux
!!\param ftoauc        tot sky toa upwd flux
!!\param ftoau0        clr sky toa upwd flux
!!\param ftoadc        toa downward (incoming) solar flux
!!\param fsfcuc        tot sky sfc upwd flux
!!\param fsfcu0        clr sky sfc upwd flux
!!\param fsfcdc        tot sky sfc dnwd flux
!!\param fsfcd0        clr sky sfc dnwd flux
!!\param sfbmc         tot sky sfc dnwd beam flux (nir/uv+vis)
!!\param sfdfc         tot sky sfc dnwd diff flux (nir/uv+vis)
!!\param sfbm0         clr sky sfc dnwd beam flux (nir/uv+vis)
!!\param sfdf0         clr sky sfc dnwd diff flux (nir/uv+vis)
!!\param suvbfc        tot sky sfc dnwd uv-b flux
!!\param suvbf0        clr sky sfc dnwd uv-b flux
!>\section spcvrtm_gen spcvrtm General Algorithm
!! @{
!-----------------------------------
      subroutine spcvrtm                                                &
     &     ( ssolar,cosz,sntz,albbm,albdf,sfluxzen,cldfmc,              &   !  ---  inputs
     &       cf1,cf0,taug,taur,tauae,ssaae,asyae,taucw,ssacw,asycw,     &
     &       nlay, nlp1,                                                &
     &       fxupc,fxdnc,fxup0,fxdn0,                                   &   !  ---  outputs
     &       ftoauc,ftoau0,ftoadc,fsfcuc,fsfcu0,fsfcdc,fsfcd0,          &
     &       sfbmc,sfdfc,sfbm0,sfdf0,suvbfc,suvbf0                      &
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
!   purpose:  computes the shortwave radiative fluxes using two-stream  !
!             method of h. barker and mcica, the monte-carlo independent!
!             column approximation, for the representation of sub-grid  !
!             cloud variability (i.e. cloud overlap).                   !
!                                                                       !
!   subprograms called:  vrtqdr                                         !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                        size  !
!    ssolar  - real, incoming solar flux at top                    1    !
!    cosz    - real, cosine solar zenith angle                     1    !
!    sntz    - real, secant solar zenith angle                     1    !
!    albbm   - real, surface albedo for direct beam radiation      2    !
!    albdf   - real, surface albedo for diffused radiation         2    !
!    sfluxzen- real, spectral distribution of incoming solar flux ngptsw!
!    cldfmc  - real, layer cloud fraction for g-point        nlay*ngptsw!
!    cf1     - real, >0: cloudy sky, otherwise: clear sky          1    !
!    cf0     - real, =1-cf1                                        1    !
!    taug    - real, spectral optical depth for gases        nlay*ngptsw!
!    taur    - real, optical depth for rayleigh scattering   nlay*ngptsw!
!    tauae   - real, aerosols optical depth                  nlay*nbdsw !
!    ssaae   - real, aerosols single scattering albedo       nlay*nbdsw !
!    asyae   - real, aerosols asymmetry factor               nlay*nbdsw !
!    taucw   - real, weighted cloud optical depth            nlay*nbdsw !
!    ssacw   - real, weighted cloud single scat albedo       nlay*nbdsw !
!    asycw   - real, weighted cloud asymmetry factor         nlay*nbdsw !
!    nlay,nlp1 - integer,  number of layers/levels                 1    !
!                                                                       !
!  output variables:                                                    !
!    fxupc   - real, tot sky upward flux                     nlp1*nbdsw !
!    fxdnc   - real, tot sky downward flux                   nlp1*nbdsw !
!    fxup0   - real, clr sky upward flux                     nlp1*nbdsw !
!    fxdn0   - real, clr sky downward flux                   nlp1*nbdsw !
!    ftoauc  - real, tot sky toa upwd flux                         1    !
!    ftoau0  - real, clr sky toa upwd flux                         1    !
!    ftoadc  - real, toa downward (incoming) solar flux            1    !
!    fsfcuc  - real, tot sky sfc upwd flux                         1    !
!    fsfcu0  - real, clr sky sfc upwd flux                         1    !
!    fsfcdc  - real, tot sky sfc dnwd flux                         1    !
!    fsfcd0  - real, clr sky sfc dnwd flux                         1    !
!    sfbmc   - real, tot sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdfc   - real, tot sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    sfbm0   - real, clr sky sfc dnwd beam flux (nir/uv+vis)       2    !
!    sfdf0   - real, clr sky sfc dnwd diff flux (nir/uv+vis)       2    !
!    suvbfc  - real, tot sky sfc dnwd uv-b flux                    1    !
!    suvbf0  - real, clr sky sfc dnwd uv-b flux                    1    !
!                                                                       !
!  internal variables:                                                  !
!    zrefb   - real, direct beam reflectivity for clear/cloudy    nlp1  !
!    zrefd   - real, diffuse reflectivity for clear/cloudy        nlp1  !
!    ztrab   - real, direct beam transmissivity for clear/cloudy  nlp1  !
!    ztrad   - real, diffuse transmissivity for clear/cloudy      nlp1  !
!    zldbt   - real, layer beam transmittance for clear/cloudy    nlp1  !
!    ztdbt   - real, lev total beam transmittance for clr/cld     nlp1  !
!                                                                       !
!  control parameters in module "physparam"                             !
!    iswmode - control flag for 2-stream transfer schemes               !
!              = 1 delta-eddington    (joseph et al., 1976)             !
!              = 2 pifm               (zdunkowski et al., 1980)         !
!              = 3 discrete ordinates (liou, 1973)                      !
!                                                                       !
!  *******************************************************************  !
!  original code description                                            !
!                                                                       !
!  method:                                                              !
!  -------                                                              !
!     standard delta-eddington, p.i.f.m., or d.o.m. layer calculations. !
!     kmodts  = 1 eddington (joseph et al., 1976)                       !
!             = 2 pifm (zdunkowski et al., 1980)                        !
!             = 3 discrete ordinates (liou, 1973)                       !
!                                                                       !
!  modifications:                                                       !
!  --------------                                                       !
!   original: h. barker                                                 !
!   revision: merge with rrtmg_sw: j.-j.morcrette, ecmwf, feb 2003      !
!   revision: add adjustment for earth/sun distance:mjiacono,aer,oct2003!
!   revision: bug fix for use of palbp and palbd: mjiacono, aer, nov2003!
!   revision: bug fix to apply delta scaling to clear sky: aer, dec2004 !
!   revision: code modified so that delta scaling is not done in cloudy !
!             profiles if routine cldprop is used; delta scaling can be !
!             applied by swithcing code below if cldprop is not used to !
!             get cloud properties. aer, jan 2005                       !
!   revision: uniform formatting for rrtmg: mjiacono, aer, jul 2006     !
!   revision: use exponential lookup table for transmittance: mjiacono, !
!             aer, aug 2007                                             !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  constant parameters:
      real (kind=kind_phys), parameter :: zcrit = 0.9999995 ! thresold for conservative scattering
      real (kind=kind_phys), parameter :: zsr3  = sqrt(3.0)
      real (kind=kind_phys), parameter :: od_lo = 0.06
      real (kind=kind_phys), parameter :: eps1  = 1.0e-8

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nlay,ngptsw), intent(in) ::      &
     &       taug, taur, cldfmc
      real (kind=kind_phys), dimension(nlay,nbdsw),  intent(in) ::      &
     &       taucw, ssacw, asycw, tauae, ssaae, asyae

      real (kind=kind_phys), dimension(ngptsw), intent(in) :: sfluxzen

      real (kind=kind_phys), dimension(2),  intent(in) :: albbm, albdf

      real (kind=kind_phys), intent(in) :: cosz, sntz, cf1, cf0, ssolar

!  ---  outputs:
      real (kind=kind_phys), dimension(nlp1,nbdsw), intent(out) ::      &
     &       fxupc, fxdnc, fxup0, fxdn0

      real (kind=kind_phys), dimension(2), intent(out) :: sfbmc, sfdfc, &
     &       sfbm0, sfdf0

      real (kind=kind_phys), intent(out) :: suvbfc, suvbf0, ftoadc,     &
     &       ftoauc, ftoau0, fsfcuc, fsfcu0, fsfcdc, fsfcd0

!  ---  locals:
      real (kind=kind_phys), dimension(nlay) :: ztaus, zssas, zasys,    &
     &       zldbt0

      real (kind=kind_phys), dimension(nlp1) :: zrefb, zrefd, ztrab,    &
     &       ztrad, ztdbt, zldbt, zfu, zfd

      real (kind=kind_phys) :: ztau1, zssa1, zasy1, ztau0, zssa0,       &
     &       zasy0, zasy3, zssaw, zasyw, zgam1, zgam2, zgam3, zgam4,    &
     &       za1, za2, zb1, zb2, zrk, zrk2, zrp, zrp1, zrm1, zrpp,      &
     &       zrkg1, zrkg3, zrkg4, zexp1, zexm1, zexp2, zexm2, zden1,    &
     &       zexp3, zexp4, ze1r45, ftind, zsolar, ztdbt0, zr1, zr2,     &
     &       zr3, zr4, zr5, zt1, zt2, zt3, zf1, zf2, zrpp1

      integer :: ib, ibd, jb, jg, k, kp, itind
!
!===> ...  begin here
!
!> -# Initialize output fluxes.

      do ib = 1, nbdsw
        do k = 1, nlp1
          fxdnc(k,ib) = f_zero
          fxupc(k,ib) = f_zero
          fxdn0(k,ib) = f_zero
          fxup0(k,ib) = f_zero
        enddo
      enddo

      ftoadc = f_zero
      ftoauc = f_zero
      ftoau0 = f_zero
      fsfcuc = f_zero
      fsfcu0 = f_zero
      fsfcdc = f_zero
      fsfcd0 = f_zero

!! --- ...  uv-b surface downward fluxes
      suvbfc  = f_zero
      suvbf0  = f_zero

!! --- ...  output surface flux components
      sfbmc(1) = f_zero
      sfbmc(2) = f_zero
      sfdfc(1) = f_zero
      sfdfc(2) = f_zero
      sfbm0(1) = f_zero
      sfbm0(2) = f_zero
      sfdf0(1) = f_zero
      sfdf0(2) = f_zero

!> -# Loop over all g-points in each band.

      lab_do_jg : do jg = 1, ngptsw

        jb = NGB(jg)
        ib = jb + 1 - nblow
        ibd = idxsfc(jb)         ! spectral band index

        zsolar = ssolar * sfluxzen(jg)

!> -# Set up toa direct beam and surface values (beam and diff).

        ztdbt(nlp1) = f_one
        ztdbt0   = f_one

        zldbt(1) = f_zero
        if (ibd /= 0) then
          zrefb(1) = albbm(ibd)
          zrefd(1) = albdf(ibd)
        else
          zrefb(1) = 0.5 * (albbm(1) + albbm(2))
          zrefd(1) = 0.5 * (albdf(1) + albdf(2))
        endif
        ztrab(1) = f_zero
        ztrad(1) = f_zero

!> -# Compute clear-sky optical parameters, layer reflectance and
!!    transmittance.
!    - Set up toa direct beam and surface values (beam and diff)
!    - Delta scaling for clear-sky condition
!    - General two-stream expressions for physparam::iswmode
!    - Compute homogeneous reflectance and transmittance for both
!      conservative and non-conservative scattering
!    - Pre-delta-scaling clear and cloudy direct beam transmittance
!    - Call swflux() to compute the upward and downward radiation fluxes

        do k = nlay, 1, -1
          kp = k + 1

          ztau0 = max( ftiny, taur(k,jg)+taug(k,jg)+tauae(k,ib) )
          zssa0 = taur(k,jg) + tauae(k,ib)*ssaae(k,ib)
          zasy0 = asyae(k,ib)*ssaae(k,ib)*tauae(k,ib)
          zssaw = min( oneminus, zssa0 / ztau0 )
          zasyw = zasy0 / max( ftiny, zssa0 )

!>  - Saving clear-sky quantities for later total-sky usage.
          ztaus(k) = ztau0
          zssas(k) = zssa0
          zasys(k) = zasy0

!>  - Delta scaling for clear-sky condition.
          za1 = zasyw * zasyw
          za2 = zssaw * za1

          ztau1 = (f_one - za2) * ztau0
          zssa1 = (zssaw - za2) / (f_one - za2)
!org      zasy1 = (zasyw - za1) / (f_one - za1)   ! this line is replaced by the next
          zasy1 = zasyw / (f_one + zasyw)         ! to reduce truncation error
          zasy3 = 0.75 * zasy1

!>  - Perform general two-stream expressions:
!!\n control parameters in module "physparam" 
!!\n iswmode - control flag for 2-stream transfer schemes 
!!\n           = 1 delta-eddington (joseph et al., 1976) 
!!\n           = 2 pifm (zdunkowski et al., 1980) 
!!\n           = 3 discrete ordinates (liou, 1973)
          if ( iswmode == 1 ) then
            zgam1 = 1.75 - zssa1 * (f_one + zasy3)
            zgam2 =-0.25 + zssa1 * (f_one - zasy3)
            zgam3 = 0.5  - zasy3 * cosz
          elseif ( iswmode == 2 ) then               ! pifm
            zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
            zgam2 = 0.75* zssa1 * (f_one- zasy1)
            zgam3 = 0.5 - zasy3 * cosz
          elseif ( iswmode == 3 ) then               ! discrete ordinates
            zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
            zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
            zgam3 = (1.0 - zsr3 * zasy1 * cosz) * 0.5
          endif
          zgam4 = f_one - zgam3

!>  - Compute homogeneous reflectance and transmittance.

          if ( zssaw >= zcrit ) then    ! for conservative scattering
            za1 = zgam1 * cosz - zgam3
            za2 = zgam1 * ztau1

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

            zb1 = min ( ztau1*sntz , 500.0 )
            if ( zb1 <= od_lo ) then
              zb2 = f_one - zb1 + 0.5*zb1*zb1
            else
              ftind = zb1 / (bpade + zb1)
              itind = ftind*NTBMX + 0.5
              zb2 = exp_tbl(itind)
            endif

!      ...  collimated beam
            zrefb(kp) = max(f_zero, min(f_one,                          &
     &                  (za2 - za1*(f_one - zb2))/(f_one + za2) ))
            ztrab(kp) = max(f_zero, min(f_one, f_one-zrefb(kp) ))

!      ...  isotropic incidence
            zrefd(kp) = max(f_zero, min(f_one, za2/(f_one + za2) ))
            ztrad(kp) = max(f_zero, min(f_one, f_one-zrefd(kp) ))

          else                          ! for non-conservative scattering
            za1 = zgam1*zgam4 + zgam2*zgam3
            za2 = zgam1*zgam3 + zgam2*zgam4
            zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
            zrk2= 2.0 * zrk

            zrp  = zrk * cosz
            zrp1 = f_one + zrp
            zrm1 = f_one - zrp
            zrpp1= f_one - zrp*zrp
            zrpp = sign( max(flimit, abs(zrpp1)), zrpp1 )    ! avoid numerical singularity
            zrkg1= zrk + zgam1
            zrkg3= zrk * zgam3
            zrkg4= zrk * zgam4

            zr1  = zrm1 * (za2 + zrkg3)
            zr2  = zrp1 * (za2 - zrkg3)
            zr3  = zrk2 * (zgam3 - za2*cosz)
            zr4  = zrpp * zrkg1
            zr5  = zrpp * (zrk - zgam1)

            zt1  = zrp1 * (za1 + zrkg4)
            zt2  = zrm1 * (za1 - zrkg4)
            zt3  = zrk2 * (zgam4 + za1*cosz)

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

            zb1 = min ( zrk*ztau1, 500.0 )
            if ( zb1 <= od_lo ) then
              zexm1 = f_one - zb1 + 0.5*zb1*zb1
            else
              ftind = zb1 / (bpade + zb1)
              itind = ftind*NTBMX + 0.5
              zexm1 = exp_tbl(itind)
            endif
            zexp1 = f_one / zexm1

            zb2 = min ( sntz*ztau1, 500.0 )
            if ( zb2 <= od_lo ) then
              zexm2 = f_one - zb2 + 0.5*zb2*zb2
            else
              ftind = zb2 / (bpade + zb2)
              itind = ftind*NTBMX + 0.5
              zexm2 = exp_tbl(itind)
            endif
            zexp2 = f_one / zexm2
            ze1r45 = zr4*zexp1 + zr5*zexm1

!      ...  collimated beam
            if (ze1r45>=-eps1 .and. ze1r45<=eps1) then
              zrefb(kp) = eps1
              ztrab(kp) = zexm2
            else
              zden1 = zssa1 / ze1r45
              zrefb(kp) = max(f_zero, min(f_one,                        &
     &                    (zr1*zexp1 - zr2*zexm1 - zr3*zexm2)*zden1 ))
              ztrab(kp) = max(f_zero, min(f_one, zexm2*(f_one           &
     &                  - (zt1*zexp1 - zt2*zexm1 - zt3*zexp2)*zden1) ))
            endif

!      ...  diffuse beam
            zden1 = zr4 / (ze1r45 * zrkg1)
            zrefd(kp) = max(f_zero, min(f_one,                          &
     &                  zgam2*(zexp1 - zexm1)*zden1 ))
            ztrad(kp) = max(f_zero, min(f_one, zrk2*zden1 ))
          endif    ! end if_zssaw_block

!>  - Calculate direct beam transmittance. use exponential lookup table
!! for transmittance, or expansion of exponential for low optical depth.

          zr1 = ztau1 * sntz
          if ( zr1 <= od_lo ) then
            zexp3 = f_one - zr1 + 0.5*zr1*zr1
          else
            ftind = zr1 / (bpade + zr1)
            itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
            zexp3 = exp_tbl(itind)
          endif

          ztdbt(k)  = zexp3 * ztdbt(kp)
          zldbt(kp) = zexp3

!>  - Calculate pre-delta-scaling clear and cloudy direct beam transmittance.
!           (must use 'orig', unscaled cloud optical depth)

          zr1 = ztau0 * sntz
          if ( zr1 <= od_lo ) then
            zexp4 = f_one - zr1 + 0.5*zr1*zr1
          else
            ftind = zr1 / (bpade + zr1)
            itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
            zexp4 = exp_tbl(itind)
          endif

          zldbt0(k) = zexp4
          ztdbt0 = zexp4 * ztdbt0
        enddo    ! end do_k_loop

!> -# Call vrtqdr(), to compute the upward and downward radiation fluxes.
        call vrtqdr                                                     &
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!> -# Compute upward and downward fluxes at levels.
        do k = 1, nlp1
          fxup0(k,ib) = fxup0(k,ib) + zsolar*zfu(k)
          fxdn0(k,ib) = fxdn0(k,ib) + zsolar*zfd(k)
        enddo

!> -# Compute surface downward beam/diffuse flux components.
        zb1 = zsolar*ztdbt0
        zb2 = zsolar*(zfd(1) - ztdbt0)

        if (ibd /= 0) then
          sfbm0(ibd) = sfbm0(ibd) + zb1
          sfdf0(ibd) = sfdf0(ibd) + zb2
        else
          zf1 = 0.5 * zb1
          zf2 = 0.5 * zb2
          sfbm0(1) = sfbm0(1) + zf1
          sfdf0(1) = sfdf0(1) + zf2
          sfbm0(2) = sfbm0(2) + zf1
          sfdf0(2) = sfdf0(2) + zf2
        endif
!       sfbm0(ibd) = sfbm0(ibd) + zsolar*ztdbt0
!       sfdf0(ibd) = sfdf0(ibd) + zsolar*(zfd(1) - ztdbt0)

!> -# Compute total sky optical parameters, layer reflectance and
!!    transmittance.
!    - Set up toa direct beam and surface values (beam and diff)
!    - Delta scaling for total-sky condition
!    - General two-stream expressions for physparam::iswmode
!    - Compute homogeneous reflectance and transmittance for
!      conservative scattering and non-conservative scattering
!    - Pre-delta-scaling clear and cloudy direct beam transmittance
!    - Call swflux() to compute the upward and downward radiation fluxes

        if ( cf1 > eps ) then

!>  - Set up toa direct beam and surface values (beam and diff).
          ztdbt0 = f_one
          zldbt(1) = f_zero

          do k = nlay, 1, -1
            kp = k + 1
            if ( cldfmc(k,jg) > ftiny ) then      ! it is a cloudy-layer

              ztau0 = ztaus(k) + taucw(k,ib)
              zssa0 = zssas(k) + ssacw(k,ib)
              zasy0 = zasys(k) + asycw(k,ib)
              zssaw = min(oneminus, zssa0 / ztau0)
              zasyw = zasy0 / max(ftiny, zssa0)

!>  - Perform delta scaling for total-sky condition.
              za1 = zasyw * zasyw
              za2 = zssaw * za1

              ztau1 = (f_one - za2) * ztau0
              zssa1 = (zssaw - za2) / (f_one - za2)
!org          zasy1 = (zasyw - za1) / (f_one - za1)
              zasy1 = zasyw / (f_one + zasyw)
              zasy3 = 0.75 * zasy1

!>  - Perform general two-stream expressions.
              if ( iswmode == 1 ) then
                zgam1 = 1.75 - zssa1 * (f_one + zasy3)
                zgam2 =-0.25 + zssa1 * (f_one - zasy3)
                zgam3 = 0.5  - zasy3 * cosz
              elseif ( iswmode == 2 ) then               ! pifm
                zgam1 = 2.0 - zssa1 * (1.25 + zasy3)
                zgam2 = 0.75* zssa1 * (f_one- zasy1)
                zgam3 = 0.5 - zasy3 * cosz
              elseif ( iswmode == 3 ) then               ! discrete ordinates
                zgam1 = zsr3 * (2.0 - zssa1 * (1.0 + zasy1)) * 0.5
                zgam2 = zsr3 * zssa1 * (1.0 - zasy1) * 0.5
                zgam3 = (1.0 - zsr3 * zasy1 * cosz) * 0.5
              endif
              zgam4 = f_one - zgam3

!>  - Compute homogeneous reflectance and transmittance for both convertive 
!! and non-convertive scattering.

              if ( zssaw >= zcrit ) then    ! for conservative scattering
                za1 = zgam1 * cosz - zgam3
                za2 = zgam1 * ztau1

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

                zb1 = min ( ztau1*sntz , 500.0 )
                if ( zb1 <= od_lo ) then
                  zb2 = f_one - zb1 + 0.5*zb1*zb1
                else
                  ftind = zb1 / (bpade + zb1)
                  itind = ftind*NTBMX + 0.5
                  zb2 = exp_tbl(itind)
                endif

!      ...  collimated beam
                zrefb(kp) = max(f_zero, min(f_one,                      &
     &                      (za2 - za1*(f_one - zb2))/(f_one + za2) ))
                ztrab(kp) = max(f_zero, min(f_one, f_one-zrefb(kp)))

!      ...  isotropic incidence
                zrefd(kp) = max(f_zero, min(f_one, za2 / (f_one+za2) ))
                ztrad(kp) = max(f_zero, min(f_one, f_one - zrefd(kp) ))

              else                          ! for non-conservative scattering
                za1 = zgam1*zgam4 + zgam2*zgam3
                za2 = zgam1*zgam3 + zgam2*zgam4
                zrk = sqrt ( (zgam1 - zgam2) * (zgam1 + zgam2) )
                zrk2= 2.0 * zrk

                zrp  = zrk * cosz
                zrp1 = f_one + zrp
                zrm1 = f_one - zrp
                zrpp1= f_one - zrp*zrp
                zrpp = sign( max(flimit, abs(zrpp1)), zrpp1 )    ! avoid numerical singularity
                zrkg1= zrk + zgam1
                zrkg3= zrk * zgam3
                zrkg4= zrk * zgam4

                zr1  = zrm1 * (za2 + zrkg3)
                zr2  = zrp1 * (za2 - zrkg3)
                zr3  = zrk2 * (zgam3 - za2*cosz)
                zr4  = zrpp * zrkg1
                zr5  = zrpp * (zrk - zgam1)

                zt1  = zrp1 * (za1 + zrkg4)
                zt2  = zrm1 * (za1 - zrkg4)
                zt3  = zrk2 * (zgam4 + za1*cosz)

!  --- ...  use exponential lookup table for transmittance, or expansion
!           of exponential for low optical depth

                zb1 = min ( zrk*ztau1, 500.0 )
                if ( zb1 <= od_lo ) then
                  zexm1 = f_one - zb1 + 0.5*zb1*zb1
                else
                  ftind = zb1 / (bpade + zb1)
                  itind = ftind*NTBMX + 0.5
                  zexm1 = exp_tbl(itind)
                endif
                zexp1 = f_one / zexm1

                zb2 = min ( ztau1*sntz, 500.0 )
                if ( zb2 <= od_lo ) then
                  zexm2 = f_one - zb2 + 0.5*zb2*zb2
                else
                  ftind = zb2 / (bpade + zb2)
                  itind = ftind*NTBMX + 0.5
                  zexm2 = exp_tbl(itind)
                endif
                zexp2 = f_one / zexm2
                ze1r45 = zr4*zexp1 + zr5*zexm1

!      ...  collimated beam
                if ( ze1r45>=-eps1 .and. ze1r45<=eps1 ) then
                  zrefb(kp) = eps1
                  ztrab(kp) = zexm2
                else
                  zden1 = zssa1 / ze1r45
                  zrefb(kp) = max(f_zero, min(f_one,                    &
     &                        (zr1*zexp1-zr2*zexm1-zr3*zexm2)*zden1 ))
                  ztrab(kp) = max(f_zero, min(f_one, zexm2*(f_one -     &
     &                        (zt1*zexp1-zt2*zexm1-zt3*zexp2)*zden1) ))
                endif

!      ...  diffuse beam
                zden1 = zr4 / (ze1r45 * zrkg1)
                zrefd(kp) = max(f_zero, min(f_one,                      &
     &                      zgam2*(zexp1 - zexm1)*zden1 ))
                ztrad(kp) = max(f_zero, min(f_one, zrk2*zden1 ))
              endif    ! end if_zssaw_block

!  --- ...  direct beam transmittance. use exponential lookup table
!           for transmittance, or expansion of exponential for low
!           optical depth

              zr1 = ztau1 * sntz
              if ( zr1 <= od_lo ) then
                zexp3 = f_one - zr1 + 0.5*zr1*zr1
              else
                ftind = zr1 / (bpade + zr1)
                itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
                zexp3 = exp_tbl(itind)
              endif

              zldbt(kp) = zexp3
              ztdbt(k)  = zexp3 * ztdbt(kp)

!  --- ...  pre-delta-scaling clear and cloudy direct beam transmittance
!           (must use 'orig', unscaled cloud optical depth)

              zr1 = ztau0 * sntz
              if ( zr1 <= od_lo ) then
                zexp4 = f_one - zr1 + 0.5*zr1*zr1
              else
                ftind = zr1 / (bpade + zr1)
                itind = max(0, min(NTBMX, int(0.5+NTBMX*ftind) ))
                zexp4 = exp_tbl(itind)
              endif

              ztdbt0 = zexp4 * ztdbt0

            else     ! if_cldfmc_block  ---  it is a clear layer

!  --- ...  direct beam transmittance
              ztdbt(k) = zldbt(kp) * ztdbt(kp)

!>  -  Calculate pre-delta-scaling clear and cloudy direct beam transmittance.
              ztdbt0 = zldbt0(k) * ztdbt0

            endif    ! end if_cldfmc_block
          enddo   ! end do_k_loop

!> -# Call vrtqdr(), to  perform vertical quadrature

          call vrtqdr                                                   &
!  ---  inputs:
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       &
     &       nlay, nlp1,                                                &
!  ---  outputs:
     &       zfu, zfd                                                   &
     &     )

!  --- ...  compute upward and downward fluxes at levels
          do k = 1, nlp1
            fxupc(k,ib) = fxupc(k,ib) + zsolar*zfu(k)
            fxdnc(k,ib) = fxdnc(k,ib) + zsolar*zfd(k)
          enddo

!> -# Process and save outputs.
!!  - surface downward beam/diffused flux components
          zb1 = zsolar*ztdbt0
          zb2 = zsolar*(zfd(1) - ztdbt0)

          if (ibd /= 0) then
           sfbmc(ibd) = sfbmc(ibd) + zb1
           sfdfc(ibd) = sfdfc(ibd) + zb2
          else
            zf1 = 0.5 * zb1
            zf2 = 0.5 * zb2
            sfbmc(1) = sfbmc(1) + zf1
            sfdfc(1) = sfdfc(1) + zf2
            sfbmc(2) = sfbmc(2) + zf1
            sfdfc(2) = sfdfc(2) + zf2
          endif
!         sfbmc(ibd) = sfbmc(ibd) + zsolar*ztdbt0
!         sfdfc(ibd) = sfdfc(ibd) + zsolar*(zfd(1) - ztdbt0)

        endif      ! end if_cf1_block

      enddo  lab_do_jg

!  --- ...  end of g-point loop

      do ib = 1, nbdsw
        ftoadc = ftoadc + fxdn0(nlp1,ib)
        ftoau0 = ftoau0 + fxup0(nlp1,ib)
        fsfcu0 = fsfcu0 + fxup0(1,ib)
        fsfcd0 = fsfcd0 + fxdn0(1,ib)
      enddo

!>  - uv-b surface downward flux
      ibd = nuvb - nblow + 1
      suvbf0 = fxdn0(1,ibd)

      if ( cf1 <= eps ) then       ! clear column, set total-sky=clear-sky fluxes
        do ib = 1, nbdsw
          do k = 1, nlp1
            fxupc(k,ib) = fxup0(k,ib)
            fxdnc(k,ib) = fxdn0(k,ib)
          enddo
        enddo

        ftoauc = ftoau0
        fsfcuc = fsfcu0
        fsfcdc = fsfcd0

!>  - surface downward beam/diffused flux components
        sfbmc(1) = sfbm0(1)
        sfdfc(1) = sfdf0(1)
        sfbmc(2) = sfbm0(2)
        sfdfc(2) = sfdf0(2)

!>  -  uv-b surface downward flux
        suvbfc = suvbf0
      else                        ! cloudy column, compute total-sky fluxes
        do ib = 1, nbdsw
          ftoauc = ftoauc + fxupc(nlp1,ib)
          fsfcuc = fsfcuc + fxupc(1,ib)
          fsfcdc = fsfcdc + fxdnc(1,ib)
        enddo

!! --- ...  uv-b surface downward flux
        suvbfc = fxdnc(1,ibd)
      endif    ! end if_cf1_block

      return
!...................................
      end subroutine spcvrtm
!! @}
!-----------------------------------

!>\ingroup module_radsw_main
!> This subroutine is called by spcvrtc() and spcvrtm(), and computes
!! the upward and downward radiation fluxes.
!!\param zrefb           layer direct beam reflectivity
!!\param zrefd           layer diffuse reflectivity
!!\param ztrab           layer direct beam transmissivity
!!\param ztrad           layer diffuse transmissivity
!!\param zldbt           layer mean beam transmittance
!!\param ztdbt           total beam transmittance at levels
!!\param NLAY, NLP1      number of layers/levels
!!\param zfu             upward flux at layer interface
!!\param zfd             downward flux at layer interface
!!\section General_vrtqdr vrtqdr General Algorithm
!> @{
!-----------------------------------
      subroutine vrtqdr                                                 &
     &     ( zrefb,zrefd,ztrab,ztrad,zldbt,ztdbt,                       & ! inputs
     &       NLAY, NLP1,                                                &
     &       zfu, zfd                                                   & ! outputs:
     &     )

!  ===================  program usage description  ===================  !
!                                                                       !
!   purpose:  computes the upward and downward radiation fluxes         !
!                                                                       !
!   interface:  "vrtqdr" is called by "spcvrc" and "spcvrm"             !
!                                                                       !
!   subroutines called : none                                           !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  input variables:                                                     !
!    zrefb(NLP1)     - layer direct beam reflectivity                   !
!    zrefd(NLP1)     - layer diffuse reflectivity                       !
!    ztrab(NLP1)     - layer direct beam transmissivity                 !
!    ztrad(NLP1)     - layer diffuse transmissivity                     !
!    zldbt(NLP1)     - layer mean beam transmittance                    !
!    ztdbt(NLP1)     - total beam transmittance at levels               !
!    NLAY, NLP1      - number of layers/levels                          !
!                                                                       !
!  output variables:                                                    !
!    zfu  (NLP1)     - upward flux at layer interface                   !
!    zfd  (NLP1)     - downward flux at layer interface                 !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer, intent(in) :: nlay, nlp1

      real (kind=kind_phys), dimension(nlp1), intent(in) :: zrefb,      &
     &       zrefd, ztrab, ztrad, ztdbt, zldbt

!  ---  outputs:
      real (kind=kind_phys), dimension(nlp1), intent(out) :: zfu, zfd

!  ---  locals:
      real (kind=kind_phys), dimension(nlp1) :: zrupb,zrupd,zrdnd,ztdn

      real (kind=kind_phys) :: zden1

      integer :: k, kp
!
!===> ... begin here
!

!> -# Link lowest layer with surface.
        zrupb(1) = zrefb(1)        ! direct beam
        zrupd(1) = zrefd(1)        ! diffused

!> -# Pass from bottom to top.
        do k = 1, nlay
          kp = k + 1

          zden1 = f_one / ( f_one - zrupd(k)*zrefd(kp) )
          zrupb(kp) = zrefb(kp) + ( ztrad(kp) *                         &
     &                ( (ztrab(kp) - zldbt(kp))*zrupd(k) +              &
     &                zldbt(kp)*zrupb(k)) ) * zden1
          zrupd(kp) = zrefd(kp) + ztrad(kp)*ztrad(kp)*zrupd(k)*zden1
        enddo

!> -# Upper boundary conditions
        ztdn (nlp1) = f_one
        zrdnd(nlp1) = f_zero
        ztdn (nlay) = ztrab(nlp1)
        zrdnd(nlay) = zrefd(nlp1)

!> -# Pass from top to bottom
        do k = nlay, 2, -1
          zden1 = f_one / (f_one - zrefd(k)*zrdnd(k))
          ztdn (k-1) = ztdbt(k)*ztrab(k) + ( ztrad(k) *                 &
     &                 ( (ztdn(k) - ztdbt(k)) + ztdbt(k) *              &
     &                 zrefb(k)*zrdnd(k) )) * zden1
          zrdnd(k-1) = zrefd(k) + ztrad(k)*ztrad(k)*zrdnd(k)*zden1
        enddo

!> -# Up and down-welling fluxes at levels.
        do k = 1, nlp1
          zden1 = f_one / (f_one - zrdnd(k)*zrupd(k))
          zfu(k) = ( ztdbt(k)*zrupb(k) +                                &
     &             (ztdn(k) - ztdbt(k))*zrupd(k) ) * zden1
          zfd(k) = ztdbt(k) + ( ztdn(k) - ztdbt(k) +                    &
     &             ztdbt(k)*zrupb(k)*zrdnd(k) ) * zden1
        enddo

      return
!...................................
      end subroutine vrtqdr
!-----------------------------------
!> @}

!>\ingroup module_radsw_main
!> This subroutine calculates optical depths for gaseous absorption and
!! rayleigh scattering
!!\n subroutine called taumol## (## = 16-29)
!!\param colamt           column amounts of absorbing gases the index
!!                        are for h2o, co2, o3, n2o, ch4, and o2,
!!                        respectively \f$(mol/cm^2)\f$
!!\param colmol           total column amount (dry air+water vapor)
!!\param facij            for each layer, these are factors that are
!!                        needed to compute the interpolation factors
!!                        that multiply the appropriate reference
!!                        k-values.  a value of 0/1 for i,j indicates
!!                        that the corresponding factor multiplies
!!                        reference k-value for the lower/higher of the
!!                        two appropriate temperatures, and altitudes,
!!                        respectively.
!!\param jp               the index of the lower (in altitude) of the
!!                        two appropriate ref pressure levels needed
!!                        for interpolation.
!!\param jt, jt1          the indices of the lower of the two approp
!!                        ref temperatures needed for interpolation
!!                        (for pressure levels jp and jp+1, respectively)
!!\param laytrop          tropopause layer index
!!\param forfac           scale factor needed to foreign-continuum.
!!\param forfrac          factor needed for temperature interpolation
!!\param indfor           index of the lower of the two appropriate
!!                        reference temperatures needed for
!!                        foreign-continuum interpolation
!!\param selffac          scale factor needed to h2o self-continuum.
!!\param selffrac         factor needed for temperature interpolation
!!                        of reference h2o self-continuum data
!!\param indself          index of the lower of the two appropriate
!!                        reference temperatures needed for the
!!                        self-continuum interpolation
!!\param nlay             number of vertical layers
!!\param sfluxzen         spectral distribution of incoming solar flux
!!\param taug             spectral optical depth for gases
!!\param taur             opt depth for rayleigh scattering
!>\section gen_al_taumol taumol General Algorithm
!! @{
!-----------------------------------
      subroutine taumol                                                 &
     &     ( colamt,colmol,fac00,fac01,fac10,fac11,jp,jt,jt1,laytrop,   & !  ---  inputs
     &       forfac,forfrac,indfor,selffac,selffrac,indself, nlay,      &
     &       sfluxzen, taug, taur                                       &  !  ---  outputs
     &     )

!  ==================   program usage description   ==================  !
!                                                                       !
!  description:                                                         !
!    calculate optical depths for gaseous absorption and rayleigh       !
!    scattering.                                                        !
!                                                                       !
!  subroutines called: taugb## (## = 16 - 29)                           !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                         size !
!    colamt  - real, column amounts of absorbing gases the index        !
!                    are for h2o, co2, o3, n2o, ch4, and o2,            !
!                    respectively (molecules/cm**2)          nlay*maxgas!
!    colmol  - real, total column amount (dry air+water vapor)     nlay !
!    facij   - real, for each layer, these are factors that are         !
!                    needed to compute the interpolation factors        !
!                    that multiply the appropriate reference k-         !
!                    values.  a value of 0/1 for i,j indicates          !
!                    that the corresponding factor multiplies           !
!                    reference k-value for the lower/higher of the      !
!                    two appropriate temperatures, and altitudes,       !
!                    respectively.                                 naly !
!    jp      - real, the index of the lower (in altitude) of the        !
!                    two appropriate ref pressure levels needed         !
!                    for interpolation.                            nlay !
!    jt, jt1 - integer, the indices of the lower of the two approp      !
!                    ref temperatures needed for interpolation (for     !
!                    pressure levels jp and jp+1, respectively)    nlay !
!    laytrop - integer, tropopause layer index                       1  !
!    forfac  - real, scale factor needed to foreign-continuum.     nlay !
!    forfrac - real, factor needed for temperature interpolation   nlay !
!    indfor  - integer, index of the lower of the two appropriate       !
!                    reference temperatures needed for foreign-         !
!                    continuum interpolation                       nlay !
!    selffac - real, scale factor needed to h2o self-continuum.    nlay !
!    selffrac- real, factor needed for temperature interpolation        !
!                    of reference h2o self-continuum data          nlay !
!    indself - integer, index of the lower of the two appropriate       !
!                    reference temperatures needed for the self-        !
!                    continuum interpolation                       nlay !
!    nlay    - integer, number of vertical layers                    1  !
!                                                                       !
!  output:                                                              !
!    sfluxzen- real, spectral distribution of incoming solar flux ngptsw!
!    taug    - real, spectral optical depth for gases        nlay*ngptsw!
!    taur    - real, opt depth for rayleigh scattering       nlay*ngptsw!
!                                                                       !
!  ===================================================================  !
!  ************     original subprogram description    ***************  !
!                                                                       !
!                  optical depths developed for the                     !
!                                                                       !
!                rapid radiative transfer model (rrtm)                  !
!                                                                       !
!            atmospheric and environmental research, inc.               !
!                        131 hartwell avenue                            !
!                        lexington, ma 02421                            !
!                                                                       !
!                                                                       !
!                           eli j. mlawer                               !
!                         jennifer delamere                             !
!                         steven j. taubman                             !
!                         shepard a. clough                             !
!                                                                       !
!                                                                       !
!                                                                       !
!                       email:  mlawer@aer.com                          !
!                       email:  jdelamer@aer.com                        !
!                                                                       !
!        the authors wish to acknowledge the contributions of the       !
!        following people:  patrick d. brown, michael j. iacono,        !
!        ronald e. farren, luke chen, robert bergstrom.                 !
!                                                                       !
!  *******************************************************************  !
!                                                                       !
!  taumol                                                               !
!                                                                       !
!    this file contains the subroutines taugbn (where n goes from       !
!    16 to 29).  taugbn calculates the optical depths and Planck        !
!    fractions per g-value and layer for band n.                        !
!                                                                       !
!  output:  optical depths (unitless)                                   !
!           fractions needed to compute planck functions at every layer !
!           and g-value                                                 !
!                                                                       !
!  modifications:                                                       !
!                                                                       !
! revised: adapted to f90 coding, j.-j.morcrette, ecmwf, feb 2003       !
! revised: modified for g-point reduction, mjiacono, aer, dec 2003      !
! revised: reformatted for consistency with rrtmg_lw, mjiacono, aer,    !
!          jul 2006                                                     !
!                                                                       !
!  *******************************************************************  !
!  ======================  end of description block  =================  !

!  ---  inputs:
      integer,               intent(in) :: nlay, laytrop

      integer, dimension(nlay), intent(in) :: indfor, indself,          &
     &       jp, jt, jt1

      real (kind=kind_phys), dimension(nlay),  intent(in) :: colmol,    &
     &       fac00, fac01, fac10, fac11, forfac, forfrac, selffac,      &
     &       selffrac

      real (kind=kind_phys), dimension(nlay,maxgas),intent(in) :: colamt

!  ---  outputs:
      real (kind=kind_phys), dimension(ngptsw), intent(out) :: sfluxzen

      real (kind=kind_phys), dimension(nlay,ngptsw), intent(out) ::     &
     &       taug, taur

!  ---  locals:
      real (kind=kind_phys) :: fs, speccomb, specmult, colm1, colm2

      integer, dimension(nlay,nblow:nbhgh) :: id0, id1

      integer :: ibd, j, jb, js, k, klow, khgh, klim, ks, njb, ns
!
!===> ... begin here
!
!  --- ...  loop over each spectral band

      do jb = nblow, nbhgh

!  --- ...  indices for layer optical depth

        do k = 1, laytrop
          id0(k,jb) = ((jp(k)-1)*5 + (jt (k)-1)) * nspa(jb)
          id1(k,jb) = ( jp(k)   *5 + (jt1(k)-1)) * nspa(jb)
        enddo

        do k = laytrop+1, nlay
          id0(k,jb) = ((jp(k)-13)*5 + (jt (k)-1)) * nspb(jb)
          id1(k,jb) = ((jp(k)-12)*5 + (jt1(k)-1)) * nspb(jb)
        enddo

!  --- ...  calculate spectral flux at toa

        ibd = ibx(jb)
        njb = ng (jb)
        ns  = ngs(jb)

        select case (jb)

          case (16, 20, 23, 25, 26, 29)

            do j = 1, njb
              sfluxzen(ns+j) = sfluxref01(j,1,ibd)
            enddo

          case (27)

            do j = 1, njb
              sfluxzen(ns+j) = scalekur * sfluxref01(j,1,ibd)
            enddo

          case default

            if (jb==17 .or. jb==28) then

              ks = nlay
              lab_do_k1 : do k = laytrop, nlay-1
                if (jp(k)<layreffr(jb) .and. jp(k+1)>=layreffr(jb)) then
                  ks = k + 1
                  exit lab_do_k1
                endif
              enddo  lab_do_k1

              colm1 = colamt(ks,ix1(jb))
              colm2 = colamt(ks,ix2(jb))
              speccomb = colm1 + strrat(jb)*colm2
              specmult = specwt(jb) * min( oneminus, colm1/speccomb )
              js = 1 + int( specmult )
              fs = mod(specmult, f_one)

              do j = 1, njb
                sfluxzen(ns+j) = sfluxref02(j,js,ibd)                   &
     &           + fs * (sfluxref02(j,js+1,ibd) - sfluxref02(j,js,ibd))
              enddo

            else

              ks = laytrop
              lab_do_k2 : do k = 1, laytrop-1
                if (jp(k)<layreffr(jb) .and. jp(k+1)>=layreffr(jb)) then
                  ks = k + 1
                  exit lab_do_k2
                endif
              enddo  lab_do_k2

              colm1 = colamt(ks,ix1(jb))
              colm2 = colamt(ks,ix2(jb))
              speccomb = colm1 + strrat(jb)*colm2
              specmult = specwt(jb) * min( oneminus, colm1/speccomb )
              js = 1 + int( specmult )
              fs = mod(specmult, f_one)

              do j = 1, njb
                sfluxzen(ns+j) = sfluxref03(j,js,ibd)                   &
     &           + fs * (sfluxref03(j,js+1,ibd) - sfluxref03(j,js,ibd))
              enddo

            endif

        end select

      enddo

!> - Call taumol## (##: 16-29) to calculate layer optical depth.

!>  - call taumol16()
      call taumol16
!>  - call taumol17()
      call taumol17
!>  - call taumol18()
      call taumol18
!>  - call taumol19()
      call taumol19
!>  - call taumol20()
      call taumol20
!>  - call taumol21()
      call taumol21
!>  - call taumol22()
      call taumol22
!>  - call taumol23()
      call taumol23
!>  - call taumol24()
      call taumol24
!>  - call taumol25()
      call taumol25
!>  - call taumol26()
      call taumol26
!>  - call taumol27()
      call taumol27
!>  - call taumol28()
      call taumol28
!>  - call taumol29()
      call taumol29


! =================
      contains
! =================

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 16:  2600-3250
!! cm-1 (low - h2o,ch4; high - ch4)
!-----------------------------------
      subroutine taumol16
!...................................

!  ------------------------------------------------------------------  !
!     band 16:  2600-3250 cm-1 (low - h2o,ch4; high - ch4)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb16

!  ---  locals:

      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ... compute the optical depth by interpolating in ln(pressure),
!          temperature, and appropriate species.  below laytrop, the water
!          vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG16
          taur(k,NS16+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        speccomb = colamt(k,1) + strrat(16)*colamt(k,5)
        specmult = 8.0 * min( oneminus, colamt(k,1)/speccomb )

        js = 1 + int( specmult )
        fs = mod( specmult, f_one )
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,16) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 9
        ind04 = ind01 + 10
        ind11 = id1(k,16) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 9
        ind14 = ind11 + 10
        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG16
          taug(k,NS16+j) = speccomb                                     &
     &        *( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)        &
     &        +  fac010 * absa(ind03,j) + fac110 * absa(ind04,j)        &
     &        +  fac001 * absa(ind11,j) + fac101 * absa(ind12,j)        &
     &        +  fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )      &
     &        + colamt(k,1) * (selffac(k) * (selfref(inds,j)            &
     &        + selffrac(k) * (selfref(indsp,j)-selfref(inds,j)))       &
     &        + forfac(k) * (forref(indf,j) + forfrac(k)                &
     &        * (forref(indfp,j) - forref(indf,j))))
        enddo
      enddo

      do k = laytrop+1, nlay
        ind01 = id0(k,16) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,16) + 1
        ind12 = ind11 + 1

        do j = 1, NG16
          taug(k,NS16+j) = colamt(k,5)                                  &
     &      * ( fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)         &
     &      +   fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j) )
        enddo
      enddo

      return
!...................................
      end subroutine taumol16
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 17:  3250-4000
!! cm-1 (low - h2o,co2; high - h2o,co2)
!-----------------------------------
      subroutine taumol17
!...................................

!  ------------------------------------------------------------------  !
!     band 17:  3250-4000 cm-1 (low - h2o,co2; high - h2o,co2)         !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb17

!  ---  locals:
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG17
          taur(k,NS17+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        speccomb = colamt(k,1) + strrat(17)*colamt(k,2)
        specmult = 8.0 * min(oneminus, colamt(k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,17) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 9
        ind04 = ind01 + 10
        ind11 = id1(k,17) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 9
        ind14 = ind11 + 10

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG17
          taug(k,NS17+j) = speccomb                                     &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &        + colamt(k,1) * (selffac(k) * (selfref(inds,j)            &
     &        + selffrac(k) * (selfref(indsp,j)-selfref(inds,j)))       &
     &        + forfac(k) * (forref(indf,j) + forfrac(k)                &
     &        * (forref(indfp,j) - forref(indf,j))))
        enddo
      enddo

      do k = laytrop+1, nlay
        speccomb = colamt(k,1) + strrat(17)*colamt(k,2)
        specmult = 4.0 * min(oneminus, colamt(k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,17) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 5
        ind04 = ind01 + 6
        ind11 = id1(k,17) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 5
        ind14 = ind11 + 6

        indf = indfor(k)
        indfp= indf + 1

        do j = 1, NG17
          taug(k,NS17+j) = speccomb                                     &
     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
     &        + colamt(k,1) * forfac(k) * (forref(indf,j)               &
     &        + forfrac(k) * (forref(indfp,j) - forref(indf,j)))
        enddo
      enddo

      return
!...................................
      end subroutine taumol17
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 18:  4000-4650
!! cm-1 (low - h2o,ch4; high - ch4)
!-----------------------------------
      subroutine taumol18
!...................................

!  ------------------------------------------------------------------  !
!     band 18:  4000-4650 cm-1 (low - h2o,ch4; high - ch4)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb18

!  ---  locals:
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG18
          taur(k,NS18+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        speccomb = colamt(k,1) + strrat(18)*colamt(k,5)
        specmult = 8.0 * min(oneminus, colamt(k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,18) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 9
        ind04 = ind01 + 10
        ind11 = id1(k,18) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 9
        ind14 = ind11 + 10

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG18
          taug(k,NS18+j) = speccomb                                     &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &        + colamt(k,1) * (selffac(k) * (selfref(inds,j)            &
     &        + selffrac(k) * (selfref(indsp,j)-selfref(inds,j)))       &
     &        + forfac(k) * (forref(indf,j) + forfrac(k)                &
     &        * (forref(indfp,j) - forref(indf,j))))
        enddo
      enddo

      do k = laytrop+1, nlay
        ind01 = id0(k,18) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,18) + 1
        ind12 = ind11 + 1

        do j = 1, NG18
          taug(k,NS18+j) = colamt(k,5)                                  &
     &        * ( fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &        +   fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j) )
        enddo
      enddo

      return
!...................................
      end subroutine taumol18
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 19:  4650-5150
!! cm-1 (low - h2o,co2; high - co2)
!-----------------------------------
      subroutine taumol19
!...................................

!  ------------------------------------------------------------------  !
!     band 19:  4650-5150 cm-1 (low - h2o,co2; high - co2)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb19

!  ---  locals:
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG19
          taur(k,NS19+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        speccomb = colamt(k,1) + strrat(19)*colamt(k,2)
        specmult = 8.0 * min(oneminus, colamt(k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,19) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 9
        ind04 = ind01 + 10
        ind11 = id1(k,19) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 9
        ind14 = ind11 + 10

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG19
          taug(k,NS19+j) = speccomb                                     &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &        + colamt(k,1) * (selffac(k) * (selfref(inds,j)            &
     &        + selffrac(k) * (selfref(indsp,j)-selfref(inds,j)))       &
     &        + forfac(k) * (forref(indf,j) + forfrac(k)                &
     &        * (forref(indfp,j) - forref(indf,j))))
        enddo
      enddo

      do k = laytrop+1, nlay
        ind01 = id0(k,19) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,19) + 1
        ind12 = ind11 + 1

        do j = 1, NG19
          taug(k,NS19+j) = colamt(k,2)                                  &
     &        * ( fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &        +   fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j) )
        enddo
      enddo

!...................................
      end subroutine taumol19
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 20:  5150-6150
!! cm-1 (low - h2o; high - h2o)
!-----------------------------------
      subroutine taumol20
!...................................

!  ------------------------------------------------------------------  !
!     band 20:  5150-6150 cm-1 (low - h2o; high - h2o)                 !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb20

!  ---  locals:
      real (kind=kind_phys) :: tauray

      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG20
          taur(k,NS20+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        ind01 = id0(k,20) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,20) + 1
        ind12 = ind11 + 1

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG20
          taug(k,NS20+j) = colamt(k,1)                                  &
     &        * ( (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)      &
     &        +    fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j))     &
     &        +   selffac(k) * (selfref(inds,j) + selffrac(k)           &
     &        *   (selfref(indsp,j) - selfref(inds,j)))                 &
     &        +   forfac(k) * (forref(indf,j) + forfrac(k)              &
     &        *   (forref(indfp,j) - forref(indf,j))) )                 &
     &        + colamt(k,5) * absch4(j)
        enddo
      enddo

      do k = laytrop+1, nlay
        ind01 = id0(k,20) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,20) + 1
        ind12 = ind11 + 1

        indf = indfor(k)
        indfp= indf + 1

        do j = 1, NG20
          taug(k,NS20+j) = colamt(k,1)                                  &
     &        * ( fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &        +   fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j)       &
     &        +   forfac(k) * (forref(indf,j) + forfrac(k)              &
     &        *   (forref(indfp,j) - forref(indf,j))) )                 &
     &        + colamt(k,5) * absch4(j)
        enddo
      enddo

      return
!...................................
      end subroutine taumol20
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 21:  6150-7700
!! cm-1 (low - h2o,co2; high - h2o,co2)
!-----------------------------------
      subroutine taumol21
!...................................

!  ------------------------------------------------------------------  !
!     band 21:  6150-7700 cm-1 (low - h2o,co2; high - h2o,co2)         !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb21

!  ---  locals:
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG21
          taur(k,NS21+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        speccomb = colamt(k,1) + strrat(21)*colamt(k,2)
        specmult = 8.0 * min(oneminus, colamt(k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,21) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 9
        ind04 = ind01 + 10
        ind11 = id1(k,21) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 9
        ind14 = ind11 + 10

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG21
          taug(k,NS21+j) = speccomb                                     &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &        + colamt(k,1) * (selffac(k) * (selfref(inds,j)            &
     &        + selffrac(k) * (selfref(indsp,j) - selfref(inds,j)))     &
     &        + forfac(k) * (forref(indf,j) + forfrac(k)                &
     &        * (forref(indfp,j) - forref(indf,j))))
        enddo
      enddo

      do k = laytrop+1, nlay
        speccomb = colamt(k,1) + strrat(21)*colamt(k,2)
        specmult = 4.0 * min(oneminus, colamt(k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,21) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 5
        ind04 = ind01 + 6
        ind11 = id1(k,21) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 5
        ind14 = ind11 + 6

        indf = indfor(k)
        indfp= indf + 1

        do j = 1, NG21
          taug(k,NS21+j) = speccomb                                     &
     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )     &
     &        + colamt(k,1) * forfac(k) * (forref(indf,j)               &
     &        + forfrac(k) * (forref(indfp,j) - forref(indf,j)))
        enddo
      enddo

!...................................
      end subroutine taumol21
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 22:  7700-8050
!! cm-1 (low - h2o,o2; high - o2)
!-----------------------------------
      subroutine taumol22
!...................................

!  ------------------------------------------------------------------  !
!     band 22:  7700-8050 cm-1 (low - h2o,o2; high - o2)               !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb22

!  ---  locals:
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111,  &
     &       o2adj, o2cont, o2tem

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!
!  --- ...  the following factor is the ratio of total o2 band intensity (lines
!           and mate continuum) to o2 band intensity (line only). it is needed
!           to adjust the optical depths since the k's include only lines.

      o2adj = 1.6
      o2tem = 4.35e-4 / (350.0*2.0)


!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG22
          taur(k,NS22+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        o2cont   = o2tem * colamt(k,6)
        speccomb = colamt(k,1) + strrat(22)*colamt(k,6)
        specmult = 8.0 * min(oneminus, colamt(k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,22) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 9
        ind04 = ind01 + 10
        ind11 = id1(k,22) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 9
        ind14 = ind11 + 10

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG22
          taug(k,NS22+j) = speccomb                                     &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &        + colamt(k,1) * (selffac(k) * (selfref(inds,j)            &
     &        + selffrac(k) * (selfref(indsp,j)-selfref(inds,j)))       &
     &        + forfac(k) * (forref(indf,j) + forfrac(k)                &
     &        * (forref(indfp,j) - forref(indf,j)))) + o2cont
        enddo
      enddo

      do k = laytrop+1, nlay
        o2cont = o2tem * colamt(k,6)

        ind01 = id0(k,22) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,22) + 1
        ind12 = ind11 + 1

        do j = 1, NG22
          taug(k,NS22+j) = colamt(k,6) * o2adj                          &
     &        * ( fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &        +   fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j) )     &
     &        + o2cont
        enddo
      enddo

      return
!...................................
      end subroutine taumol22
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 23:  8050-12850
!! cm-1 (low - h2o; high - nothing)
!-----------------------------------
      subroutine taumol23
!...................................

!  ------------------------------------------------------------------  !
!     band 23:  8050-12850 cm-1 (low - h2o; high - nothing)            !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb23

!  ---  locals:
      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        do j = 1, NG23
          taur(k,NS23+j) = colmol(k) * rayl(j)
        enddo
      enddo

      do k = 1, laytrop
        ind01 = id0(k,23) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,23) + 1
        ind12 = ind11 + 1

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG23
          taug(k,NS23+j) = colamt(k,1) * (givfac                        &
     &        * ( fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &        +   fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j) )     &
     &        + selffac(k) * (selfref(inds,j) + selffrac(k)             &
     &        * (selfref(indsp,j) - selfref(inds,j)))                   &
     &        + forfac(k) * (forref(indf,j) + forfrac(k)                &
     &        * (forref(indfp,j) - forref(indf,j))))
        enddo
      enddo

      do k = laytrop+1, nlay
        do j = 1, NG23
          taug(k,NS23+j) = f_zero
        enddo
      enddo

!...................................
      end subroutine taumol23
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 24:  12850-16000
!! cm-1 (low - h2o,o2; high - o2)
!-----------------------------------
      subroutine taumol24
!...................................

!  ------------------------------------------------------------------  !
!     band 24:  12850-16000 cm-1 (low - h2o,o2; high - o2)             !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb24

!  ---  locals:
      real (kind=kind_phys) :: speccomb, specmult, fs, fs1,             &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: inds, indf, indsp, indfp, j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, laytrop
        speccomb = colamt(k,1) + strrat(24)*colamt(k,6)
        specmult = 8.0 * min(oneminus, colamt(k,1) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,24) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 9
        ind04 = ind01 + 10
        ind11 = id1(k,24) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 9
        ind14 = ind11 + 10

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG24
          taug(k,NS24+j) = speccomb                                     &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )     &
     &        + colamt(k,3) * abso3a(j) + colamt(k,1)                   &
     &        * (selffac(k) * (selfref(inds,j) + selffrac(k)            &
     &        * (selfref(indsp,j) - selfref(inds,j)))                   &
     &        + forfac(k) * (forref(indf,j) + forfrac(k)                &
     &        * (forref(indfp,j) - forref(indf,j))))

          taur(k,NS24+j) = colmol(k)                                    &
     &           * (rayla(j,js) + fs*(rayla(j,js+1) - rayla(j,js)))
        enddo
      enddo

      do k = laytrop+1, nlay
        ind01 = id0(k,24) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,24) + 1
        ind12 = ind11 + 1

        do j = 1, NG24
          taug(k,NS24+j) = colamt(k,6)                                  &
     &        * ( fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &        +   fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j) )     &
     &        + colamt(k,3) * abso3b(j)

          taur(k,NS24+j) = colmol(k) * raylb(j)
        enddo
      enddo

      return
!...................................
      end subroutine taumol24
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 25:  16000-22650
!! cm-1 (low - h2o; high - nothing)
!-----------------------------------
      subroutine taumol25
!...................................

!  ------------------------------------------------------------------  !
!     band 25:  16000-22650 cm-1 (low - h2o; high - nothing)           !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb25

!  ---  locals:
      integer :: ind01, ind02, ind11, ind12
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        do j = 1, NG25
          taur(k,NS25+j) = colmol(k) * rayl(j)
        enddo
      enddo

      do k = 1, laytrop
        ind01 = id0(k,25) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,25) + 1
        ind12 = ind11 + 1

        do j = 1, NG25
          taug(k,NS25+j) = colamt(k,1)                                  &
     &        * ( fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &        +   fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j) )     &
     &        + colamt(k,3) * abso3a(j)
        enddo
      enddo

      do k = laytrop+1, nlay
        do j = 1, NG25
          taug(k,NS25+j) = colamt(k,3) * abso3b(j)
        enddo
      enddo

      return
!...................................
      end subroutine taumol25
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 26:  22650-29000
!! cm-1 (low - nothing; high - nothing)
!-----------------------------------
      subroutine taumol26
!...................................

!  ------------------------------------------------------------------  !
!     band 26:  22650-29000 cm-1 (low - nothing; high - nothing)       !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb26

!  ---  locals:
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        do j = 1, NG26
          taug(k,NS26+j) = f_zero
          taur(k,NS26+j) = colmol(k) * rayl(j)
        enddo
      enddo

      return
!...................................
      end subroutine taumol26
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 27:  29000-38000
!! cm-1 (low - o3; high - o3)
!-----------------------------------
      subroutine taumol27
!...................................

!  ------------------------------------------------------------------  !
!     band 27:  29000-38000 cm-1 (low - o3; high - o3)                 !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb27
!
!  ---  locals:
      integer :: ind01, ind02, ind11, ind12
      integer :: j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        do j = 1, NG27
          taur(k,NS27+j) = colmol(k) * rayl(j)
        enddo
      enddo

      do k = 1, laytrop
        ind01 = id0(k,27) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,27) + 1
        ind12 = ind11 + 1

        do j = 1, NG27
          taug(k,NS27+j) = colamt(k,3)                                  &
     &        * ( fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)       &
     &        +   fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j) )
        enddo
      enddo

      do k = laytrop+1, nlay
        ind01 = id0(k,27) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,27) + 1
        ind12 = ind11 + 1

        do j = 1, NG27
          taug(k,NS27+j) = colamt(k,3)                                  &
     &        * ( fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &        +   fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j) )
        enddo
      enddo

      return
!...................................
      end subroutine taumol27
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 28:  38000-50000
!! cm-1 (low - o3,o2; high - o3,o2)
!-----------------------------------
      subroutine taumol28
!...................................

!  ------------------------------------------------------------------  !
!     band 28:  38000-50000 cm-1 (low - o3,o2; high - o3,o2)           !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb28

!  ---  locals:
      real (kind=kind_phys) :: speccomb, specmult, tauray, fs, fs1,     &
     &       fac000,fac001,fac010,fac011, fac100,fac101,fac110,fac111

      integer :: ind01, ind02, ind03, ind04, ind11, ind12, ind13, ind14
      integer :: j, js, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG28
          taur(k,NS28+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        speccomb = colamt(k,3) + strrat(28)*colamt(k,6)
        specmult = 8.0 * min(oneminus, colamt(k,3) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,28) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 9
        ind04 = ind01 + 10
        ind11 = id1(k,28) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 9
        ind14 = ind11 + 10

        do j = 1, NG28
          taug(k,NS28+j) = speccomb                                     &
     &        * ( fac000 * absa(ind01,j) + fac100 * absa(ind02,j)       &
     &        +   fac010 * absa(ind03,j) + fac110 * absa(ind04,j)       &
     &        +   fac001 * absa(ind11,j) + fac101 * absa(ind12,j)       &
     &        +   fac011 * absa(ind13,j) + fac111 * absa(ind14,j) )
        enddo
      enddo

      do k = laytrop+1, nlay
        speccomb = colamt(k,3) + strrat(28)*colamt(k,6)
        specmult = 4.0 * min(oneminus, colamt(k,3) / speccomb)

        js = 1 + int(specmult)
        fs = mod(specmult, f_one)
        fs1= f_one - fs
        fac000 = fs1 * fac00(k)
        fac010 = fs1 * fac10(k)
        fac100 = fs  * fac00(k)
        fac110 = fs  * fac10(k)
        fac001 = fs1 * fac01(k)
        fac011 = fs1 * fac11(k)
        fac101 = fs  * fac01(k)
        fac111 = fs  * fac11(k)

        ind01 = id0(k,28) + js
        ind02 = ind01 + 1
        ind03 = ind01 + 5
        ind04 = ind01 + 6
        ind11 = id1(k,28) + js
        ind12 = ind11 + 1
        ind13 = ind11 + 5
        ind14 = ind11 + 6

        do j = 1, NG28
          taug(k,NS28+j) = speccomb                                     &
     &        * ( fac000 * absb(ind01,j) + fac100 * absb(ind02,j)       &
     &        +   fac010 * absb(ind03,j) + fac110 * absb(ind04,j)       &
     &        +   fac001 * absb(ind11,j) + fac101 * absb(ind12,j)       &
     &        +   fac011 * absb(ind13,j) + fac111 * absb(ind14,j) )
        enddo
      enddo

      return
!...................................
      end subroutine taumol28
!-----------------------------------

!>\ingroup module_radsw_main
!> The subroutine computes the optical depth in band 29:  820-2600
!! cm-1 (low - h2o; high - co2)
!-----------------------------------
      subroutine taumol29
!...................................

!  ------------------------------------------------------------------  !
!     band 29:  820-2600 cm-1 (low - h2o; high - co2)                  !
!  ------------------------------------------------------------------  !
!
      use module_radsw_kgb29

!  ---  locals:
      real (kind=kind_phys) :: tauray

      integer :: ind01, ind02, ind11, ind12
      integer :: inds, indf, indsp, indfp, j, k

!
!===> ... begin here
!

!  --- ...  compute the optical depth by interpolating in ln(pressure),
!           temperature, and appropriate species.  below laytrop, the water
!           vapor self-continuum is interpolated (in temperature) separately.

      do k = 1, nlay
        tauray = colmol(k) * rayl

        do j = 1, NG29
          taur(k,NS29+j) = tauray
        enddo
      enddo

      do k = 1, laytrop
        ind01 = id0(k,29) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,29) + 1
        ind12 = ind11 + 1

        inds = indself(k)
        indf = indfor (k)
        indsp= inds + 1
        indfp= indf + 1

        do j = 1, NG29
          taug(k,NS29+j) = colamt(k,1)                                  &
     &        * ( (fac00(k)*absa(ind01,j) + fac10(k)*absa(ind02,j)      &
     &        +    fac01(k)*absa(ind11,j) + fac11(k)*absa(ind12,j) )    &
     &        +  selffac(k) * (selfref(inds,j) + selffrac(k)            &
     &        *  (selfref(indsp,j) - selfref(inds,j)))                  &
     &        +  forfac(k) * (forref(indf,j) + forfrac(k)               &
     &        *  (forref(indfp,j) - forref(indf,j))))                   &
     &        +  colamt(k,2) * absco2(j)
        enddo
      enddo

      do k = laytrop+1, nlay
        ind01 = id0(k,29) + 1
        ind02 = ind01 + 1
        ind11 = id1(k,29) + 1
        ind12 = ind11 + 1

        do j = 1, NG29
          taug(k,NS29+j) = colamt(k,2)                                  &
     &        * ( fac00(k)*absb(ind01,j) + fac10(k)*absb(ind02,j)       &
     &        +   fac01(k)*absb(ind11,j) + fac11(k)*absb(ind12,j) )     &
     &        + colamt(k,1) * absh2o(j)
        enddo
      enddo

      return
!...................................
      end subroutine taumol29
!-----------------------------------

!...................................
      end subroutine taumol
!-----------------------------------
!! @}

!
!........................................!
      end module rrtmg_sw                !
!========================================!
