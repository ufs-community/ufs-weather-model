!>  \file radiation_clouds.f
!!  This file contains routines to compute cloud related quantities
!!  for radiation computations.
!              module_radiation_clouds description             !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!    the 'radiation_clouds.f' contains:                                !
!                                                                      !
!       'module_radiation_clouds' ---  compute cloud related quantities!
!                for radiation computations                            !
!                                                                      !
!    the following are the externally accessable subroutines:          !
!                                                                      !
!       'cld_init'           --- initialization routine                !
!          inputs:                                                     !
!           ( si, NLAY, me )                                           !
!          outputs:                                                    !
!           ( none )                                                   !
!                                                                      !
!       'progcld1'           --- zhao/moorthi prognostic cloud scheme  !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                   !
!            xlat,xlon,slmsk,                                          !
!            IX, NLAY, NLP1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!       'progcld2'           --- ferrier prognostic cloud microphysics !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                   !
!            xlat,xlon,slmsk, f_ice,f_rain,r_rime,flgmin,              !
!            IX, NLAY, NLP1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!       'progcld3'           --- zhao/moorthi prognostic cloud + pdfcld!
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw, cnvw,cnvc,        !
!            xlat,xlon,slmsk,                                          !
!            ix, nlay, nlp1,                                           !
!            deltaq,sup,kdt,me,                                        !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)
!                                                                      !
!       'progclduni'           --- for unified clouds with MG microphys!
!          inputs:                                                     !
!           (plyr,plvl,tlyr,tvly,clw,ciw,                              !
!            xlat,xlon,slmsk,                                          !
!            IX, NLAY, NLP1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!       'diagcld1'           --- diagnostic cloud calc routine         !
!          inputs:                                                     !
!           (plyr,plvl,tlyr,rhly,vvel,cv,cvt,cvb,                      !
!            xlat,xlon,slmsk,                                          !
!            IX, NLAY, NLP1,                                           !
!          outputs:                                                    !
!            clouds,clds,mtop,mbot)                                    !
!                                                                      !
!    internal accessable only subroutines:                             !
!       'gethml'             --- get diagnostic hi, mid, low clouds    !
!                                                                      !
!       'rhtable'            --- rh lookup table for diag cloud scheme !
!                                                                      !
!                                                                      !
!    cloud array description:                                          !
!                ---  for prognostic cloud: icldflg=1  ---             !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud liq water path                !
!          clouds(:,:,3)  -  mean effective radius for liquid cloud    !
!          clouds(:,:,4)  -  layer cloud ice water path                !
!          clouds(:,:,5)  -  mean effective radius for ice cloud       !
!          clouds(:,:,6)  -  layer rain drop water path                !
!          clouds(:,:,7)  -  mean effective radius for rain drop       !
!   **     clouds(:,:,8)  -  layer snow flake water path               !
!          clouds(:,:,9)  -  mean effective radius for snow flake      !
!   ** fu's scheme need to be normalized by snow density (g/m**3/1.0e6)!
!                ---  for diagnostic cloud: icldflg=0  ---             !
!          clouds(:,:,1)  -  layer total cloud fraction                !
!          clouds(:,:,2)  -  layer cloud optical depth                 !
!          clouds(:,:,3)  -  layer cloud single scattering albedo      !
!          clouds(:,:,4)  -  layer cloud asymmetry factor              !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module physparam'           in 'physparam.f'                  !
!       'module physcons'            in 'physcons.f'                   !
!       'module module_microphysics' in 'module_bfmicrophysics.f'      !
!                                                                      !
!                                                                      !
! program history log:                                                 !
!      nov 1992,   y.h., k.a.c, a.k. - cloud parameterization          !
!        'cldjms' patterned after slingo and slingo's work (jgr,       !
!        1992), stratiform clouds are allowed in any layer except      !
!        the surface and upper stratosphere. the relative humidity     !
!        criterion may cery in different model layers.                 !
!      mar 1993,   kenneth campana   - created original crhtab tables  !
!        for critical rh look up references.
!      feb 1994,   kenneth campana   - modified to use only one table  !
!        for all forecast hours.                                       !
!      oct 1995,   kenneth campana   - tuned cloud rh curves           !
!        rh-cld relation from tables created using mitchell-hahn       !
!        tuning technique on airforce rtneph observations.             !
!      nov 1995,   kenneth campana   - the bl relationships used       !
!        below llyr, except in marine stratus regions.                 !
!      apr 1996,   kenneth campana   - save bl cld amt in cld(,5)      !
!      aug 1997,   kenneth campana   - smooth out last bunch of bins   !
!        of the rh look up tables                                      !
!      dec 1998,   s. moorthi        - added prognostic cloud method   !
!      apr 2003,   yu-tai hou        - rewritten in frotran 90         !
!        modulized form 'module_rad_clouds' from combining the original!
!        subroutines 'cldjms', 'cldprp', and 'gcljms'. and seperated   !
!        prognostic and diagnostic methods into two packages.          !
!      --- 2003,   s. moorthi        - adapted b. ferrier's prognostic !
!        cloud scheme to ncep gfs model as additional option.          !
!      apr 2004,   yu-tai hou        - separated calculation of the    !
!        averaged h,m,l,bl cloud amounts from each of the cld schemes  !
!        to become an shared individule subprogram 'gethml'.           !
!      may 2004,   yu-tai hou        - rewritten ferrier's scheme as a !
!        separated program 'progcld2' in the cloud module.             !
!      apr 2005,   yu-tai hou        - modified cloud array and module !
!        structures.                                                   !
!      dec 2008,   yu-tai hou        - changed low-cld calculation,    !
!        now cantains clds from sfc layer and upward to the low/mid    !
!        boundary (include bl-cld). h,m,l clds domain boundaries are   !
!        adjusted for better agreement with observations.              !
!      jan 2011,   yu-tai hou        - changed virtual temperature     !
!        as input variable instead of originally computed inside the   !
!        two prognostic cld schemes 'progcld1' and 'progcld2'.         !
!      aug 2012,   yu-tai hou        - modified subroutine cld_init    !
!        to pass all fixed control variables at the start. and set     !
!        their correponding internal module variables to be used by    !
!        module subroutines.                                           !
!      nov 2012,   yu-tai hou        - modified control parameters     !
!        thru module 'physparam'.                                      !
!      apr 2013,   h.lin/y.hou       - corrected error in calculating  !
!        llyr for different vertical indexing directions.              !
!      jul 2013, r. sun/h. pan       - modified to use pdf cloud and   !
!        convective cloud cover and water for radiation                !
!                                                                      !
!      jul 2014 s. moorthi - merging with gfs version                  !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!

! \ingroup RRTMG
!> \defgroup module_radiation_clouds RRTMG Clouds Module
!! \brief This module computes cloud related quantities for radiation
!! computations.
!!
!! Knowledge of cloud properties and their vertical structure is
!! important for meteorological studies due to their impact on both the
!! Earth's radiation budget and adiabatic heating within the atmosphere.
!! Cloud properties in the US National Oceanic and Atmospheric
!! Administration National Centers for Environmental Prediction Global
!! Forecast System (GFS) include (i) cloud liquid/ice water path; (ii)
!! the fraction of clouds; (iii) effective radius of water/ice droplet:
!!\version NCEP-Radiation_clouds    v5.1  Nov 2012
!!
!! This module has six externally accessible subroutines:
!!  - cld_init():initialization routine
!!  - progcld1():zhao/moorthi prognostic cloud scheme
!!  - progcld2():Ferrier prognostic cloud microphysics
!!  - progcld3():Zhao/Moorthi prognostic cloud + pdfcld
!!  - progclduni():unified clouds with MG microphys
!!  - diagcld1():diagnostic cloud calculation routine
!!
!!  and two internally accessable only subroutines:
!!  - gethml():get diagnostic hi, mid, low,total,BL clouds
!!  - rhtable():rh lookup table for diag cloud scheme
!!
!> \section gen_al_clouds  RRTMG Clouds Module General Algorithm
!! @{
!! -# Cloud Liquid/Ice Water Path (LWP,IWP)
!!\n We define the fraction of liquid and ice cloud as:
!!\n Fraction of ice cloud (F): \f$F=(273.16K-T)/20\f$
!!\n LWP = total cloud condensate path X (1-F)
!!\n IWP = total clod condensate path X F
!!
!! -# GFS Cloud Fraction
!! \n The cloud fraction in a given grid box of the GFS model is
!! computed using the parameterization scheme of 
!! \cite xu_and_randall_1996 :
!!  \f[
!!  \sigma =RH^{k_{1}}\left[1-exp\left(-\frac{k_{2}q_{l}}{\left[\left(1-RH\right)q_{s}\right]^{k_{3}}}\right)\right]
!!  \f]
!!  Where \f$RH\f$ is relative humidity, \f$q_{l}\f$ is the cloud
!!  condensate, \f$q_{s}\f$ is saturation specific humidity,
!! \f$k_{1}(=0.25)\f$, \f$k_{2}(=100)\f$, \f$k_{3}(=0.49)\f$ are the
!! empirical parameters. The cloud condensate is partitioned into
!!  cloud water and ice in radiation based on temperature. Cloud drop
!! effective radius ranges 5-10 microns over land depending on
!!  temperature. Ice crystal radius is function of ice water content
!! (\cite heymsfield_and_mcfarquhar_1996).
!!  Maximum-randomly cloud overlapping is used in both long-wave
!! radiation and short-wave radiation. Convective clouds are not
!! considered in radiation.
!!\n
!! -# The parameterization of effective radius of water/ice droplet
!!    (\f$r_{e}\f$)
!>\n Two methods has been used to parameterize cloud properties in the
!! GFS model. The first method makes use of a diagnostic cloud scheme,
!! in which cloud properties are determined based on model-predicted
!! temperature, pressure, and boundary layer circulation from
!! \cite harshvardhan_et_al_1989 . The
!! diagnostic scheme is now replaced with a prognostic scheme that uses
!! cloud condensate information instead (NCEP Office Note 441).
!! \n For the parameterization of effective radius,\f$r_{ew}\f$, of
!! water droplet, we fix \f$r_{ew}\f$ to a value of \f$10\mu m\f$ over
!! the oceans. Over the land, \f$r_{ew}\f$ is defined as:
!!\f[
!! r_{ew} = 5+5 \times F
!!\f]
!! Thus, the effective radius of cloud water droplets will reach to a
!! minimum values of \f$5\mu m\f$ when F=0, and to a maximum value of
!! \f$10\mu m\f$ when the ice fraction is increasing.
!! \n For ice clouds, following \cite heymsfield_and_mcfarquhar_1996,
!! we have made the effective ice droplet radius to be an empirical
!! function of ice water concentration (IWC) and environmental temperature as:
!! \f[
!! r_{ei}=\begin{cases}(1250/9.917)IWC^{0.109} & T <-50^0C \\(1250/9.337)IWC^{0.080} & -50^0C \leq T<-40^0C\\(1250/9.208)IWC^{0.055} & -40^0C\leq T <-30^0C\\(1250/9.387)IWC^{0.031} & -30^0C \leq T\end{cases}
!! \f]
!! where IWC and IWP satisfy:
!! \f[
!! IWP_{\triangle Z}=\int_{\triangle Z} IWCdZ
!! \f]
!! @}
!========================================!
      module module_radiation_clouds     !
!........................................!
!
      use physparam,           only : icldflg, icmphys, iovrsw, iovrlw, &
     &                                lcrick, lcnorm, lnoprec,          &
     &                                ivflip, kind_phys, kind_io4
      use physcons,            only : con_fvirt, con_ttp, con_rocp,     &
     &                                con_t0c, con_pi, con_g, con_rd,   &
     &                                con_thgni
      use module_microphysics, only : rsipath2
      use module_iounitdef,    only : NICLTUN
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGCLD='NCEP-Radiation_clouds    v5.1  Nov 2012 '
!    &   VTAGCLD='NCEP-Radiation_clouds    v5.0  Aug 2012 '

!  ---  set constant parameters
      real (kind=kind_phys), parameter :: gfac=1.0e5/con_g              &
     &,                                   gord=con_g/con_rd
! number of fields in cloud array
      integer, parameter, public :: NF_CLDS = 9
! number of cloud vertical domains
      integer, parameter, public :: NK_CLDS = 3

! pressure limits of cloud domain interfaces (low,mid,high) in mb (0.1kPa)
      real (kind=kind_phys), save :: ptopc(NK_CLDS+1,2)

!org  data ptopc / 1050., 642., 350., 0.0,  1050., 750., 500., 0.0 /
      data ptopc / 1050., 650., 400., 0.0,  1050., 750., 500., 0.0 /

!     real (kind=kind_phys), parameter :: climit = 0.01
      real (kind=kind_phys), parameter :: climit = 0.001, climit2=0.05
      real (kind=kind_phys), parameter :: ovcst  = 1.0 - 1.0e-8

! default liq radius to 10 micron
      real (kind=kind_phys), parameter :: reliq_def = 10.0
! default ice radius to 50 micron
      real (kind=kind_phys), parameter :: reice_def = 50.0
! default rain radius to 1000 micron
      real (kind=kind_phys), parameter :: rrain_def = 1000.0
! default snow radius to 250 micron
      real (kind=kind_phys), parameter :: rsnow_def = 250.0

! rh in one percent interval
      integer, parameter :: NBIN=100
! =1,2 for eastern and western hemispheres
      integer, parameter :: NLON=2
! =1,4 for 60n-30n,30n-equ,equ-30s,30s-60s
      integer, parameter :: NLAT=4
! =1,4 for bl,low,mid,hi cld type
      integer, parameter :: MCLD=4
! =1,2 for land,sea
      integer, parameter :: NSEAL=2

! default cld single scat albedo
      real (kind=kind_phys), parameter :: cldssa_def = 0.99
! default cld asymmetry factor
      real (kind=kind_phys), parameter :: cldasy_def = 0.84

!  ---  xlabdy: lat bndry between tuning regions, +/- xlim for transition
!       xlobdy: lon bndry between tuning regions
! lat bndry between tuning regions
      real (kind=kind_phys)            :: xlabdy(3)
! lon bndry between tuning regions
      real (kind=kind_phys)            :: xlobdy(3)
! +/- xlim for transition
      real (kind=kind_phys), parameter :: xlim=5.0

      data xlabdy / 30.0,  0.0, -30.0 /,  xlobdy / 0.0, 180., 360. /

! low cloud vertical velocity adjustment boundaries in mb/sec
      real (kind=kind_phys), parameter :: vvcld1= 0.0003e0
! low cloud vertical velocity adjustment boundaries in mb/sec
      real (kind=kind_phys), parameter :: vvcld2=-0.0005e0

!  ---  those data will be set up by "cld_init"
!     rhcl : tuned rh relation table for diagnostic cloud scheme

! tuned relative humidity relation table for diagnostic cloud scheme
      real (kind=kind_phys) :: rhcl(NBIN,NLON,NLAT,MCLD,NSEAL)

! upper limit of boundary layer clouds
      integer  :: llyr   = 2
! maximum-random cloud overlapping method
      integer  :: iovr   = 1

      public progcld1, progcld2, progcld3, progclduni, diagcld1,        &
     &       cld_init


! =================
      contains
! =================

!> \ingroup module_radiation_clouds
!> This subroutine is an initialization program for cloud-radiation
!! calculations and sets up boundary layer cloud top.
!!\param si              model vertical sigma layer interface
!!\param NLAY            vertical layer number
!!\param me              print control flag
!>\section gen_cld_init cld_init General Algorithm
!> @{
      subroutine cld_init                                               &
     &     ( si, NLAY, me ) !  ---  inputs
!  ---  outputs:
!          ( none )

!  ===================================================================  !
!                                                                       !
! abstract: cld_init is an initialization program for cloud-radiation   !
!   calculations. it sets up boundary layer cloud top.                  !
!                                                                       !
!                                                                       !
! inputs:                                                               !
!   si (L+1)        : model vertical sigma layer interface              !
!   NLAY            : vertical layer number                             !
!   me              : print control flag                                !
!                                                                       !
!  outputs: (none)                                                      !
!           to module variables                                         !
!                                                                       !
!  external module variables: (in physparam)                            !
!   icldflg         : cloud optical property scheme control flag        !
!                     =0: model use diagnostic cloud method             !
!                     =1: model use prognostic cloud method             !
!   icmphys         : cloud microphysics scheme control flag            !
!                     =1: zhao/carr/sundqvist microphysics cloud        !
!                     =2: ferrier microphysics cloud scheme             !
!                     =3: zhao/carr/sundqvist microphysics cloud +pdfcld!
!   iovrsw/iovrlw   : sw/lw control flag for cloud overlapping scheme   !
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!   ivflip          : control flag for direction of vertical index      !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!  usage:       call cld_init                                           !
!                                                                       !
!  subroutines called:    rhtable                                       !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, me

      real (kind=kind_phys), intent(in) :: si(:)

!  ---  outputs: (none)

!  ---  locals:
      integer :: k, kl, ier

!
!===> ...  begin here
!
!  ---  set up module variables

      iovr    = max( iovrsw, iovrlw )    !cld ovlp used for diag HML cld output

      if (me == 0) print *, VTAGCLD      !print out version tag

      if ( icldflg == 0 ) then
        if (me == 0) print *,' - Using Diagnostic Cloud Method'

!> -# Call rhtable() to set up tuned relative humidity table.

        call rhtable( me, ier )

        if (ier < 0) then
          write(6,99) ier
  99      format(3x,' *** Error in finding tuned RH table ***'          &
     &,         /3x,'     STOP at calling subroutine RHTABLE !!'/)
          stop 99
        endif
      else
        if (me == 0) then
          print *,' - Using Prognostic Cloud Method'
          if (icmphys == 1) then
            print *,'   --- Zhao/Carr/Sundqvist microphysics'
          elseif (icmphys == 2) then
            print *,'   --- Ferrier cloud microphysics'
          elseif (icmphys == 3) then
            print *,'   --- zhao/carr/sundqvist + pdf cloud'
          else
            print *,'  !!! ERROR in cloud microphysc specification!!!', &
     &              '  icmphys (NP3D) =',icmphys
            stop
          endif
        endif
      endif

!> -# Compute the top of BL cld (llyr), which is the topmost non
!!    cld(low) layer for stratiform (at or above lowest 0.1 of the
!!     atmosphere).

      if ( ivflip == 0 ) then    ! data from toa to sfc
        lab_do_k0 : do k = NLAY, 2, -1
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k0
        enddo  lab_do_k0

        llyr = kl
      else                      ! data from sfc to top
        lab_do_k1 : do k = 2, NLAY
          kl = k
          if (si(k) < 0.9e0) exit lab_do_k1
        enddo  lab_do_k1

        llyr = kl - 1
      endif                     ! end_if_ivflip

!
      return
!...................................
      end subroutine cld_init
!-----------------------------------
!> @}

!> \ingroup module_radiation_clouds
!> This subroutine computes cloud related quantities using
!! zhao/moorthi's prognostic cloud microphysics scheme.
!!\param plyr    (IX,NLAY), model layer mean pressure in mb (100Pa)
!!\param plvl    (IX,NLP1), model level pressure in mb (100Pa)
!!\param tlyr    (IX,NLAY), model layer mean temperature in K
!!\param tvly    (IX,NLAY), model layer virtual temperature in K
!!\param qlyr    (IX,NLAY), layer specific humidity in gm/gm
!!\param qstl    (IX,NLAY), layer saturate humidity in gm/gm
!!\param rhly    (IX,NLAY), layer relative humidity \f$ (=qlyr/qstl) \f$
!!\param clw     (IX,NLAY), layer cloud condensate amount
!!\param xlat    (IX), grid latitude in radians, default to pi/2 ->
!!               -pi/2 range, otherwise see in-line comment
!!\param xlon    (IX), grid longitude in radians  (not used)
!!\param slmsk   (IX), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param IX           horizontal dimention
!!\param NLAY,NLP1    vertical layer/level dimensions
!!\param clouds      (IX,NLAY,NF_CLDS), cloud profiles
!!\n                 (:,:,1) - layer total cloud fraction
!!\n                 (:,:,2) - layer cloud liq water path \f$(g/m^2)\f$
!!\n                 (:,:,3) - mean eff radius for liq cloud (micron)
!!\n                 (:,:,4) - layer cloud ice water path \f$(g/m^2)\f$
!!\n                 (:,:,5) - mean eff radius for ice cloud (micron)
!!\n                 (:,:,6) - layer rain drop water path (not assigned)
!!\n                 (:,:,7) - mean eff radius for rain drop (micron)
!!\n                 (:,:,8) - layer snow flake water path (not assigned)
!!\n                 (:,:,9) - mean eff radius for snow flake (micron)
!!\n  *** fu's scheme need to be normalized by snow density \f$ (g/m^3/1.0e6)\f$
!!\param clds       (IX,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop       (IX,3), vertical indices for low, mid, hi cloud tops
!!\param mbot       (IX,3), vertical indices for low, mid, hi cloud bases
!>\section gen_progcld1 progcld1 General Algorithm
!> @{
      subroutine progcld1                                               &
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &    !  ---  inputs:
     &       xlat,xlon,slmsk, IX, NLAY, NLP1,                           &
     &       uni_cld, lmfshal, lmfdeep2, cldcov,                        &
     &       clouds,clds,mtop,mbot                                      &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld1    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!   uni_cld         : logical - true for cloud fraction from shoc       !
!   lmfshal         : logical - true for mass flux shallow convection   !
!   lmfdeep2        : logical - true for mass flux deep convection      !
!   cldcov          : layer cloud fraction (used when uni_cld=.true.    !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      logical, intent(in)  :: uni_cld, lmfshal, lmfdeep2

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw, cldcov

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!  ---  constant values
!     real (kind=kind_phys), parameter :: xrc3 = 200.
      real (kind=kind_phys), parameter :: xrc3 = 100.

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, IX
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, NLAY-1
          do i = 1, IX
            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!> -# Find top pressure for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!> -# Compute cloud liquid/ice condensate path in \f$ g/m^2 \f$ .

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            clwt     = max(0.0, clwf(i,k)) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            clwt     = max(0.0, clwf(i,k)) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      endif                            ! end_if_ivflip

!> -# Compute effective liquid cloud droplet radius over land.

      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

      if (uni_cld) then     ! use unified sgs clouds generated outside
        do k = 1, NLAY
          do i = 1, IX
            cldtot(i,k) = cldcov(i,k)
          enddo
        enddo

      else

!> -# Calculate layer cloud fraction.

      if ( ivflip == 0 ) then              ! input data from toa to sfc

        clwmin = 0.0
        if (.not. lmfshal) then
          do k = NLAY, 1, -1
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1
!             tem1  = 1000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = NLAY, 1, -1
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )
              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      else                                 ! input data from sfc to toa

        clwmin = 0.0
        if (.not. lmfshal) then
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then

              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1

!             tem1  = 1000.0 / tem1

              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = 1, NLAY
          do i = 1, IX
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
!           clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clwf(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)  !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clwf(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      endif                                ! end_if_flip
      endif                                ! if (uni_cld) then

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Compute effective ice cloud droplet radius following 
!!    \cite heymsfield_and_mcfarquhar_1996.

      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo


!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!!    cloud fractions and clouds top/bottom layer indices for low, mid,
!!    and high clouds.
!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld1
!-----------------------------------
!> @}

!> \ingroup module_radiation_clouds
!> This subroutine computes cloud related quantities using ferrier's
!! prognostic cloud microphysics scheme.
!!\param plyr    (IX,NLAY), model layer mean pressure in mb (100Pa)
!!\param plvl    (IX,NLP1), model level pressure in mb (100Pa)
!!\param tlyr    (IX,NLAY), model layer mean temperature in K
!!\param tvly    (IX,NLAY), model layer virtual temperature in K
!!\param qlyr    (IX,NLAY), layer specific humidity in gm/gm
!!\param qstl    (IX,NLAY), layer saturate humidity in gm/gm
!!\param rhly    (IX,NLAY), layer relative humidity (=qlyr/qstl)
!!\param clw     (IX,NLAY), layer cloud condensate amount
!!\param f_ice   (IX,NLAY), fraction of layer cloud ice  (ferrier micro-phys)
!!\param f_rain  (IX,NLAY), fraction of layer rain water (ferrier micro-phys)
!!\param r_rime  (IX,NLAY), mass ratio of total ice to unrimed ice (>=1)
!!\param flgmin  (IX), minimum large ice fraction
!!\param xlat    (IX), grid latitude in radians, default to pi/2 ->
!!               -pi/2 range, otherwise see in-line comment
!!\param xlon    (IX), grid longitude in radians  (not used)
!!\param slmsk   (IX), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param IX      horizontal dimention
!!\param NLAY,NLP1    vertical layer/level dimensions
!!\param clouds      (IX,NLAY,NF_CLDS), cloud profiles
!!\n                 (:,:,1) - layer total cloud fraction
!!\n                 (:,:,2) - layer cloud liq water path  \f$(g/m^2)\f$
!!\n                 (:,:,3) - mean eff radius for liq cloud (micron)
!!\n                 (:,:,4) - layer cloud ice water path  \f$(g/m^2)\f$
!!\n                 (:,:,5) - mean eff radius for ice cloud (micron)
!!\n                 (:,:,6) - layer rain drop water path  \f$(g/m^2)\f$
!!\n                 (:,:,7) - mean eff radius for rain drop (micron)
!!\n                 (:,:,8) - layer snow flake water path \f$(g/m^2)\f$
!!\n                 (:,:,9) - mean eff radius for snow flake (micron)
!!\n  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6)
!!\param clds  (IX,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop  (IX,3), vertical indices for low, mid, hi cloud tops
!!\param mbot  (IX,3), vertical indices for low, mid, hi cloud bases
!>\section gen_progcld2 progcld2 General Algorithm
!> @{
      subroutine progcld2                                               &
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,                    &    !  ---  inputs:
     &       xlat,xlon,slmsk, f_ice,f_rain,r_rime,flgmin,               &
     &       IX, NLAY, NLP1, lmfshal, lmfdeep2,                         &
     &       clouds,clds,mtop,mbot                                      &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld2    computes cloud related quantities using    !
!   ferrier's prognostic cloud microphysics scheme.                     !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld2                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   qlyr  (IX,NLAY) : layer specific humidity in gm/gm                  !
!   qstl  (IX,NLAY) : layer saturate humidity in gm/gm                  !
!   rhly  (IX,NLAY) : layer relative humidity (=qlyr/qstl)              !
!   clw   (IX,NLAY) : layer cloud condensate amount                     !
!   f_ice (IX,NLAY) : fraction of layer cloud ice  (ferrier micro-phys) !
!   f_rain(IX,NLAY) : fraction of layer rain water (ferrier micro-phys) !
!   r_rime(IX,NLAY) : mass ratio of total ice to unrimed ice (>=1)      !
!   flgmin(IX)      : minimim large ice fraction                        !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         (g/m**2)      !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        (g/m**2)      !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! external module variables:                                            !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!   lnoprec         : precip effect in radiation flag (ferrier scheme)  !
!                     =t: snow/rain has no impact on radiation          !
!                     =f: snow/rain has impact on radiation             !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  constants

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      logical, intent(in)  :: lmfshal, lmfdeep2

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, qlyr, qstl, rhly, clw, f_ice, f_rain, r_rime

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk
      real (kind=kind_phys), dimension(:), intent(in) :: flgmin

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, tem2d, clw2,       &
     &       qcwat, qcice, qrain, fcice, frain, rrime, rsden, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id

!  ---  constant values
!     real (kind=kind_phys), parameter :: xrc3 = 200.
      real (kind=kind_phys), parameter :: xrc3 = 100.

!
!===> ... begin here
!
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          fcice (i,k) = max(0.0, min(1.0, f_ice(i,k)))
          frain (i,k) = max(0.0, min(1.0, f_rain(i,k)))
          rrime (i,k) = max(1.0, r_rime(i,k))
          tem2d (i,k) = tlyr(i,k) - con_t0c
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, IX
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, NLAY-1
          do i = 1, IX
            clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, NLAY
          do i = 1, IX
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

!> -# Find top pressure (ptopc) for each cloud domain for given latitude.
!    - ptopc(k,i): top pressure of each cld domain (k=1-4 are sfc,l,m,
!!     h; i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!> -# Seperate cloud condensate into liquid, ice, and rain types, and
!! save the liquid+ice condensate in array clw2 for later calculation
!!  of cloud fraction.

      do k = 1, NLAY
        do i = 1, IX
          if (tem2d(i,k) > -40.0) then
            qcice(i,k) = clwf(i,k) * fcice(i,k)
            tem1       = clwf(i,k) - qcice(i,k)
            qrain(i,k) = tem1 * frain(i,k)
            qcwat(i,k) = tem1 - qrain(i,k)
            clw2 (i,k) = qcwat(i,k) + qcice(i,k)
          else
            qcice(i,k) = clwf(i,k)
            qrain(i,k) = 0.0
            qcwat(i,k) = 0.0
            clw2 (i,k) = clwf(i,k)
          endif
        enddo
      enddo

!> -# Call module_microphysics::rsipath2(), in Ferrier's scheme, to
!! compute layer's cloud liquid, ice, rain, and snow water condensate
!! path and the partical effective radius for liquid droplet, rain drop,
!! and snow flake.
      call  rsipath2                                                    &
!  ---  inputs:
     &     ( plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime,        &
     &       IX, NLAY, ivflip, flgmin,                                  &
!  ---  outputs:
     &       cwp, cip, crp, csp, rew, rer, res, rsden                   &
     &     )


      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, NLAY
          do i = 1, IX
            tem2d(i,k) = (con_g * plyr(i,k))                            &
     &                 / (con_rd* (plvl(i,k+1) - plvl(i,k)))
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, NLAY
          do i = 1, IX
            tem2d(i,k) = (con_g * plyr(i,k))                            &
     &                 / (con_rd* (plvl(i,k) - plvl(i,k+1)))
          enddo
        enddo
      endif                            ! end_if_ivflip

!> -# Calculate layer cloud fraction.

      if ( ivflip == 0 ) then              ! input data from toa to sfc

        clwmin = 0.0
        if (.not. lmfshal) then
          do k = NLAY, 1, -1
          do i = 1, IX
!           clwt = 1.0e-7 * (plyr(i,k)*0.001)
!           clwt = 1.0e-6 * (plyr(i,k)*0.001)
            clwt = 2.0e-6 * (plyr(i,k)*0.001)
!           clwt = 5.0e-6 * (plyr(i,k)*0.001)
!           clwt = 5.0e-6

            if (clw2(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

!             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
!             tem1  = 100.0 / tem1

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1
!             tem1  = 2400.0 / tem1
!cnt          tem1  = 2500.0 / tem1
!             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
!             tem1  = 2000.0 / tem1
!             tem1  = 1000.0 / tem1
!             tem1  = 100.0 / tem1

              value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = NLAY, 1, -1
          do i = 1, IX
!           clwt = 1.0e-6 * (plyr(i,k)*0.001)
            clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clw2(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)    !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      else                                 ! input data from sfc to toa

        clwmin = 0.0e-6
        if (.not. lmfshal) then
          do k = 1, NLAY
          do i = 1, IX
!           clwt = 1.0e-7 * (plyr(i,k)*0.001)
!           clwt = 1.0e-6 * (plyr(i,k)*0.001)
            clwt = 2.0e-6 * (plyr(i,k)*0.001)
!           clwt = 5.0e-6 * (plyr(i,k)*0.001)
!           clwt = 5.0e-6

            if (clw2(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )

!             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
!             tem1  = 100.0 / tem1

              tem1  = min(max(sqrt(sqrt(onemrh*qstl(i,k))),0.0001),1.0)
              tem1  = 2000.0 / tem1
!             tem1  = 2400.0 / tem1
!cnt          tem1  = 2500.0 / tem1
!             tem1  = min(max(sqrt(onemrh*qstl(i,k)),0.0001),1.0)
!             tem1  = 2000.0 / tem1
!             tem1  = 1000.0 / tem1
!             tem1  = 100.0 / tem1

              value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        else
          do k = 1, NLAY
          do i = 1, IX
!           clwt = 1.0e-6 * (plyr(i,k)*0.001)
            clwt = 2.0e-6 * (plyr(i,k)*0.001)

            if (clw2(i,k) > clwt) then
              onemrh= max( 1.e-10, 1.0-rhly(i,k) )
              clwm  = clwmin / max( 0.01, plyr(i,k)*0.001 )
!
              tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)   !jhan
              if (lmfdeep2) then
                tem1  = xrc3 / tem1
              else
                tem1  = 100.0 / tem1
              endif
!
              value = max( min( tem1*(clw2(i,k)-clwm), 50.0 ), 0.0 )
              tem2  = sqrt( sqrt(rhly(i,k)) )

              cldtot(i,k) = max( tem2*(1.0-exp(-value)), 0.0 )
            endif
          enddo
          enddo
        endif

      endif                                ! end_if_flip

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

!     When lnoprec = .true. snow/rain has no impact on radiation
      if ( lnoprec ) then
        do k = 1, NLAY
          do i = 1, IX
            crp(i,k) = 0.0
            csp(i,k) = 0.0
          enddo
        enddo
      endif
!
      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Calculate effective ice cloud droplet radius.

      do k = 1, NLAY
        do i = 1, IX
          tem1 = tlyr(i,k) - con_ttp
          tem2 = cip(i,k)

          if (tem2 > 0.0) then
            tem3 = tem2d(i,k) * tem2 / tvly(i,k)

            if (tem1 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem1 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem1 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif

!           if (lprnt .and. k == l) print *,' reiL=',rei(i,k),' icec=', &
!    &        icec,' cip=',cip(i,k),' tem=',tem,' delt=',delt

            rei(i,k)   = max(10.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!!!!        rei(i,k)   = max(30.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(50.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(100.0, min(rei(i,k), 300.0))
          endif
        enddo
      enddo
!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
          clouds(i,k,6) = crp(i,k)
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = csp(i,k)               !ncar scheme
          clouds(i,k,8) = csp(i,k) * rsden(i,k)  !fu's scheme
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo


!> -# Call gethml(), to compute low, mid, high, total, and boundary
!! layer cloud fractions and clouds top/bottom layer indices for low,
!! mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld2
!-----------------------------------
!> @}

!> \ingroup module_radiation_clouds
!> This subroutine computes cloud related quantities using
!! zhao/moorthi's prognostic cloud microphysics scheme + pdfcld.
!!\param plyr       (ix,nlay), model layer mean pressure in mb (100pa)
!!\param plvl       (ix,nlp1), model level pressure in mb (100pa)
!!\param tlyr       (ix,nlay), model layer mean temperature in K
!!\param tvly       (ix,nlay), model layer virtual temperature in K
!!\param qlyr       (ix,nlay), layer specific humidity in gm/gm
!!\param qstl       (ix,nlay), layer saturate humidity in gm/gm
!!\param rhly       (ix,nlay), layer relative humidity (=qlyr/qstl)
!!\param clw        (ix,nlay), layer cloud condensate amount
!!\param xlat       (ix), grid latitude in radians, default to pi/2 ->
!!                  -pi/2 range, otherwise see in-line comment
!!\param xlon       (ix), grid longitude in radians  (not used)
!!\param slmsk      (ix), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param ix         horizontal dimention
!!\param nlay,nlp1  vertical layer/level dimensions
!!\param deltaq     (ix,nlay), half total water distribution width
!!\param sup        supersaturation
!!\param kdt
!!\param me         print control flag
!!\param clouds    (ix,nlay,nf_clds), cloud profiles
!!\n               (:,:,1) - layer total cloud fraction
!!\n               (:,:,2) - layer cloud liq water path (g/m**2)
!!\n               (:,:,3) - mean eff radius for liq cloud (micron)
!!\n               (:,:,4) - layer cloud ice water path (g/m**2)
!!\n               (:,:,5) - mean eff radius for ice cloud (micron)
!!\n               (:,:,6) - layer rain drop water path         not assigned
!!\n               (:,:,7) - mean eff radius for rain drop (micron)
!!\n               (:,:,8) - layer snow flake water path        not assigned
!!\n               (:,:,9) - mean eff radius for snow flake(micron)
!!\param clds      (ix,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop      (ix,3), vertical indices for low, mid, hi cloud tops
!!\param mbot      (ix,3), vertical indices for low, mid, hi cloud bases
!>\section gen_progcld3 progcld3 General Algorithm
!> @{
      subroutine progcld3                                               &
     &     ( plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,clw,cnvw,cnvc,          &    !  ---  inputs:
     &       xlat,xlon,slmsk,                                           &
     &       ix, nlay, nlp1,                                            &
     &       deltaq,sup,kdt,me,                                         &
     &       clouds,clds,mtop,mbot                                      &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progcld3    computes cloud related quantities using    !
!   zhao/moorthi's prognostic cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progcld3                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (ix,nlay) : model layer mean pressure in mb (100pa)           !
!   plvl  (ix,nlp1) : model level pressure in mb (100pa)                !
!   tlyr  (ix,nlay) : model layer mean temperature in k                 !
!   tvly  (ix,nlay) : model layer virtual temperature in k              !
!   qlyr  (ix,nlay) : layer specific humidity in gm/gm                  !
!   qstl  (ix,nlay) : layer saturate humidity in gm/gm                  !
!   rhly  (ix,nlay) : layer relative humidity (=qlyr/qstl)              !
!   clw   (ix,nlay) : layer cloud condensate amount                     !
!   xlat  (ix)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (ix)      : grid longitude in radians  (not used)             !
!   slmsk (ix)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   ix              : horizontal dimention                              !
!   nlay,nlp1       : vertical layer/level dimensions                   !
!   cnvw  (ix,nlay) : layer convective cloud condensate                 !
!   cnvc  (ix,nlay) : layer convective cloud cover                      !
!   deltaq(ix,nlay) : half total water distribution width               !
!   sup             : supersaturation                                   !

!                                                                       !
! output variables:                                                     !
!   clouds(ix,nlay,nf_clds) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (ix,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (ix,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (ix,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lcrick          : control flag for eliminating crick                !
!                     =t: apply layer smoothing to eliminate crick      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: ix, nlay, nlp1,kdt

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,    &
     &       tlyr, tvly, qlyr, qstl, rhly, clw
!     &       tlyr, tvly, qlyr, qstl, rhly, clw, cnvw, cnvc
!      real (kind=kind_phys), dimension(:,:), intent(in) :: deltaq
      real (kind=kind_phys), dimension(:,:) :: deltaq, cnvw, cnvc
      real (kind=kind_phys) qtmp,qsc,rhs
      real (kind=kind_phys), intent(in) :: sup
      real (kind=kind_phys), parameter :: epsq = 1.0e-12

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,    &
     &       slmsk
      integer :: me

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(ix,nlay) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real (kind=kind_phys) :: ptop1(ix,nk_clds+1)

      real (kind=kind_phys) :: clwmin, clwm, clwt, onemrh, value,       &
     &       tem1, tem2, tem3

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, nlay
        do i = 1, ix
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
          clwf(i,k)   = 0.0
        enddo
      enddo
!
      if ( lcrick ) then
        do i = 1, ix
          clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
          clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
        enddo
        do k = 2, nlay-1
          do i = 1, ix
            clwf(i,k) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
          enddo
        enddo
      else
        do k = 1, nlay
          do i = 1, ix
            clwf(i,k) = clw(i,k)
          enddo
        enddo
      endif

      if(kdt==1) then
        do k = 1, nlay
          do i = 1, ix
            deltaq(i,k) = (1.-0.95)*qstl(i,k)
          enddo
        enddo
      endif

!> -# Find top pressure (ptopc) for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,l,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, ix
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!> -# Calculate liquid/ice condensate path in \f$ g/m^2 \f$

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, nlay
          do i = 1, ix
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            clwt  = max(0.0,(clwf(i,k)+cnvw(i,k))) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, nlay
          do i = 1, ix
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            clwt = max(0.0,(clwf(i,k)+cnvw(i,k))) * gfac * delp(i,k)
            cip(i,k) = clwt * tem2d(i,k)
            cwp(i,k) = clwt - cip(i,k)
          enddo
        enddo
      endif                            ! end_if_ivflip

!> -# Calculate effective liquid cloud droplet radius over land.

      do i = 1, ix
        if (nint(slmsk(i)) == 1) then
          do k = 1, nlay
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo

!> -# Calculate layer cloud fraction.

      if ( ivflip == 0 ) then              ! input data from toa to sfc
          do k = nlay, 1, -1
          do i = 1, ix
            tem1 = tlyr(i,k) - 273.16
            if(tem1 < con_thgni) then  ! for pure ice, has to be consistent with gscond
              qsc = sup * qstl(i,k)
              rhs = sup
            else
              qsc = qstl(i,k)
              rhs = 1.0
            endif
           if(rhly(i,k) >= rhs) then
             cldtot(i,k) = 1.0
           else
             qtmp = qlyr(i,k) + clwf(i,k) - qsc
             if(deltaq(i,k) > epsq) then
               if(qtmp < -deltaq(i,k) .or. clwf(i,k) < epsq) then
!               if(qtmp < -deltaq(i,k)) then
                 cldtot(i,k) = 0.0
               elseif(qtmp >= deltaq(i,k)) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.5*qtmp/deltaq(i,k) + 0.5
                 cldtot(i,k) = max(cldtot(i,k),0.0)
                 cldtot(i,k) = min(cldtot(i,k),1.0)
               endif
             else
               if(qtmp.gt.0) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.0
               endif
             endif
           endif
           cldtot(i,k) = cnvc(i,k)+(1-cnvc(i,k))*cldtot(i,k)
           cldtot(i,k) = max(cldtot(i,k),0.)
           cldtot(i,k) = min(cldtot(i,k),1.)
          enddo
          enddo
      else                                 ! input data from sfc to toa
          do k = 1, nlay
          do i = 1, ix
            tem1 = tlyr(i,k) - 273.16
            if(tem1 < con_thgni) then  ! for pure ice, has to be consistent with gscond
              qsc = sup * qstl(i,k)
              rhs = sup
            else
              qsc = qstl(i,k)
              rhs = 1.0
            endif
           if(rhly(i,k) >= rhs) then
             cldtot(i,k) = 1.0
           else
             qtmp = qlyr(i,k) + clwf(i,k) - qsc
             if(deltaq(i,k) > epsq) then
!              if(qtmp <= -deltaq(i,k) .or. cwmik < epsq) then
               if(qtmp <= -deltaq(i,k)) then
                 cldtot(i,k) = 0.0
               elseif(qtmp >= deltaq(i,k)) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.5*qtmp/deltaq(i,k) + 0.5
                 cldtot(i,k) = max(cldtot(i,k),0.0)
                 cldtot(i,k) = min(cldtot(i,k),1.0)
               endif
             else
               if(qtmp > 0.) then
                 cldtot(i,k) = 1.0
               else
                 cldtot(i,k) = 0.0
               endif
             endif
           endif
           cldtot(i,k) = cnvc(i,k) + (1-cnvc(i,k))*cldtot(i,k)
           cldtot(i,k) = max(cldtot(i,k),0.)
           cldtot(i,k) = min(cldtot(i,k),1.)

          enddo
          enddo
      endif                                ! end_if_flip

      do k = 1, nlay
        do i = 1, ix
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, nlay
          do i = 1, ix
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Calculate effective ice cloud droplet radius following 
!! \cite heymsfield_and_mcfarquhar_1996.

      do k = 1, nlay
        do i = 1, ix
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
!           tem3 = gord * cip(i,k) * (plyr(i,k)/delp(i,k)) / tvly(i,k)
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0, min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, nlay
        do i = 1, ix
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo


!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!! cloud fractions and clouds top/bottom layer indices for low, mid,
!! and high clouds.
!       the three cloud domain boundaries are defined by ptopc.  the cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.


      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       ix,nlay,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progcld3
!-----------------------------------
!> @}

!> \ingroup module_radiation_clouds
!> This subroutine computes cloud related quantities using
!! zhao/moorthi's prognostic cloud microphysics scheme.
!!\param plyr    (IX,NLAY), model layer mean pressure in mb (100Pa)
!!\param plvl    (IX,NLP1), model level pressure in mb (100Pa)
!!\param tlyr    (IX,NLAY), model layer mean temperature in K
!!\param tvly    (IX,NLAY), model layer virtual temperature in K
!!\param qlyr    (IX,NLAY), layer specific humidity in gm/gm
!!\param clw     (IX,NLAY), layer cloud liquid water amount
!!\param ciw     (IX,NLAY), layer cloud ice water amount
!!\param xlat    (IX), grid latitude in radians, default to pi/2 ->
!!               -pi/2 range, otherwise see in-line comment
!!\param xlon    (IX), grid longitude in radians  (not used)
!!\param slmsk   (IX), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param IX           horizontal dimention
!!\param NLAY,NLP1    vertical layer/level dimensions
!!\param clouds      (IX,NLAY,NF_CLDS), cloud profiles
!!\n                 (:,:,1) - layer total cloud fraction
!!\n                 (:,:,2) - layer cloud liq water path \f$(g/m^2)\f$
!!\n                 (:,:,3) - mean eff radius for liq cloud (micron)
!!\n                 (:,:,4) - layer cloud ice water path \f$(g/m^2)\f$
!!\n                 (:,:,5) - mean eff radius for ice cloud (micron)
!!\n                 (:,:,6) - layer rain drop water path (not assigned)
!!\n                 (:,:,7) - mean eff radius for rain drop (micron)
!!\n                 (:,:,8) - layer snow flake water path (not assigned)
!!\n                 (:,:,9) - mean eff radius for snow flake (micron)
!!\n  *** fu's scheme need to be normalized by snow density \f$ (g/m^3/1.0e6)\f$
!!\param clds       (IX,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop       (IX,3), vertical indices for low, mid, hi cloud tops
!!\param mbot       (IX,3), vertical indices for low, mid, hi cloud bases
!>\section gen_progclduni progclduni General Algorithm
!> @{
      subroutine progclduni                                             &
     &     ( plyr,plvl,tlyr,tvly,clw,ciw,                               &    !  ---  inputs:
     &       xlat,xlon,slmsk, IX, NLAY, NLP1, cldcov,                   &
     &       clouds,clds,mtop,mbot                                      &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    progclduni    computes cloud related quantities using    !
!   for unified cloud microphysics scheme.                !
!                                                                       !
! abstract:  this program computes cloud fractions from cloud           !
!   condensates, calculates liquid/ice cloud droplet effective radius,  !
!   and computes the low, mid, high, total and boundary layer cloud     !
!   fractions and the vertical indices of low, mid, and high cloud      !
!   top and base.  the three vertical cloud domains are set up in the   !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call progclduni                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   tvly  (IX,NLAY) : model layer virtual temperature in k              !
!   clw   (IX,NLAY) : layer cloud liquid water amount                   !
!   ciw   (IX,NLAY) : layer cloud ice water amount                      !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians  (not used)             !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud liq water path         (g/m**2)      !
!      clouds(:,:,3) - mean eff radius for liq cloud      (micron)      !
!      clouds(:,:,4) - layer cloud ice water path         (g/m**2)      !
!      clouds(:,:,5) - mean eff radius for ice cloud      (micron)      !
!      clouds(:,:,6) - layer rain drop water path         not assigned  !
!      clouds(:,:,7) - mean eff radius for rain drop      (micron)      !
!  *** clouds(:,:,8) - layer snow flake water path        not assigned  !
!      clouds(:,:,9) - mean eff radius for snow flake     (micron)      !
!  *** fu's scheme need to be normalized by snow density (g/m**3/1.0e6) !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! module variables:                                                     !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   lmfshal         : mass-flux shallow conv scheme flag                !
!   lmfdeep2        : scale-aware mass-flux deep conv scheme flag       !
!   lcrick          : control flag for eliminating CRICK                !
!                     =t: apply layer smoothing to eliminate CRICK      !
!                     =f: do not apply layer smoothing                  !
!   lcnorm          : control flag for in-cld condensate                !
!                     =t: normalize cloud condensate                    !
!                     =f: not normalize cloud condensate                !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, tvly, clw, ciw, cldcov

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d
!    &       cwp, cip, crp, csp, rew, rei, res, rer, delp, tem2d, clwf

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: tem1, tem2, tem3

      integer :: i, k, id, nf

!
!===> ... begin here
!
      do nf=1,nf_clds
        do k=1,nlay
          do i=1,ix
            clouds(i,k,nf) = 0.0
          enddo
        enddo
      enddo
!     clouds(:,:,:) = 0.0

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          cwp   (i,k) = 0.0
          cip   (i,k) = 0.0
          crp   (i,k) = 0.0
          csp   (i,k) = 0.0
          rew   (i,k) = reliq_def            ! default liq radius to 10 micron
          rei   (i,k) = reice_def            ! default ice radius to 50 micron
          rer   (i,k) = rrain_def            ! default rain radius to 1000 micron
          res   (i,k) = rsnow_def            ! default snow radius to 250 micron
          tem2d (i,k) = min( 1.0, max( 0.0, (con_ttp-tlyr(i,k))*0.05 ) )
!         clwf(i,k)   = 0.0
        enddo
      enddo
!
!     if ( lcrick ) then
!       do i = 1, IX
!         clwf(i,1)    = 0.75*clw(i,1)    + 0.25*clw(i,2)
!         clwf(i,nlay) = 0.75*clw(i,nlay) + 0.25*clw(i,nlay-1)
!       enddo
!       do k = 2, NLAY-1
!         do i = 1, IX
!           clwf(i,K) = 0.25*clw(i,k-1) + 0.5*clw(i,k) + 0.25*clw(i,k+1)
!         enddo
!       enddo
!     else
!       do k = 1, NLAY
!         do i = 1, IX
!           clwf(i,k) = clw(i,k)
!         enddo
!       enddo
!     endif

!> -# Find top pressure for each cloud domain for given latitude.
!     ptopc(k,i): top presure of each cld domain (k=1-4 are sfc,L,m,h;
!  ---  i=1,2 are low-lat (<45 degree) and pole regions)

      do id = 1, 4
        tem1 = ptopc(id,2) - ptopc(id,1)

        do i =1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,id) = ptopc(id,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!> -# Compute cloud liquid/ice condensate path in \f$ g/m^2 \f$ .

      if ( ivflip == 0 ) then          ! input data from toa to sfc
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k+1) - plvl(i,k)
            tem1      = gfac * delp(i,k)
            cip(i,k)  = ciw(i,k) * tem1
            cwp(i,k)  = clw(i,k) * tem1
          enddo
        enddo
      else                             ! input data from sfc to toa
        do k = 1, NLAY
          do i = 1, IX
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
            tem1      = gfac * delp(i,k)
            cip(i,k)  = ciw(i,k) * tem1
            cwp(i,k)  = clw(i,k) * tem1
          enddo
        enddo
      endif                            ! end_if_ivflip

!> -# Compute effective liquid cloud droplet radius over land.

      do i = 1, IX
        if (nint(slmsk(i)) == 1) then
          do k = 1, NLAY
            rew(i,k) = 5.0 + 5.0 * tem2d(i,k)
          enddo
        endif
      enddo
      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = cldcov(i,k)
        enddo
      enddo

      do k = 1, NLAY
        do i = 1, IX
          if (cldtot(i,k) < climit) then
            cldtot(i,k) = 0.0
            cwp(i,k)    = 0.0
            cip(i,k)    = 0.0
            crp(i,k)    = 0.0
            csp(i,k)    = 0.0
          endif
        enddo
      enddo

      if ( lcnorm ) then
        do k = 1, NLAY
          do i = 1, IX
            if (cldtot(i,k) >= climit) then
              tem1 = 1.0 / max(climit2, cldtot(i,k))
              cwp(i,k) = cwp(i,k) * tem1
              cip(i,k) = cip(i,k) * tem1
              crp(i,k) = crp(i,k) * tem1
              csp(i,k) = csp(i,k) * tem1
            endif
          enddo
        enddo
      endif

!> -# Compute effective ice cloud droplet radius following 
!! \cite heymsfield_and_mcfarquhar_1996.

      do k = 1, NLAY
        do i = 1, IX
          tem2 = tlyr(i,k) - con_ttp

          if (cip(i,k) > 0.0) then
            tem3 = gord * cip(i,k) * plyr(i,k) / (delp(i,k)*tvly(i,k))

            if (tem2 < -50.0) then
              rei(i,k) = (1250.0/9.917) * tem3 ** 0.109
            elseif (tem2 < -40.0) then
              rei(i,k) = (1250.0/9.337) * tem3 ** 0.08
            elseif (tem2 < -30.0) then
              rei(i,k) = (1250.0/9.208) * tem3 ** 0.055
            else
              rei(i,k) = (1250.0/9.387) * tem3 ** 0.031
            endif
!           rei(i,k)   = max(20.0, min(rei(i,k), 300.0))
!           rei(i,k)   = max(10.0, min(rei(i,k), 100.0))
            rei(i,k)   = max(10.0, min(rei(i,k), 150.0))
!           rei(i,k)   = max(5.0,  min(rei(i,k), 130.0))
          endif
        enddo
      enddo

!
      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = cldtot(i,k)
          clouds(i,k,2) = cwp(i,k)
          clouds(i,k,3) = rew(i,k)
          clouds(i,k,4) = cip(i,k)
          clouds(i,k,5) = rei(i,k)
!         clouds(i,k,6) = 0.0
          clouds(i,k,7) = rer(i,k)
!         clouds(i,k,8) = 0.0
          clouds(i,k,9) = rei(i,k)
        enddo
      enddo


!> -# Call gethml() to compute low,mid,high,total, and boundary layer
!!    cloud fractions and clouds top/bottom layer indices for low, mid,
!!    and high clouds.
!  ---  compute low, mid, high, total, and boundary layer cloud fractions
!       and clouds top/bottom layer indices for low, mid, and high clouds.
!       The three cloud domain boundaries are defined by ptopc.  The cloud
!       overlapping method is defined by control flag 'iovr', which may
!       be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX,NLAY,                                                   &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )


!
      return
!...................................
      end subroutine progclduni
!-----------------------------------
!> @}

!> \ingroup module_radiation_clouds
!> This subroutine computes cloud fractions for radiation calculations.
!!\param plyr      (IX,NLAY), model layer mean pressure in mb (100Pa)
!!\param plvl      (IX,NLP1), model level pressure in mb (100Pa)
!!\param tlyr      (IX,NLAY), model layer mean temperature in K
!!\param rhly      (IX,NLAY), layer relative humidity
!!\param vvel      (IX,NLAY), layer mean vertical velocity in mb/sec
!!\param cv        (IX), fraction of convective cloud
!!\param cvt, cvb  (IX), conv cloud top/bottom pressure in mb
!!\param xlat      (IX), grid latitude in radians, default to pi/2 ->
!!                 -pi/2 range, otherwise see in-line comment
!!\param xlon      (IX), grid longitude in radians, ok for both 0->2pi
!!                 or -pi -> +pi ranges
!!\param slmsk     (IX), sea/land mask array (sea:0,land:1,sea-ice:2)
!!\param IX        horizontal dimention
!!\param NLAY,NLP1       vertical layer/level dimensions
!!\param clouds    (IX,NLAY,NF_CLDS), cloud profiles
!!\n               (:,:,1) - layer total cloud fraction
!!\n               (:,:,2) - layer cloud optical depth
!!\n               (:,:,3) - layer cloud single scattering albedo
!!\n               (:,:,4) - layer cloud asymmetry factor
!!\param clds      (IX,5), fraction of clouds for low, mid, hi, tot, bl
!!\param mtop      (IX,3), vertical indices for low, mid, hi cloud tops
!!\param mbot      (IX,3), vertical indices for low, mid, hi cloud bases
!!\section gen_diagcld1 diagcld1 General Algorithm
!> @{
      subroutine diagcld1                                               &
     &     ( plyr,plvl,tlyr,rhly,vvel,cv,cvt,cvb,                       &    !  ---  inputs:
     &       xlat,xlon,slmsk,                                           &
     &       IX, NLAY, NLP1,                                            &
     &       clouds,clds,mtop,mbot                                      &    !  ---  outputs:
     &      )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:    diagcld1    computes cloud fractions for radiation     !
!   calculations.                                                       !
!                                                                       !
! abstract:  clouds are diagnosed from layer relative humidity, and     !
!   estimate cloud optical depth from temperature and layer thickness.  !
!   then computes the low, mid, high, total and boundary layer cloud    !
!   fractions and the vertical indices of low, mid, and high cloud top  !
!   and base.  the three vertical cloud domains are set up in the       !
!   initial subroutine "cld_init".                                      !
!                                                                       !
! usage:         call diagcld1                                          !
!                                                                       !
! subprograms called:   gethml                                          !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IX,NLP1) : model level pressure in mb (100Pa)                !
!   tlyr  (IX,NLAY) : model layer mean temperature in k                 !
!   rhly  (IX,NLAY) : layer relative humidity                           !
!   vvel  (IX,NLAY) : layer mean vertical velocity in mb/sec            !
!   clw   (IX,NLAY) : layer cloud condensate amount         (not used)  !
!   xlat  (IX)      : grid latitude in radians, default to pi/2 -> -pi/2!
!                     range, otherwise see in-line comment              !
!   xlon  (IX)      : grid longitude in radians, ok for both 0->2pi or  !
!                     -pi -> +pi ranges                                 !
!   slmsk (IX)      : sea/land mask array (sea:0,land:1,sea-ice:2)      !
!   cv    (IX)      : fraction of convective cloud                      !
!   cvt, cvb (IX)   : conv cloud top/bottom pressure in mb              !
!   IX              : horizontal dimention                              !
!   NLAY,NLP1       : vertical layer/level dimensions                   !
!                                                                       !
! output variables:                                                     !
!   clouds(IX,NLAY,NF_CLDS) : cloud profiles                            !
!      clouds(:,:,1) - layer total cloud fraction                       !
!      clouds(:,:,2) - layer cloud optical depth                        !
!      clouds(:,:,3) - layer cloud single scattering albedo             !
!      clouds(:,:,4) - layer cloud asymmetry factor                     !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! external module variables:                                            !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer,  intent(in) :: IX, NLAY, NLP1

      real (kind=kind_phys), dimension(:,:), intent(in) :: plvl, plyr,  &
     &       tlyr, rhly, vvel

      real (kind=kind_phys), dimension(:),   intent(in) :: xlat, xlon,  &
     &       slmsk, cv, cvt, cvb

!  ---  outputs
      real (kind=kind_phys), dimension(:,:,:), intent(out) :: clouds

      real (kind=kind_phys), dimension(:,:),   intent(out) :: clds

      integer,               dimension(:,:),   intent(out) :: mtop,mbot

!  ---  local variables:
      real (kind=kind_phys), dimension(IX,NLAY) :: cldtot, cldcnv,      &
     &       cldtau, taufac, dthdp, tem2d

      real (kind=kind_phys) :: ptop1(IX,NK_CLDS+1)

      real (kind=kind_phys) :: cr1, cr2, crk, pval, cval, omeg, value,  &
     &       tem1, tem2

      integer, dimension(IX):: idom, kcut

!  ---  for rh-cl calculation
      real (kind=kind_phys) :: xlatdg, xlondg, xlnn, xlss, xrgt, xlft,  &
     &       rhcla(NBIN,NLON,MCLD,NSEAL),  rhcld(IX,NBIN,MCLD)

      integer :: ireg, ib, ic, id, id1, il, is, nhalf

      integer :: i, j, k, klowt
!     integer :: klowb

      logical :: notstop

!
!===> ... begin here
!
      clouds(:,:,:) = 0.0

      tem1 = 180.0 / con_pi

      lab_do_i_IX : do i = 1, IX

        xlatdg = xlat(i) * tem1                    ! if xlat in pi/2 -> -pi/2 range
!       xlatdg = 90.0 - xlat(i)*tem1               ! if xlat in 0 -> pi range

        xlondg = xlon(i) * tem1
        if (xlondg < 0.0) xlondg = xlondg + 360.0  ! if in -180->180, chg to 0->360 range

        ireg = 4

!> -# Get rh-cld relation for this lat.

        lab_do_j : do j = 1, 3
          if (xlatdg > xlabdy(j)) then
            ireg = j
            exit lab_do_j
          endif
        enddo  lab_do_j

        do is = 1, NSEAL
          do ic = 1, MCLD
            do il = 1, NLON
              do ib = 1, NBIN
                rhcla(ib,il,ic,is) = rhcl(ib,il,ireg,ic,is)
              enddo
            enddo
          enddo
        enddo

!> -# Linear transition between latitudinal regions.
        do j = 1, 3
          xlnn = xlabdy(j) + xlim
          xlss = xlabdy(j) - xlim

          if (xlatdg < xlnn .and. xlatdg > xlss) then
            do is = 1, NSEAL
              do ic = 1, MCLD
                do il = 1, NLON
                  do ib = 1, NBIN
                    rhcla(ib,il,ic,is) = rhcl(ib,il,j+1,ic,is)          &
     &                + (rhcl(ib,il,j,ic,is)-rhcl(ib,il,j+1,ic,is))     &
     &                * (xlatdg-xlss) / (xlnn-xlss)
                  enddo
                enddo
              enddo
            enddo
          endif

        enddo        ! end_j_loop

!> -# Get rh-cld relationship for each grid point, interpolating
!!    longitudinally between regions if necessary.

        if (slmsk(i) < 1.0) then
          is = 2
        else
          is = 1
        endif

!  ---  which hemisphere (e,w)

        if (xlondg > 180.e0) then
          il = 2
        else
          il = 1
        endif

        do ic = 1, MCLD
          do ib = 1, NBIN
            rhcld(i,ib,ic) = rhcla(ib,il,ic,is)
          enddo

          lab_do_k : do k = 1, 3
            tem2 = abs(xlondg - xlobdy(k))

            if (tem2 < xlim) then
              id = il
              id1= id + 1
              if (id1 > NLON) id1 = 1

              xlft = xlobdy(k) - xlim
              xrgt = xlobdy(k) + xlim

              do ib = 1, NBIN
                rhcld(i,ib,ic) = rhcla(ib,id1,ic,is)                    &
     &            + (rhcla(ib,id,ic,is) - rhcla(ib,id1,ic,is))          &
     &            * (xlondg-xrgt)/(xlft-xrgt)
              enddo
              exit lab_do_k
            endif

          enddo  lab_do_k

        enddo   ! end_do_ic_loop
      enddo  lab_do_i_IX

!> -# Find top pressure for each cloud domain.

      do j = 1, 4
        tem1 = ptopc(j,2) - ptopc(j,1)

        do i = 1, IX
          tem2 = xlat(i) / con_pi        ! if xlat in pi/2 -> -pi/2 range
!         tem2 = 0.5 - xlat(i)/con_pi    ! if xlat in 0 -> pi range

          ptop1(i,j) = ptopc(j,1) + tem1*max( 0.0, 4.0*abs(tem2)-1.0 )
        enddo
      enddo

!> -# Compute stratiform cloud optical depth.

      do k = 1, NLAY
        do i = 1, IX
          tem1 = tlyr(i,k) - con_ttp
          if (tem1 <= -10.0) then
            cldtau(i,k) = max( 0.1e-3, 2.0e-6*(tem1+82.5)**2 )
          else
            cldtau(i,k) = min( 0.08, 6.949e-3*tem1+0.08 )
          endif
        enddo
      enddo

!> -# Calculate potential temperature and its lapse rate.

      do k = 1, NLAY
        do i = 1, IX
          cldtot(i,k) = 0.0
          cldcnv(i,k) = 0.0
          tem1        = (plyr(i,k)*0.001) ** (-con_rocp)
          tem2d(i,k)  = tem1 * tlyr(i,k)
        enddo
      enddo

      do k = 1, NLAY-1
        do i = 1, IX
          dthdp(i,k) = (tem2d(i,k+1)-tem2d(i,k))/(plyr(i,k+1)-plyr(i,k))
        enddo
      enddo
!
!> -# Diagnostic method to find cloud amount cldtot, cldcnv.
!

      if ( ivflip == 0 ) then                 ! input data from toa to sfc

!>    - Find the lowest low cloud top sigma level, computed for each
!!      lat cause domain definition changes with latitude.

!       klowb = 1
        klowt = 1
        do k = 1, NLAY
          do i = 1, IX
!           if (plvl(i,k) < ptop1(i,2))  klowb = k
            if (plvl(i,k) < ptop1(i,2))  klowt = max(klowt,k)
            taufac(i,k) = plvl(i,k+1) - plvl(i,k)
          enddo
        enddo

        do i = 1, IX

!>    - Find the stratosphere cut off layer for high cloud (about
!!      250mb). It is assumed to be above the layerwith dthdp less than
!!      -0.25 in the high cloud domain.

          kcut(i) = 2
          lab_do_kcut0 : do k = klowt-1, 2, -1
            if (plyr(i,k) <= ptop1(i,3) .and.                           &
     &          dthdp(i,k) < -0.25e0) then
              kcut(i) = k
              exit lab_do_kcut0
            endif
          enddo  lab_do_kcut0

!>    - Put convective cloud into 'cldcnv', no merge at this point.

          if (cv(i) >= climit .and. cvt(i) < cvb(i)) then
            id  = NLAY
            id1 = NLAY

            lab_do_k_cvt0 : do k = 2, NLAY
              if (cvt(i) <= plyr(i,k)) then
                id = k - 1
                exit lab_do_k_cvt0
              endif
            enddo  lab_do_k_cvt0

            lab_do_k_cvb0 : do k = NLAY-1, 1, -1
              if (cvb(i) >= plyr(i,k)) then
                id1 = k + 1
                exit lab_do_k_cvb0
              endif
            enddo  lab_do_k_cvb0

            tem1 = plyr(i,id1) - plyr(i,id)
            do k = id, id1
              cldcnv(i,k) = cv(i)
              taufac(i,k) = taufac(i,k) * max( 0.25, 1.0-0.125*tem1 )
              cldtau(i,k) = 0.06
            enddo
          endif

        enddo                ! end_do_i_loop

!>    - Calculate stratiform cloud and put into array 'cldtot' using
!!      the cld-rh relationship from table look-up, where tables
!!      obtained using k.mitchell frequency distribution tuning.
!bl       (observations are daily means from us af rtneph).....k.a.c.
!bl       tables created without lowest 10 percent of atmos.....k.a.c.
!      (observations are synoptic using -6,+3 window from rtneph)
!       tables are created with lowest 10-percent-of-atmos, and are
!  ---  now used..  25 october 1995 ... kac.

        do k = NLAY-1, 2, -1

          if (k < llyr) then
            do i = 1, IX
              idom(i) = 0
            enddo

            do i = 1, IX
              lab_do_ic0 : do ic = 2, 4
                if(plyr(i,k) >= ptop1(i,ic)) then
                  idom(i) = ic
                  exit lab_do_ic0
                endif
              enddo  lab_do_ic0
            enddo
          else
            do i = 1, IX
              idom(i) = 1
            enddo
          endif

          do i = 1, IX
            id = idom(i)
            nhalf = (NBIN + 1) / 2

            if (id <= 0 .or. k < kcut(i)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) <= rhcld(i,1,id)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) >= rhcld(i,NBIN,id)) then
              cldtot(i,k) = 1.0
            else
              ib = nhalf
              crk = rhly(i,k)

              notstop = .true.
              do while ( notstop )
                nhalf = (nhalf + 1) / 2
                cr1 = rhcld(i,ib,  id)
                cr2 = rhcld(i,ib+1,id)

                if (crk <= cr1) then
                  ib = max( ib-nhalf, 1 )
                elseif (crk > cr2) then
                  ib = min( ib+nhalf, NBIN-1 )
                else
                  cldtot(i,k) = 0.01 * (ib + (crk - cr1)/(cr2 - cr1))
                  notstop = .false.
                endif
              enddo      ! end_do_while
            endif
          enddo          ! end_do_i_loop

        enddo            ! end_do_k_loop

!>    - Compute vertical velocity adjustment on low clouds.

        value = vvcld1 - vvcld2
        do k = klowt, llyr+1
          do i = 1, IX

            omeg = vvel(i,k)
            cval = cldtot(i,k)
            pval = plyr(i,k)

! --- vertical velocity adjustment on low clouds

            if (cval >= climit .and. pval >= ptop1(i,2)) then
              if (omeg >= vvcld1) then
                cldtot(i,k) = 0.0
              elseif (omeg > vvcld2) then
                tem1 = (vvcld1 - omeg) / value
                cldtot(i,k) = cldtot(i,k) * sqrt(tem1)
              endif
            endif

          enddo     ! end_do_i_loop
        enddo       ! end_do_k_loop

      else                                    ! input data from sfc to toa

!  ---  find the lowest low cloud top sigma level, computed for each lat cause
!       domain definition changes with latitude...

!       klowb = NLAY
        klowt = NLAY
        do k = NLAY, 1, -1
          do i = 1, IX
!           if (plvl(i,k) < ptop1(i,2))  klowb = k
            if (plvl(i,k) < ptop1(i,2))  klowt = min(klowt,k)
            taufac(i,k) = plvl(i,k) - plvl(i,k+1)       ! dp for later cal cldtau use
          enddo
        enddo

        do i = 1, IX

!  ---  find the stratosphere cut off layer for high cloud (about 250mb).
!       it is assumed to be above the layer with dthdp less than -0.25 in
!       the high cloud domain

          kcut(i) = NLAY - 1
          lab_do_kcut1 : do k = klowt+1, NLAY-1
            if (plyr(i,k) <= ptop1(i,3) .and.                           &
     &          dthdp(i,k) < -0.25e0) then
              kcut(i) = k
              exit lab_do_kcut1
            endif
          enddo  lab_do_kcut1

!  ---  put convective cloud into 'cldcnv', no merge at this point..

          if (cv(i) >= climit .and. cvt(i) < cvb(i)) then
            id  = 1
            id1 = 1

            lab_do_k_cvt : do k = NLAY-1, 1, -1
              if (cvt(i) <= plyr(i,k)) then
                id = k + 1
                exit lab_do_k_cvt
              endif
            enddo  lab_do_k_cvt

            lab_do_k_cvb : do k = 2, NLAY
              if (cvb(i) >= plyr(i,k)) then
                id1 = k - 1
                exit lab_do_k_cvb
              endif
            enddo  lab_do_k_cvb

            tem1 = plyr(i,id1) - plyr(i,id)
            do k = id1, id
              cldcnv(i,k) = cv(i)
              taufac(i,k) = taufac(i,k) * max( 0.25, 1.0-0.125*tem1 )
              cldtau(i,k) = 0.06
            enddo
          endif

        enddo     ! end_do_i_loop

!  ---  calculate stratiform cloud and put into array 'cldtot' using
!       the cloud-rel.humidity relationship from table look-up..where
!       tables obtained using k.mitchell frequency distribution tuning
!bl       (observations are daily means from us af rtneph).....k.a.c.
!bl       tables created without lowest 10 percent of atmos.....k.a.c.
!      (observations are synoptic using -6,+3 window from rtneph)
!       tables are created with lowest 10-percent-of-atmos, and are
!  ---  now used..  25 october 1995 ... kac.

        do k = 2, NLAY-1

          if (k > llyr) then
            do i = 1, IX
              idom(i) = 0
            enddo

            do i = 1, IX
              lab_do_ic1 : do ic = 2, 4
                if(plyr(i,k) >= ptop1(i,ic)) then
                  idom(i) = ic
                  exit lab_do_ic1
                endif
              enddo  lab_do_ic1
            enddo
          else
            do i = 1, IX
              idom(i) = 1
            enddo
          endif

          do i = 1, IX
            id = idom(i)
            nhalf = (NBIN + 1) / 2

            if (id <= 0 .or. k > kcut(i)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) <= rhcld(i,1,id)) then
              cldtot(i,k) = 0.0
            elseif (rhly(i,k) >= rhcld(i,NBIN,id)) then
              cldtot(i,k) = 1.0
            else
              ib = nhalf
              crk = rhly(i,k)

              notstop = .true.
              do while ( notstop )
                nhalf = (nhalf + 1) / 2
                cr1 = rhcld(i,ib,  id)
                cr2 = rhcld(i,ib+1,id)

                if (crk <= cr1) then
                  ib = max( ib-nhalf, 1 )
                elseif (crk > cr2) then
                  ib = min( ib+nhalf, NBIN-1 )
                else
                  cldtot(i,k) = 0.01 * (ib + (crk - cr1)/(cr2 - cr1))
                  notstop = .false.
                endif
              enddo      ! end_do_while
            endif
          enddo          ! end_do_i_loop

        enddo            ! end_do_k_loop

! --- vertical velocity adjustment on low clouds

        value = vvcld1 - vvcld2
        do k = llyr-1, klowt
          do i = 1, IX

            omeg = vvel(i,k)
            cval = cldtot(i,k)
            pval = plyr(i,k)

! --- vertical velocity adjustment on low clouds

            if (cval >= climit .and. pval >= ptop1(i,2)) then
              if (omeg >= vvcld1) then
                cldtot(i,k) = 0.0
              elseif (omeg > vvcld2) then
                tem1 = (vvcld1 - omeg) / value
                cldtot(i,k) = cldtot(i,k) * sqrt(tem1)
              endif
            endif

          enddo     ! end_do_i_loop
        enddo       ! end_do_k_loop

      endif                                   ! end_if_ivflip

!>    - Calculate diagnostic cloud optical depth.

!     cldtau = cldtau * taufac

      where (cldtot < climit)
        cldtot = 0.0
      endwhere
      where (cldcnv < climit)
        cldcnv = 0.0
      endwhere

      where (cldtot < climit .and. cldcnv < climit)
        cldtau = 0.0
      endwhere

      do k = 1, NLAY
        do i = 1, IX
          clouds(i,k,1) = max(cldtot(i,k), cldcnv(i,k))
          clouds(i,k,2) = cldtau(i,k) * taufac(i,k)
          clouds(i,k,3) = cldssa_def
          clouds(i,k,4) = cldasy_def
        enddo
      enddo

!> -# Call gethml(), to compute low, mid, high, total, and boundary
!!    layer cloud fractions and cloud top/bottom layer indices for low,
!!    mid, and high clouds.
!         the three cloud domain boundaries are defined by ptopc.  the cloud
!         overlapping method is defined by control flag 'iovr', which may
!         be different for lw and sw radiation programs.

      call gethml                                                       &
!  ---  inputs:
     &     ( plyr, ptop1, cldtot, cldcnv,                               &
     &       IX, NLAY,                                                  &
!  ---  outputs:
     &       clds, mtop, mbot                                           &
     &     )

!
      return
!...................................
      end subroutine diagcld1
!-----------------------------------
!> @}

!> \ingroup module_radiation_clouds
!> This subroutine computes high, mid, low, total, and boundary cloud
!! fractions and cloud top/bottom layer indices for model diagnostic
!! output. The three cloud domain boundaries are defined by ptopc. The
!! cloud overlapping method is defined by control flag 'iovr', which is
!! also used by LW and SW radiation programs.
!> \param plyr    (IX,NLAY), model layer mean pressure in mb (100Pa)
!> \param ptop1   (IX,4), pressure limits of cloud domain interfaces
!!                    (sfc,low,mid,high) in mb (100Pa)
!> \param cldtot  (IX,NLAY), total or stratiform cloud profile in fraction
!> \param cldcnv  (IX,NLAY), convective cloud (for diagnostic scheme only)
!> \param IX      horizontal dimension
!> \param NLAY    vertical layer dimensions
!> \param clds   (IX,5), fraction of clouds for low, mid, hi, tot, bl
!> \param mtop   (IX,3),vertical indices for low, mid, hi cloud tops
!> \param mbot   (IX,3),vertical indices for low, mid, hi cloud bases
!!
!>\section detail gethml General Algorithm
!! @{
      subroutine gethml                                                 &
     &     ( plyr, ptop1, cldtot, cldcnv,                               &       !  ---  inputs:
     &       IX, NLAY,                                                  &
     &       clds, mtop, mbot                                           &       !  ---  outputs:
     &     )

!  ===================================================================  !
!                                                                       !
! abstract: compute high, mid, low, total, and boundary cloud fractions !
!   and cloud top/bottom layer indices for model diagnostic output.     !
!   the three cloud domain boundaries are defined by ptopc.  the cloud  !
!   overlapping method is defined by control flag 'iovr', which is also !
!   used by lw and sw radiation programs.                               !
!                                                                       !
! usage:         call gethml                                            !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
! attributes:                                                           !
!   language:   fortran 90                                              !
!   machine:    ibm-sp, sgi                                             !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IX,NLAY) : model layer mean pressure in mb (100Pa)           !
!   ptop1 (IX,4)    : pressure limits of cloud domain interfaces        !
!                     (sfc,low,mid,high) in mb (100Pa)                  !
!   cldtot(IX,NLAY) : total or straiform cloud profile in fraction      !
!   cldcnv(IX,NLAY) : convective cloud (for diagnostic scheme only)     !
!   IX              : horizontal dimention                              !
!   NLAY            : vertical layer dimensions                         !
!                                                                       !
! output variables:                                                     !
!   clds  (IX,5)    : fraction of clouds for low, mid, hi, tot, bl      !
!   mtop  (IX,3)    : vertical indices for low, mid, hi cloud tops      !
!   mbot  (IX,3)    : vertical indices for low, mid, hi cloud bases     !
!                                                                       !
! external module variables:  (in physparam)                            !
!   ivflip          : control flag of vertical index direction          !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!                                                                       !
! internal module variables:                                            !
!   iovr            : control flag for cloud overlap                    !
!                     =0 random overlapping clouds                      !
!                     =1 max/ran overlapping clouds                     !
!                                                                       !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none!

!  ---  inputs:
      integer, intent(in) :: IX, NLAY

      real (kind=kind_phys), dimension(:,:), intent(in) :: plyr, ptop1, &
     &       cldtot, cldcnv

!  ---  outputs
      real (kind=kind_phys), dimension(:,:), intent(out) :: clds

      integer,               dimension(:,:), intent(out) :: mtop, mbot

!  ---  local variables:
      real (kind=kind_phys) :: cl1(IX), cl2(IX)
      real (kind=kind_phys) :: pcur, pnxt, ccur, cnxt

      integer, dimension(IX):: idom, kbt1, kth1, kbt2, kth2
      integer :: i, k, id, id1, kstr, kend, kinc

!
!===> ... begin here
!
      clds(:,:) = 0.0

      do i = 1, IX
        cl1(i) = 1.0
        cl2(i) = 1.0
      enddo

!  ---  total and bl clouds, where cl1, cl2 are fractions of clear-sky view
!       layer processed from surface and up

!> -# Calculate total and BL cloud fractions (maximum-random cloud
!!    overlapping is operational).

      if ( ivflip == 0 ) then                   ! input data from toa to sfc
        kstr = NLAY
        kend = 1
        kinc = -1
      else                                      ! input data from sfc to toa
        kstr = 1
        kend = NLAY
        kinc = 1
      endif                                     ! end_if_ivflip

      if ( iovr == 0 ) then                     ! random overlap

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) cl1(i) = cl1(i) * (1.0 - ccur)
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = 1.0 - cl1(i)          ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = 1.0 - cl1(i)              ! save total cloud
        enddo

      else                                      ! max/ran overlap

        do k = kstr, kend, kinc
          do i = 1, IX
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))
            if (ccur >= climit) then             ! cloudy layer
              cl2(i) = min( cl2(i), (1.0 - ccur) )
            else                                ! clear layer
              cl1(i) = cl1(i) * cl2(i)
              cl2(i) = 1.0
            endif
          enddo

          if (k == llyr) then
            do i = 1, IX
              clds(i,5) = 1.0 - cl1(i) * cl2(i) ! save bl cloud
            enddo
          endif
        enddo

        do i = 1, IX
          clds(i,4) = 1.0 - cl1(i) * cl2(i)     ! save total cloud
        enddo

      endif                                     ! end_if_iovr

!  ---  high, mid, low clouds, where cl1, cl2 are cloud fractions
!       layer processed from one layer below llyr and up
!  ---  change! layer processed from surface to top, so low clouds will
!       contains both bl and low clouds.

!> -# Calculte high, mid, low cloud fractions and vertical indices of
!!    cloud tops/bases.
      if ( ivflip == 0 ) then                   ! input data from toa to sfc

        do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = NLAY
          kbt2(i) = NLAY
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = NLAY
          mtop(i,1) = NLAY
          mbot(i,2) = NLAY - 1
          mtop(i,2) = NLAY - 1
          mbot(i,3) = NLAY - 1
          mtop(i,3) = NLAY - 1
        enddo

!org    do k = llyr-1, 1, -1
        do k = NLAY, 1, -1
          do i = 1, IX
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k > 1) then
              pnxt = plyr(i,k-1)
              cnxt = min( ovcst, max( cldtot(i,k-1), cldcnv(i,k-1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if ( iovr == 0 ) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i) )      &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i) )      &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = k - 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif   ! end_if_cnxt_or_pnxt
            endif     ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = min( kbt1(i), kbt1(i)-kth1(i)+1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = k - 1
              kth1(i) = 0

              if (id1 <= NK_CLDS) then
                mbot(i,id1) = kbt1(i)
                mtop(i,id1) = kbt1(i)
              endif
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      else                                      ! input data from sfc to toa

        do i = 1, IX
          cl1 (i) = 0.0
          cl2 (i) = 0.0
          kbt1(i) = 1
          kbt2(i) = 1
          kth1(i) = 0
          kth2(i) = 0
          idom(i) = 1
          mbot(i,1) = 1
          mtop(i,1) = 1
          mbot(i,2) = 2
          mtop(i,2) = 2
          mbot(i,3) = 2
          mtop(i,3) = 2
        enddo

!org    do k = llyr+1, NLAY
        do k = 1, NLAY
          do i = 1, IX
            id = idom(i)
            id1= id + 1

            pcur = plyr(i,k)
            ccur = min( ovcst, max( cldtot(i,k), cldcnv(i,k) ))

            if (k < NLAY) then
              pnxt = plyr(i,k+1)
              cnxt = min( ovcst, max( cldtot(i,k+1), cldcnv(i,k+1) ))
            else
              pnxt = -1.0
              cnxt = 0.0
            endif

            if (pcur < ptop1(i,id1)) then
              id = id + 1
              id1= id1 + 1
              idom(i) = id
            endif

            if (ccur >= climit) then
              if (kth2(i) == 0) kbt2(i) = k
              kth2(i) = kth2(i) + 1

              if ( iovr == 0 ) then
                cl2(i) = cl2(i) + ccur - cl2(i)*ccur
              else
                cl2(i) = max( cl2(i), ccur )
              endif

              if (cnxt < climit .or. pnxt < ptop1(i,id1)) then
                kbt1(i) = nint( (cl1(i)*kbt1(i) + cl2(i)*kbt2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                kth1(i) = nint( (cl1(i)*kth1(i) + cl2(i)*kth2(i))       &
     &                  / (cl1(i) + cl2(i)) )
                cl1 (i) = cl1(i) + cl2(i) - cl1(i)*cl2(i)

                kbt2(i) = k + 1
                kth2(i) = 0
                cl2 (i) = 0.0
              endif     ! end_if_cnxt_or_pnxt
            endif       ! end_if_ccur

            if (pnxt < ptop1(i,id1)) then
              clds(i,id) = cl1(i)
              mtop(i,id) = max( kbt1(i), kbt1(i)+kth1(i)-1 )
              mbot(i,id) = kbt1(i)

              cl1 (i) = 0.0
              kbt1(i) = min(k+1, nlay)
              kth1(i) = 0

              if (id1 <= NK_CLDS) then
                mbot(i,id1) = kbt1(i)
                mtop(i,id1) = kbt1(i)
              endif
            endif     ! end_if_pnxt

          enddo       ! end_do_i_loop
        enddo         ! end_do_k_loop

      endif                                     ! end_if_ivflip

!
      return
!...................................
      end subroutine gethml
!-----------------------------------
!! @}

!> \ingroup module_radiation_clouds
!> cld-rh relations obtained from mitchell-hahn procedure.
!\section rhtable_gen rhtable General Algorithm
      subroutine rhtable                                                &
     &     ( me                                                         &         !  ---  inputs:
     &,      ier )                                                                !  ---  outputs:

!  ===================================================================  !
!                                                                       !
! abstract: cld-rh relations obtained from mitchell-hahn procedure,     !
!   here read cld/rh tuning tables for day 0,1,...,5 and merge into 1   !
!   file.                                                               !
!                                                                       !
! inputs:                                                               !
!   me              : check print control flag                          !
!                                                                       !
! outputs:                                                              !
!   ier             : error flag                                        !
!                                                                       !
!  ===================================================================  !
!
      implicit none!

!  ---  inputs:
      integer, intent(in) :: me

!  ---  output:
      integer, intent(out) :: ier

!  ---  locals:
      real (kind=kind_phys), dimension(NBIN,NLON,NLAT,MCLD,NSEAL) ::    &
     &      rhfd, rtnfd, rhcf, rtncf, rhcla

      real (kind=kind_io4), dimension(NBIN,NLON,NLAT,MCLD,NSEAL) ::     &
     &      rhfd4, rtnfd4

      real(kind=kind_io4)  :: fhour

      real(kind=kind_phys) :: binscl, cfrac, clsat, rhsat, cstem

      integer, dimension(NLON,NLAT,MCLD,NSEAL) :: kpts, kkpts

      integer :: icdays(15), idate(4), nbdayi, isat

      integer :: i, i1, j, k, l, m, id, im, iy

!
!===> ...  begin here
!

      ier = 1

      rewind NICLTUN

      binscl = 1.0 / NBIN

!  ---  array initializations

      do m=1,NSEAL
       do l=1,MCLD
        do k=1,NLAT
         do j=1,NLON
          do i=1,NBIN
            rhcf (i,j,k,l,m) = 0.0
            rtncf(i,j,k,l,m) = 0.0
            rhcla(i,j,k,l,m) = -0.1
          enddo
         enddo
        enddo
       enddo
      enddo

      kkpts = 0

!  ---  read the data off the rotating file

      read (NICLTUN,ERR=998,END=999) nbdayi, icdays

      if (me == 0) print 11, nbdayi
  11  format('   from rhtable DAYS ON FILE =',i5)

      do i = 1, nbdayi
       id = icdays(i) / 10000
       im = (icdays(i)-id*10000) / 100
       iy = icdays(i)-id*10000-im*100
       if (me == 0) print 51, id,im,iy
  51   format('   from rhtable ARCHV DATA FROM DA,MO,YR=',3i4)
      enddo

      read (NICLTUN,ERR=998,END=999) fhour,idate

      do i1 = 1, nbdayi
        read (NICLTUN) rhfd4
        rhfd = rhfd4

        read (NICLTUN) rtnfd4
        rtnfd = rtnfd4

        read (NICLTUN) kpts

        do m=1,NSEAL
         do l=1,MCLD
          do k=1,NLAT
           do j=1,NLON
            do i=1,NBIN
              rhcf (i,j,k,l,m) = rhcf (i,j,k,l,m) + rhfd (i,j,k,l,m)
              rtncf(i,j,k,l,m) = rtncf(i,j,k,l,m) + rtnfd(i,j,k,l,m)
            enddo
           enddo
          enddo
         enddo
        enddo

        kkpts = kkpts + kpts

      enddo     ! end_do_i1_loop

      do m = 1, NSEAL
       do l = 1, MCLD
        do k = 1, NLAT
         do j = 1, NLON

!  ---  compute the cumulative frequency distribution

           do i = 2, NBIN
             rhcf (i,j,k,l,m) = rhcf (i-1,j,k,l,m) + rhcf (i,j,k,l,m)
             rtncf(i,j,k,l,m) = rtncf(i-1,j,k,l,m) + rtncf(i,j,k,l,m)
           enddo   ! end_do_i_loop

           if (kkpts(j,k,l,m) > 0) then
             do i = 1, NBIN
               rhcf (i,j,k,l,m)= rhcf (i,j,k,l,m)/kkpts(j,k,l,m)
               rtncf(i,j,k,l,m)=min(1., rtncf(i,j,k,l,m)/kkpts(j,k,l,m))
             enddo

!  ---  cause we mix calculations of rh retune with cray and ibm words
!       the last value of rhcf is close to but ne 1.0,
!  ---  so we reset it in order that the 360 loop gives complete tabl

             rhcf(NBIN,j,k,l,m) = 1.0

             do i = 1, NBIN
               lab_do_i1 : do i1 = 1, NBIN
                 if (rhcf(i1,j,k,l,m) >= rtncf(i,j,k,l,m)) then
                   rhcla(i,j,k,l,m) = i1 * binscl
                   exit  lab_do_i1
                 endif
               enddo  lab_do_i1
             enddo

           else                   ! if_kkpts
!  ---  no critical rh

             do i = 1, NBIN
               rhcf (i,j,k,l,m) = -0.1
               rtncf(i,j,k,l,m) = -0.1
             enddo

             if (me == 0) then
               print 210, k,j,m
 210           format('  NO CRIT RH FOR LAT=',I3,' AND LON BAND=',I3,   &
     &                ' LAND(=1) SEA=',I3//'  MODEL RH ',' OBS RTCLD')
               do i = 1, NBIN
                 print 203, rhcf(i,j,k,l,m), rtncf(i,j,k,l,m)
 203             format(2f10.2)
               enddo
             endif

           endif               ! if_kkpts

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

      do m = 1, NSEAL
       do l = 1, MCLD
        do k = 1, NLAT
         do j = 1, NLON

           isat = 0
           do i = 1, NBIN-1
             cfrac = binscl * (i - 1)

             if (rhcla(i,j,k,l,m) < 0.0) then
               print 1941, i,m,l,k,j
 1941          format('  NEG RHCLA FOR IT,NSL,NC,LAT,LON=',5I4          &
     &,               '...STOPPP..')
               stop
             endif

             if (rtncf(i,j,k,l,m) >= 1.0) then
               if (isat <= 0) then
                 isat  = i
                 rhsat = rhcla(i,j,k,l,m)
                 clsat = cfrac
               endif

               rhcla(i,j,k,l,m) = rhsat + (1.0 - rhsat)                 &
     &                         * (cfrac - clsat) / (1.0 - clsat)
             endif
           enddo

           rhcla(NBIN,j,k,l,m) = 1.0

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

!  ---  smooth out the table as it reaches rh=1.0, via linear interpolation
!       between location of rh ge .98 and the NBIN bin (where rh=1.0)
!       previously rh=1.0 occurred for many of the latter bins in the
!  ---  table, thereby giving a cloud value of less then 1.0 for rh=1.0

      rhcl = rhcla

      do m = 1, NSEAL
       do l = 1, MCLD
        do k = 1, NLAT
         do j = 1, NLON

           lab_do_i : do i = 1, NBIN - 2
             cfrac = binscl * i

             if (rhcla(i,j,k,l,m) >= 0.98) then
               do i1 = i, NBIN
                 cstem = binscl * i1

                 rhcl(i1,j,k,l,m) = rhcla(i,j,k,l,m)                    &
     &                    + (rhcla(NBIN,j,k,l,m) - rhcla(i,j,k,l,m))    &
     &                    * (cstem - cfrac) / (1.0 - cfrac)
               enddo
               exit  lab_do_i
             endif
           enddo  lab_do_i

         enddo    ! end_do_j_loop
        enddo     ! end_do_k_loop
       enddo      ! end_do_l_loop
      enddo       ! end_do_m_loop

      if (me == 0) then
        print *,'completed rhtable for cloud tuninig tables'
      endif
      return

 998  print 988
 988  format(' from rhtable ERROR READING TABLES')
      ier = -1
      return

 999  print 989
 989  format(' from rhtable E.O.F READING TABLES')
      ier = -1
      return

!...................................
      end subroutine rhtable
!-----------------------------------


!
!........................................!
      end module module_radiation_clouds !
!========================================!
