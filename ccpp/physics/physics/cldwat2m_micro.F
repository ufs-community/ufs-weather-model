      module cldwat2m_micro

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for microphysics
!
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------

#ifdef NEMS_GSM
       use machine,       only : r8 => kind_phys
       use physcons, gravit => con_g,    rair    => con_rd,             &
     &               rh2o   => con_rv,   epsilon => con_eps,            &
     &               tmelt  => con_tice, cpair   => con_cp,             &
     &               latvap => con_hvap, latice  => con_hfus
       use wv_saturation, only : estblf, hlatv, tmin, hlatf, rgasv, pcf,&
     &                           epsqs, ttrice, vqsatd2,cp,             &
     &                           vqsatd2_single,polysvp,gestbl
       use funcphys,      only : fpvs, fpvsl, fpvsi
#endif

#ifdef GEOS5
       use MAPL_ConstantsMod, r8 => MAPL_R8
       use wv_saturation,    only: estblf, hlatv, tmin, hlatf, rgasv,   &
     &                             pcf, cp, epsqs, ttrice, vqsatd2,     &
     &                             vqsatd2_single,polysvp,gestbl
#endif

#ifdef CAM
       use shr_kind_mod,     only: r8=>shr_kind_r8
       use spmd_utils,       only: masterproc
       use ppgrid,           only: pcols, pver, pverp
       use physconst,        only: gravit, rair, tmelt, cpair, rh2o,
     &                             r_universal, mwh2o, rhoh2o
       use physconst,        only: latvap, latice
       use abortutils,       only: endrun
       use error_function,   only: erf,erfc
       use wv_saturation,    only: estblf, hlatv, tmin, hlatf, rgasv,   &
     &                             pcf, cp, epsqs, ttrice, vqsatd2,     &
     &                             vqsatd2_single,polysvp
       use cam_history,      only: addfld, add_default, phys_decomp,    &
     &                             outfld
       use cam_logfile,      only: iulog
       use rad_constituents, only: rad_cnst_get_clim_info,
     &                             rad_cnst_get_clim_aer_props
       use phys_control,     only: phys_getopts
       use cldwat2m_macro,   only: rhmini, rhmaxi
#endif


      implicit none

       real(r8), parameter :: zero=0.0_r8,  one=1.0_r8,  two=2.0_r8     &
     &,                       three=3.0_r8, four=4.0_r8, five=5.0_r8    &
     &,                       half=0.5_r8,  oneb3=one/three             &
     &,                       onebcp=one/cpair
#ifdef NEMS_GSM
!
       integer,  parameter :: iulog = 6

       real(r8), parameter :: rhmini = 0.80_r8
       real(r8), parameter :: rhmaxi = 1.1_r8
!      real(r8), parameter :: r_universal = 6.02214e26*1.38065e-23
!      real(r8), parameter :: r_universal = con_rgas * 1000.0
       real(r8), parameter :: mwh2o   = 18.016
       real(r8), parameter :: rhoh2o  = 1.000e3

       logical  :: ip = .true.
       real(r8) :: tmn = 173.16_r8, tmx = 375.16_r8, trice = 35.00_r8
#endif


#ifdef GEOS5
!++jtb
!
       real(r8), parameter :: rhmini = 0.80_r8
       real(r8), parameter :: rhmaxi = 1.1_r8
       integer,  parameter :: iulog = 6



       real(r8), parameter :: gravit = MAPL_GRAV
       real(r8), parameter :: rair = MAPL_RGAS
       real(r8), parameter :: tmelt = MAPL_TICE
       real(r8), parameter :: cpair = MAPL_CP
       real(r8), parameter :: rh2o = MAPL_RVAP
       real(r8), parameter :: r_universal = MAPL_RUNIV
       real(r8), parameter :: mwh2o = MAPL_H2OMW
       real(r8), parameter :: rhoh2o = MAPL_RHOWTR
       real(r8), parameter :: latvap = MAPL_ALHL
       real(r8), parameter :: latice = MAPL_ALHF
       real(r8), parameter :: epsilon = MAPL_VIREPS


       logical  :: ip    = .true.
       real(r8) :: tmn   = 173.16_r8
       real(r8) :: tmx   = 375.16_r8
       real(r8) :: trice = 35.00_r8
#endif
!--jtb


       private
!      save

       logical, public :: liu_in = .false.


       public :: ini_micro, mmicro_pcond,gamma,derf

!constants remaped
       real(r8), private:: g, ginv
       real(r8), private:: r
       real(r8), private:: rv
!      real(r8), private:: rr             ! not used
       real(r8), private:: cpp
       real(r8), private:: rhow, pirhow
       real(r8), private:: xxlv
       real(r8), private:: xlf, xlfocp, cpoxlf
       real(r8), private:: xxls

      real(r8), private:: rhosn, pirhosn
      real(r8), private:: rhoi,  pirhoi

      real(r8), private:: ac,bc,as,bs,ai,bi,ar,br
      real(r8), private:: ci,di,oneodi
      real(r8), private:: cs,ds
      real(r8), private:: cr,dr
      real(r8), private:: f1s,f2s
      real(r8), private:: Eii
      real(r8), private:: Ecc
      real(r8), private:: Ecr
      real(r8), private:: f1r,f2r
      real(r8), private:: DCS, ts_auto_ice
      real(r8), private:: qsmall
      real(r8), private:: qvsmall
      real(r8), private:: bimm,aimm
      real(r8), private:: rhosu
      real(r8), private:: mi0
      real(r8), private:: rin
      real(r8), private:: qcvar
      real(r8), private:: pi

! Additional constants to help speed up code

      real(r8), private:: cons1,  cons2,  cons3,  cons4,  cons5
     &,                   cons6,  cons7,  cons8,  cons9,  cons10
     &,                   cons11, cons12, cons13, cons14, cons15
     &,                   cons16, cons17, cons18, cons19, cons20
     &,                   cons21, cons22, cons23, cons24, cons25
     &,                   cons27, cons28

      real(r8), private:: lammini, lammaxi, lamminr, lammaxr
     &,                   lammins, lammaxs

! parameters for snow/rain fraction for convective clouds
      real(r8), private, parameter :: tmax_fsnow = tmelt
     &,                               tmin_fsnow = tmelt-5._r8

!needed for findsp
      real(r8), private:: tt0

!switch for specification of droplet and crystal number

      real(r8), private:: csmin,csmax,minrefl,mindbz


       contains

!===============================================================================

      subroutine ini_micro(Dcs_, QCVAR_, ts_auto_ice_)

!-----------------------------------------------------------------------
!
! Purpose:
! initialize constants for the morrison microphysics
! called from stratiform.F90
!
! Author: Andrew Gettelman Dec 2005
!
!-----------------------------------------------------------------------

#ifdef CAM
       use cloud_fraction, only: cldfrc_getparams
#endif
       real(r8), intent(in) :: Dcs_, QCVAR_, ts_auto_ice_


       integer k, l, m, iaer
       real(r8) surften, arg, derf

       character(len=16) :: eddy_scheme = ' '
       logical           :: history_microphysics



#ifdef CAM 

       call phys_getopts(eddy_scheme_out = eddy_scheme,
     & history_microphysics_out = history_microphysics )


       call addfld ('QRAIN ','kg/kg ',pver, 'A','Diagnostic grid-
     &mean rain mixing ratio' ,phys_decomp)
       call addfld ('QSNOW ','kg/kg ',pver, 'A','Diagnostic grid-
     &mean snow mixing ratio' ,phys_decomp)
       call addfld ('NRAIN ','m-3 ',pver, 'A','Diagnostic grid-
     &mean rain number conc' ,phys_decomp)
       call addfld ('NSNOW ','m-3 ',pver, 'A','Diagnostic grid-
     &mean snow number conc' ,phys_decomp)


       call addfld ('RERCLD ','m ',pver, 'A',
     &'Diagnostic effective radius of Liquid Cloud and Rain' ,
     &phys_decomp)

       call addfld ('DSNOW ','m ',pver, 'A','Diagnostic grid-
     &mean snow diameter' ,phys_decomp)


       call addfld ('MGFLXPRC ','kg/m2/s ',pver+1, 'A','Diagnostic grid-
     &mean rain flux at layer interface', phys_decomp)
       call addfld ('MGFLXSNW ','kg/m2/s ',pver+1, 'A','Diagnostic grid-
     &mean snow flux at layer interface', phys_decomp)


       call addfld ('REFL ','DBz ',pver, 'A',
     &'94 GHz radar reflectivity' ,phys_decomp)
       call addfld ('AREFL ','DBz ',pver, 'A',
     &'Average 94 GHz radar reflectivity' ,phys_decomp)
       call addfld ('FREFL ','fraction ',pver, 'A',
     &'Fractional occurance of radar reflectivity' ,phys_decomp)

       call addfld ('CSRFL ','DBz ',pver, 'A',
     &'94 GHz radar reflectivity (CloudSat thresholds)' ,phys_decomp)
       call addfld ('ACSRFL ','DBz ',pver, 'A',
     &'Average 94 GHz radar reflectivity (CloudSat thresholds)' ,
     &phys_decomp)
       call addfld ('FCSRFL ','fraction ',pver, 'A',
     &'Fractional occurance of radar reflectivity (CloudSat thresholds)'
     & ,phys_decomp)

       call addfld ('AREFLZ ','mm^6/m^3 ',pver, 'A',
     &'Average 94 GHz radar reflectivity' ,phys_decomp)


       call addfld ('NCAL ','#/m3 ',pver, 'A',
     &'Number Concentation Activated for Liquid',phys_decomp)
       call addfld ('NCAI ','#/m3 ',pver, 'A',
     &'Number Concentation Activated for Ice',phys_decomp)



       call addfld ('AQRAIN ','kg/kg ',pver, 'A',
     &'Average rain mixing ratio' ,phys_decomp)
       call addfld ('AQSNOW ','kg/kg ',pver, 'A',
     &'Average snow mixing ratio' ,phys_decomp)
       call addfld ('ANRAIN ','m-3 ',pver, 'A',
     &'Average rain number conc' ,phys_decomp)
       call addfld ('ANSNOW ','m-3 ',pver, 'A',
     &'Average snow number conc' ,phys_decomp)
       call addfld ('ADRAIN ','Micron ',pver, 'A',
     &'Average rain effective Diameter' ,phys_decomp)
       call addfld ('ADSNOW ','Micron ',pver, 'A',
     &'Average snow effective Diameter' ,phys_decomp)
       call addfld ('FREQR ','fraction ',pver, 'A',
     &'Fractional occurance of rain' ,phys_decomp)
       call addfld ('FREQS ','fraction ',pver, 'A',
     &'Fractional occurance of snow' ,phys_decomp)

       if ( history_microphysics) then
       call add_default ('AQSNOW ', 1, ' ')
       call add_default ('FREQR ', 1, ' ')
       call add_default ('FREQS ', 1, ' ')
       call add_default ('AQRAIN ', 1, ' ')
       call add_default ('AQSNOW ', 1, ' ')
       call add_default ('ANRAIN ', 1, ' ')
       call add_default ('ANSNOW ', 1, ' ')
       end if

#endif
!--jtb commented out when GEOS5

!declarations for morrison codes (transforms variable names)

       g    = gravit
       ginv = one / g
       r    = rair
       rv   = rh2o
!      rr   = r_universal
       cpp  = cpair
       cp  = cpair
       rhow = rhoh2o

! latent heats

       xxlv   = latvap
       xlf    = latice
       xxls   = xxlv + xlf
       xlfocp = xlf  / cpair
       cpoxlf = cpair / xlf

!      write(0,*)' xlfocp=',xlfocp,' cpoxlf=',cpoxlf
! parameters below from Reisner et al. (1998)
! density parameters (kg/m3)

       rhosn = 100._r8
       rhoi  = 500._r8
       rhow  = 1000._r8


! fall speed parameters, V = aD^b
! V is in m/s

! droplets
       ac = 3.e7_r8
       bc = two

! snow
       as = 11.72_r8
       bs = 0.41_r8

! cloud ice
       ai = 700._r8
       bi = one

! rain
       ar = 841.99667_r8
       br = 0.8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d

!      pi= 3.1415927_r8
!      pi= 3.1415926535897931_r8

       pi = four*atan(one)

       pirhow  = pi * rhow
       pirhosn = pi * rhosn
       pirhoi  = pi * rhoi

! cloud ice mass-diameter relationship

       ci     = pirhoi/6._r8
       di     = three
       oneodi = one / di

! snow mass-diameter relationship

       cs = pirhosn/6._r8
       ds = three

! drop mass-diameter relationship

       cr = pirhow/6._r8
       dr = three

! ventilation parameters for snow
! hall and prupacher

       f1s = 0.86_r8
       f2s = 0.28_r8

! collection efficiency, aggregation of cloud ice and snow

       Eii = 0.1_r8

! collection efficiency, accretion of cloud water by rain

       Ecr = one

! ventilation constants for rain

       f1r = 0.78_r8
       f2r = 0.32_r8

! autoconversion size threshold for cloud ice to snow (m)

       Dcs = Dcs_ * 1.0e-6_r8
! ice autoconversion time scale (default 180s)

       ts_auto_ice = ts_auto_ice_

! smallest mixing ratio considered in microphysics

       qsmall = 1.e-18_r8
       qvsmall = 1.e-6_r8

! immersion freezing parameters, bigg 1953

       bimm = 100._r8
       aimm = 0.66_r8

! typical air density at 850 mb

       rhosu = 85000._r8/(rair * tmelt)

! mass of new crystal due to aerosol freezing and growth (kg)

       mi0 = (four/three)*pirhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)


! radius of contact nuclei aerosol (m)

       rin = 0.1e-6_r8

! 1 / relative variance of sub-grid cloud water distribution
! see morrison and gettelman, 2008, J. Climate for details



       qcvar = QCVAR_

! freezing temperature
       tt0  = 273.15_r8


!++jtb
       tmn   = 173.16_r8
       tmx   = 375.16_r8
       trice = 35.00_r8
       ip    = .true.

       call gestbl(tmn ,tmx ,trice ,ip ,epsilon , latvap ,latice ,rh2o ,
     &             cpair ,tmelt )

!--jtb (08/18/10)

!      pi = 4._r8*atan(1.0_r8)


       csmin   = -30._r8
       csmax   =  26._r8
       mindbz  = -99._r8

       minrefl = 1.26e-10_r8

! Define constants to help speed up code (limit calls to gamma function)

       cons1   = gamma(one+di)
!      cons1   = gamma(one+miu_ice+ di)/gamma(one+miu_ice)

       cons2   = gamma(qcvar+2.47_r8)
       cons3   = gamma(qcvar)
       cons4   = gamma(one+br)
       cons5   = gamma(four+br)
       cons6   = gamma(one+ds)
       cons7   = gamma(one+bs)
       cons8   = gamma(four+bs)
       cons9   = gamma(qcvar+two)
       cons10  = gamma(qcvar+one)
       cons11  = gamma(three+bs)
       cons12  = gamma(qcvar+1.15_r8)
       cons13  = gamma((five+br)*half)
       cons14  = gamma((five+bs)*half)
       cons15  = gamma(qcvar+bc/three)





       cons18  = qcvar**2.47_r8
       cons19  = qcvar*qcvar
       cons20  = qcvar**1.15_r8
       cons21  = qcvar**(bc/three)
       cons22  = (four/three)*pirhow*(25.e-6_r8)**3
       cons23  = dcs*dcs*dcs
       cons24  = dcs*dcs
       cons25  = dcs**bs
       cons27  = xxlv*xxlv
       cons28  = xxls*xxls

       lammaxi = one / 1.e-6_r8
       lammini = one / (one*dcs)
       lammaxr = one / 20.e-6_r8
       lamminr = one / 500.e-6_r8
       lammaxs = one / 10.e-6_r8
       lammins = one / 2000.e-6_r8

       return
       end subroutine ini_micro

!===============================================================================
!microphysics routine for each timestep goes here...

      subroutine mmicro_pcond ( lchnk, ncol, deltatin, tn, ttend,
     & pcols, pver,
     & qn, qtend, cwtend, qc, qi, nc, ni,fprcp,qrn,qsnw,nrn,nsnw,
     & p, pdel, cldn, liqcldf,
     & icecldf, cldo, pint, rpdel, zm,        rate1ord_cw2pr_st, naai,
!    & icecldf, cldo, pint, rpdel, zm, omega, rate1ord_cw2pr_st, naai,
     & npccnin, rndst,nacon, rhdfda, rhu00, fice, tlat, qvlat, qctend,
     & qitend, nctend, nitend, effc, effc_fn, effi, prect, preci,
     & nevapr, evapsnow, prain, prodsnow, cmeout, deffi, pgamrad,
     & lamcrad,qsout2,qrout2, drout2, qcsevap,qisevap,qvres,
     & cmeiout, vtrmc,vtrmi,qcsedten,qisedten, prao,
     & prco,mnuccco,mnuccto,
     & msacwio,psacwso, bergso,bergo,melto,homoo,qcreso,prcio,praio,
     & qireso, mnuccro,pracso,meltsdt,frzrdt, ncal, ncai, mnuccdo,
     & nnuccto,
     & nsout2, nrout2, ncnst, ninst, nimm, miu_disp, nsoot, rnsoot,
     & ui_scale, dcrit, nnuccdo, nnuccco, nsacwio, nsubio, nprcio,
     & npraio, npccno, npsacwso, nsubco, nprao, nprc1o, tlataux,
     & nbincontactdust, lprint,xlat,xlon)
!    & nbincontactdust, ts_auto_ice,xlat,xlon)


!Author: Hugh Morrison, Andrew Gettelman, NCAR
! e-mail: morrison@ucar.edu, andrew@ucar.edu

       use wv_saturation, only: vqsatd, vqsatd_water

#ifdef CAM
       use constituents, only: pcnst
       real(r8), parameter :: ncnst = 100.0e6
       real(r8), parameter :: ninst = 0.10e6
#endif

!--jtb (08/18/10) Is this still needed? Pcnst is not used.

       integer,  intent (in) :: pcols, pver,fprcp
       real(r8), intent (in) :: ncnst
       real(r8), intent (in) :: ninst
       real(r8), intent (in) :: nimm (pcols,pver)
       real(r8), intent (in) :: miu_disp , ui_scale, dcrit
!      real(r8), intent (in) :: miu_disp , ui_scale, dcrit, ts_auto_ice
       real(r8), intent (in) :: nsoot (pcols,pver) , rnsoot (pcols,pver)
     &     ,xlon,xlat
!      integer,  intent(in)  :: ktrop_min


       logical lprint


       integer, intent(in)  :: lchnk
       integer, intent(in)  :: ncol
       real(r8), intent(in) :: deltatin
       real(r8), intent(in) :: tn(pcols,pver)
       real(r8), intent(in) :: ttend(pcols,pver)
       real(r8), intent(in) :: qn(pcols,pver)
       real(r8), intent(in) :: qtend(pcols,pver)
       real(r8), intent(in) :: cwtend(pcols,pver)

       real(r8), intent(inout) :: qc(pcols,pver)
       real(r8), intent(inout) :: qi(pcols,pver)
       real(r8), intent(inout) :: nc(pcols,pver)
       real(r8), intent(inout) :: ni(pcols,pver)
       real(r8), intent(inout) :: qrn(pcols,pver)
       real(r8), intent(inout) :: qsnw(pcols,pver)
       real(r8), intent(inout) :: nrn(pcols,pver)
       real(r8), intent(inout) :: nsnw(pcols,pver)
       real(r8), intent(in)    :: p(pcols,pver)
       real(r8), intent(in)    :: pdel(pcols,pver)
       real(r8), intent(in)    :: cldn(pcols,pver)
       real(r8), intent(in)    :: icecldf(pcols,pver)
       real(r8), intent(in)    :: liqcldf(pcols,pver)
       real(r8), intent(inout) :: cldo(pcols,pver)
       real(r8), intent(in)    :: pint(pcols,pver+1)

       real(r8), intent(in)    :: rpdel(pcols,pver)
       real(r8), intent(in)    :: zm(pcols,pver)
!      real(r8), intent(in)    :: omega(pcols,pver)

       real(r8), intent(out)   :: rate1ord_cw2pr_st(pcols,pver)

! Inputs for aerosol activation
       real(r8), intent(in)    :: naai(pcols,pver)
       real(r8), intent(inout) :: npccnin(pcols,pver)
       integer                 :: nbincontactdust
       real(r8), intent(in), dimension(pcols,pver, 10) :: rndst, nacon


       real(r8), intent(in) :: rhdfda(pcols,pver)
       real(r8), intent(in) :: rhu00(pcols,pver)
       real(r8), intent(in) :: fice(pcols,pver)

       real(r8), intent(out) :: tlat(pcols,pver)

       real(r8), intent(out) :: tlataux(pcols,pver)

       real(r8), intent(out) :: qvlat(pcols,pver)
       real(r8), intent(out) :: qctend(pcols,pver)
       real(r8), intent(out) :: qitend(pcols,pver)
       real(r8), intent(out) :: nctend(pcols,pver)
       real(r8), intent(out) :: nitend(pcols,pver)
       real(r8), intent(out) :: effc(pcols,pver)
       real(r8), intent(out) :: effc_fn(pcols,pver)
       real(r8), intent(out) :: effi(pcols,pver)
       real(r8), intent(out) :: prect(pcols)
       real(r8), intent(out) :: preci(pcols)
       real(r8), intent(out) :: nevapr(pcols,pver)
       real(r8), intent(out) :: evapsnow(pcols,pver)
       real(r8), intent(out) :: prain(pcols,pver)
       real(r8), intent(out) :: prodsnow(pcols,pver)
       real(r8), intent(out) :: cmeout(pcols,pver)
       real(r8), intent(out) :: deffi(pcols,pver)
       real(r8), intent(out) :: pgamrad(pcols,pver)
       real(r8), intent(out) :: lamcrad(pcols,pver)

       real(r8), intent(out) :: qcsevap(pcols,pver)
       real(r8), intent(out) :: qisevap(pcols,pver)
       real(r8), intent(out) :: qvres(pcols,pver)
       real(r8), intent(out) :: cmeiout(pcols,pver)
       real(r8), intent(out) :: vtrmc(pcols,pver)
       real(r8), intent(out) :: vtrmi(pcols,pver)
       real(r8), intent(out) :: qcsedten(pcols,pver)
       real(r8), intent(out) :: qisedten(pcols,pver)
! microphysical process rates for output (mixing ratio tendencies)
       real(r8), intent(out) :: prao(pcols,pver)
       real(r8), intent(out) :: prco(pcols,pver)
       real(r8), intent(out) :: mnuccco(pcols,pver)
       real(r8), intent(out) :: mnuccto(pcols,pver)
       real(r8), intent(out) :: msacwio(pcols,pver)
       real(r8), intent(out) :: psacwso(pcols,pver)
       real(r8), intent(out) :: bergso(pcols,pver)
       real(r8), intent(out) :: bergo(pcols,pver)
       real(r8), intent(out) :: melto(pcols,pver)
       real(r8), intent(out) :: homoo(pcols,pver)
       real(r8), intent(out) :: qcreso(pcols,pver)
       real(r8), intent(out) :: prcio(pcols,pver)
       real(r8), intent(out) :: praio(pcols,pver)
       real(r8), intent(out) :: qireso(pcols,pver)
       real(r8), intent(out) :: mnuccro(pcols,pver)
       real(r8), intent(out) :: pracso (pcols,pver)
       real(r8), intent(out) :: meltsdt(pcols,pver)
       real(r8), intent(out) :: frzrdt (pcols,pver)
       real(r8), intent(out) :: mnuccdo(pcols,pver)
       real(r8), intent(out) :: nnuccto(pcols,pver)

       real(r8), intent(out) :: nnuccdo(pcols,pver)
       real(r8), intent(out) :: nnuccco(pcols,pver)
       real(r8), intent(out) :: nsacwio(pcols,pver)
       real(r8), intent(out) :: nsubio(pcols,pver)
       real(r8), intent(out) :: nprcio(pcols,pver)
       real(r8), intent(out) :: npraio(pcols,pver)

       real(r8), intent(out) :: npccno(pcols,pver)
       real(r8), intent(out) :: npsacwso(pcols,pver)
       real(r8), intent(out) :: nsubco(pcols,pver)
       real(r8), intent(out) :: nprao(pcols,pver)
       real(r8), intent(out) :: nprc1o(pcols,pver)

! local workspace
! all units mks unless otherwise stated

! temporary variables for sub-stepping

       real(r8) :: t1(pcols,pver)
       real(r8) :: q1(pcols,pver)
       real(r8) :: qc1(pcols,pver)
       real(r8) :: qi1(pcols,pver)
       real(r8) :: nc1(pcols,pver)
       real(r8) :: ni1(pcols,pver)
       real(r8) :: tlat1(pcols,pver)
       real(r8) :: qvlat1(pcols,pver)
       real(r8) :: qctend1(pcols,pver)
       real(r8) :: qitend1(pcols,pver)
       real(r8) :: nctend1(pcols,pver)
       real(r8) :: nitend1(pcols,pver)
       real(r8) :: prect1(pcols)
       real(r8) :: preci1(pcols)

       real(r8) :: tlat1_aux(pcols,pver)


! hm 5/12/11
! temporary variable for old nc before updating with activation tendency
       real(r8) :: ncold(pcols,pver)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       real(r8) :: deltat, dti, deltam
       real(r8) :: omsm
       real(r8) :: dto2
       real(r8) :: mincld
       real(r8), dimension(pcols,pver) :: q, t, rho, irho, rhof, dv
     &,                                   mu, sc, kap, cldmax, cldm
     &,                                   icldm, lcldm, cme, cmei
     &,                                   cwml, cwmi, lcldn, lcldo
     &,                                   nctend_mixnuc, npccn

       real(r8), dimension(pcols)      :: icwc, calpha, cbeta, cbetah
     &,                                   cgamma, cgamah, rcgama
     &,                                   cmec1, cmec2, cmec3, cmec4

       real(r8), dimension(pver)       :: nnuccd, mnuccd
     &,                                   qcsinksum_rate1ord
     &,                                   qcsum_rate1ord

       real(r8)                        :: qtmp, dum, dum1, dum2, qcld
     &,                                   arg, alpha


       real(r8) :: qcic(pcols,pver)
       real(r8) :: qiic(pcols,pver)
       real(r8) :: qniic(pcols,pver)
       real(r8) :: qric(pcols,pver)
       real(r8) :: ncic(pcols,pver)
       real(r8) :: niic(pcols,pver)
       real(r8) :: nsic(pcols,pver)
       real(r8) :: nric(pcols,pver)
       real(r8) :: lami(pver)
       real(r8) :: n0i(pver)
       real(r8) :: lamc(pver)
       real(r8) :: n0c(pver)
       real(r8) :: lams(pver)
       real(r8) :: n0s(pver)
       real(r8) :: lamr(pver)
       real(r8) :: n0r(pver)
       real(r8) :: cdist1(pver)
! combined size of precip & cloud drops
       real(r8) :: rercld(pcols,pver)
       real(r8) :: arcld(pcols,pver)
       real(r8) :: Actmp
       real(r8) :: Artmp

       real(r8) :: pgam(pver)
       real(r8) :: lammax
       real(r8) :: lammin
       real(r8) :: nacnt
       real(r8) :: mnuccc(pver)
       real(r8) :: nnuccc(pver)

       real(r8) :: mnucct(pver)
       real(r8) :: nnucct(pver)
       real(r8) :: msacwi(pver)
       real(r8) :: nsacwi(pver)

       real(r8) :: prc(pver)
       real(r8) :: nprc(pver)
       real(r8) :: nprc1(pver)
       real(r8) :: nsagg(pver)
       real(r8) :: dc0
       real(r8) :: ds0
       real(r8) :: eci
       real(r8) :: psacws(pver)
       real(r8) :: npsacws(pver)
       real(r8) :: uni
       real(r8) :: umi
       real(r8) :: uns(pver)
       real(r8) :: ums(pver)
       real(r8) :: unr(pver)
       real(r8) :: umr(pver)
       real(r8) :: unc
       real(r8) :: umc
       real(r8) :: pracs(pver)
       real(r8) :: npracs(pver)
       real(r8) :: mnuccr(pver)
       real(r8) :: nnuccr(pver)
       real(r8) :: pra(pver)
       real(r8) :: npra(pver)
       real(r8) :: nragg(pver)
       real(r8) :: prci(pver)
       real(r8) :: nprci(pver)
       real(r8) :: prai(pver)
       real(r8) :: nprai(pver)
       real(r8) :: qvs
       real(r8) :: qvi
       real(r8) :: dqsdt
       real(r8) :: dqsidt
       real(r8) :: ab
       real(r8) :: qclr
       real(r8) :: abi,oneoabi
       real(r8) :: epss
       real(r8) :: epsr
       real(r8) :: pre(pver)
       real(r8) :: prds(pver)
       real(r8) :: qce
       real(r8) :: qie
       real(r8) :: nce
       real(r8) :: nie
       real(r8) :: ratio
       real(r8) :: dumc(pcols,pver)
       real(r8) :: dumnc(pcols,pver)
       real(r8) :: dumi(pcols,pver)
       real(r8) :: dumni(pcols,pver)
       real(r8) :: dums(pcols,pver)
       real(r8) :: dumns(pcols,pver)
       real(r8) :: dumr(pcols,pver)
       real(r8) :: dumnr(pcols,pver)
! below are parameters for cloud water and cloud ice sedimentation calculations
       real(r8) :: fr(pver)
       real(r8) :: fnr(pver)
       real(r8) :: fc(pver)
       real(r8) :: fnc(pver)
       real(r8) :: fi(pver)
       real(r8) :: fni(pver)
       real(r8) :: fs(pver)
       real(r8) :: fns(pver)
       real(r8) :: faloutr(pver)
       real(r8) :: faloutnr(pver)
       real(r8) :: faloutc(pver)
       real(r8) :: faloutnc(pver)
       real(r8) :: falouti(pver)
       real(r8) :: faloutni(pver)
       real(r8) :: falouts(pver)
       real(r8) :: faloutns(pver)
       real(r8) :: faltndr
       real(r8) :: faltndnr
       real(r8) :: faltndc
       real(r8) :: faltndnc
       real(r8) :: faltndi
       real(r8) :: faltndni
       real(r8) :: faltnds
       real(r8) :: faltndns
       real(r8) :: faltndqie
       real(r8) :: faltndqce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       real(r8) :: relhum(pcols,pver)
       real(r8) :: csigma(pcols)
       real(r8) :: rgvm
       real(r8) :: arn(pcols,pver)
       real(r8) :: asn(pcols,pver)
       real(r8) :: acn(pcols,pver)
       real(r8) :: ain(pcols,pver)
       real(r8) :: nsubi(pver)
       real(r8) :: nsubc(pver)
       real(r8) :: nsubs(pver)
       real(r8) :: nsubr(pver)
       real(r8) :: mtime
       real(r8) :: dz(pcols,pver)

!fice variable
       real(r8) :: nfice(pcols,pver)

!add variables for rain and snow flux at layer interfaces
       real(r8) :: rflx(pcols,pver+1)
       real(r8) :: sflx(pcols,pver+1)

       real(r8) :: rflx1(pcols,pver+1)
       real(r8) :: sflx1(pcols,pver+1)

! returns from function/subroutine calls
       real(r8) :: tsp(pcols,pver)
       real(r8) :: qsp(pcols,pver)
       real(r8) :: qsphy(pcols,pver)
       real(r8) :: qs(pcols)
       real(r8) :: es(pcols)
       real(r8) :: esl(pcols,pver)
       real(r8) :: esi(pcols,pver)
!      real(r8) :: gammas(pcols)

! sum of source/sink terms for diagnostic precip

       real(r8) :: qnitend(pcols,pver)
       real(r8) :: nstend(pcols,pver)
       real(r8) :: qrtend(pcols,pver)
       real(r8) :: nrtend(pcols,pver)
       real(r8) :: qrtot
       real(r8) :: nrtot
       real(r8) :: qstot
       real(r8) :: nstot

! new terms for Bergeron process

       real(r8) :: dumnnuc
       real(r8) :: ninew
       real(r8) :: qinew
       real(r8) :: qvl
       real(r8) :: epsi
       real(r8) :: prd
       real(r8) :: berg(pcols,pver)
       real(r8) :: bergs(pver)

!bergeron terms
       real(r8) :: bergtsf
       real(r8) :: rhin

! diagnostic rain/snow for output to history
! values are in-precip (local) !!!!

       real(r8) :: qrout(pcols,pver)
       real(r8) :: nrout(pcols,pver)
       real(r8) :: nsout(pcols,pver)
       real(r8) :: dsout(pcols,pver)
       real(r8) :: qsout(pcols,pver)


!averageed rain/snow for history
       real(r8) , intent(out) :: qrout2(pcols,pver)
       real(r8) , intent(out) :: qsout2(pcols,pver)
       real(r8) , intent(out) :: nrout2(pcols,pver)
       real(r8) , intent(out) :: nsout2(pcols,pver)
       real(r8) :: freqs(pcols,pver)
       real(r8) :: freqr(pcols,pver)
       real(r8) :: dumfice
       real(r8), intent(out) :: drout2(pcols,pver)
       real(r8)              :: dsout2(pcols,pver)

!ice nucleation, droplet activation
       real(r8) :: dum2i(pcols,pver)
       real(r8) :: dum2l(pcols,pver)
       real(r8) :: ncmax(pcols,pver)
       real(r8) :: nimax

!output fields for number conc
       real(r8) :: ncai(pcols,pver)
       real(r8) :: ncal(pcols,pver)

! loop array variables
       integer i,k,nstep,n, l
       integer ii,kk, m, ind_aux, km, kp

! loop variables for sub-step solution
       integer iter,it,ltrue(pcols)

! used in contact freezing via dust particles
       real(r8) tcnt, viscosity, mfp, nslipsoot, ndfaersoot
       real(r8), dimension(nbincontactdust) :: ndfaer, nslip, slip


! used in ice effective radius
       real(r8) bbi, cci, ak, iciwc, rvi, riter

! used in Bergeron processe and water vapor deposition
       real(r8) Tk, deles, Aprpr, Bprpr, Cice, qi0, Crate, qidep

! mean cloud fraction over the time step
       real(r8) cldmw(pcols,pver)

! used in secondary ice production
       real(r8) ni_secp

! variabels to check for RH after rain evap

       real(r8) :: esn
       real(r8) :: qsn
       real(r8) :: ttmp


       real(r8) :: refl(pcols,pver)
       real(r8) :: rainrt(pcols,pver)
       real(r8) :: rainrt1(pcols,pver)
       real(r8) :: csrfl(pcols,pver)
       real(r8) :: arefl(pcols,pver)
       real(r8) :: acsrfl(pcols,pver)
       real(r8) :: frefl(pcols,pver)
       real(r8) :: fcsrfl(pcols,pver)
       real(r8) :: areflz(pcols,pver)
       real(r8) :: tmp, miu_ice(pver)
       real(r8) :: cons16
       real(r8) :: cons17

       real(r8) dmc,ssmc,dstrn

       real(r8), parameter :: cdnl = 0.e6_r8

!      integer,  parameter :: auto_option=3     ! dcrit only used with auto_option=4
!      integer,  parameter :: auto_option=2     ! dcrit only used with auto_option=4
       integer,  parameter :: auto_option=1     ! dcrit only used with auto_option=4
!      integer,  parameter :: auto_option=4     ! dcrit only used with auto_option=4
       real(r8) :: beta6, xs, nssoot, nsdust, taux, psc, Bh, vaux, aux,
     &             LW, NW, tx1, tx2, tx3, tx4, tx5, omeps, esloesi
     &,            rdz, rdzi


! note: number will be adjusted as needed to keep mean size within bounds,
! even when cosntant droplet or ice number is used

! ***note: Even if constant cloud ice number is set, ice number is allowed
! to evolve based on process rates. This is needed in order to calculate
! the change in mass due to ice nucleation. All other ice microphysical
! processes are consistent with the specified constant ice number if
! this switch is turned on.

       logical :: nccons,nicons


       omeps  = one - epsqs

       nccons = .false.
       nicons = .false.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! initialize  output fields for number conc qand ice nucleation
       do k=1,pver
         miu_ice(k) = zero
       enddo
       do k=1,pver
         do i=1,ncol
           ncai(i,k)      = zero
           ncal(i,k)      = zero

!Initialize rain size
           rercld(i,k)    = zero
           arcld(i,k)     = zero

!initialize radiation output variables
           pgamrad(i,k)   = zero
           lamcrad(i,k)   = zero
           deffi (i,k)    = zero
!initialize radiation output variables
!initialize water vapor tendency term output
           qcsevap(i,k)   = zero
           qisevap(i,k)   = zero
           qvres (i,k)    = zero
           cmeiout (i,k)  = zero
           vtrmc (i,k)    = zero
           vtrmi (i,k)    = zero
           qcsedten (i,k) = zero
           qisedten (i,k) = zero

           prao(i,k)      = zero
           prco(i,k)      = zero

           mnuccco(i,k)   = zero
           mnuccto(i,k)   = zero
           msacwio(i,k)   = zero
           psacwso(i,k)   = zero
           bergso(i,k)    = zero
           bergo(i,k)     = zero
           melto(i,k)     = zero
           homoo(i,k)     = zero
           qcreso(i,k)    = zero
           prcio(i,k)     = zero
           praio(i,k)     = zero
           qireso(i,k)    = zero
           mnuccro(i,k)   = zero
           pracso (i,k)   = zero
           meltsdt(i,k)   = zero
           frzrdt (i,k)   = zero
           mnuccdo(i,k)   = zero
           nnuccto(i,k)   = zero


           nnuccdo(i,k)   = zero
           nnuccco(i,k)   = zero
           nsacwio(i,k)   = zero
           nsubio(i,k)    = zero
           nprcio(i,k)    = zero
           npraio(i,k)    = zero

           npccno(i,k)    = zero
           npsacwso(i,k)  = zero
           nsubco(i,k)    = zero
           nprao(i,k)     = zero
           nprc1o(i,k)    = zero
         enddo
       enddo

! assign variable deltat for sub-stepping...
       deltat = deltatin
       dti    = one / deltat
       deltam = one / max(deltat, 150.0_r8)

! parameters for scheme

!      omsm   = 0.99999_r8
       omsm   = 0.99999999_r8
       dto2   = 0.5_r8*deltat
       mincld = 0.00001_r8

! initialize time-varying parameters

       do k=1,pver
         do i=1,ncol

! initialize multi-level fields
           t(i,k)    = tn(i,k)
           q(i,k)    = qn(i,k)
           rho(i,k)  = p(i,k) / (r*t(i,k))
           irho(i,k) = one / rho(i,k)
           dv(i,k)   = 8.794E-5_r8*t(i,k)**1.81_r8/p(i,k)
           tx1       = t(i,k) * sqrt(t(i,k)) / (t(i,k)+120._r8)
           mu(i,k)   = 1.496E-6_r8 * tx1
           sc(i,k)   = mu(i,k)/(rho(i,k)*dv(i,k))
           kap(i,k)  = 1.414e3_r8 * 1.496e-6_r8 * tx1

! air density adjustment for fallspeed parameters
! includes air density correction factor to the
! power of 0.54 following Heymsfield and Bansemer 2007

           rhof(i,k) = (rhosu*irho(i,k))**0.54

           arn(i,k)  = ar * rhof(i,k)
           asn(i,k)  = as * rhof(i,k)
           acn(i,k)  = ac * rhof(i,k)
           ain(i,k)  = ai * rhof(i,k)

! get dz from dp and hydrostatic approx
! keep dz positive (define as layer k-1 - layer k)

           dz(i,k)   = pdel(i,k) * irho(i,k) * ginv

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! droplet activation
! hm, modify 5/12/11
! get provisional droplet number after activation. This is used for
! all microphysical process calculations, for consistency with update of
! droplet mass before microphysics

! calculate potential for droplet activation if cloud water is present
! tendency from activation (npccnin) is read in from companion routine

! hm note: npccn is no longer needed below this code - so this can
! be rewwritten and this parameters can be removed

           lcldm(i,k) = max(liqcldf(i,k),mincld)
           if (qc(i,k) >= qsmall) then
             npccn(i,k) = max( (lcldm(i,k)*npccnin(i,k)-nc(i,k))
     &                          * deltam, zero)
             ncold(i,k) = nc(i,k)
             nc(i,k)    = nc(i,k) + npccn(i,k)*deltat
           else
             npccn(i,k) = zero
           end if

! initialization
           t1(i,k)  = t(i,k)
           q1(i,k)  = q(i,k)
           qc1(i,k) = qc(i,k)
           qi1(i,k) = qi(i,k)
           nc1(i,k) = nc(i,k)
           ni1(i,k) = ni(i,k)

! initialize tendencies to zero
           tlat1(i,k)     = zero
           tlat1_aux(i,k) = zero

           qvlat1(i,k)  = zero
           qctend1(i,k) = zero
           qitend1(i,k) = zero
           nctend1(i,k) = zero
           nitend1(i,k) = zero

! initialize precip output
           qrout(i,k) = zero
           qsout(i,k) = zero
           nrout(i,k) = zero
           nsout(i,k) = zero
           dsout(i,k) = zero

! initialize variables for trop_mozart
           nevapr(i,k)   = zero
           evapsnow(i,k) = zero
           prain(i,k)    = zero
           prodsnow(i,k) = zero
           cmeout(i,k)   = zero

! for refl calc
           rainrt1(i,k) = zero

! initialize precip fraction and output tendencies
           cldmax(i,k) = mincld

!initialize aerosol number
!          naer2(i,k,:) = zero
           dum2l(i,k) = zero
           dum2i(i,k) = zero 
           ncmax(i,k) = zero

! for debug purpose
!          prect(1:ncol)=0._r8
!          preci(1:ncol)=0._r8
!          tlat(1:ncol,1:pver)=0._r8
!          qvlat(1:ncol,1:pver)=0._r8
!          qctend(1:ncol,1:pver)=0._r8
!          qitend(1:ncol,1:pver)=0._r8
!          nctend(1:ncol,1:pver)=0._r8

         enddo
       enddo

! initialize avg precip rate
       do i=1,ncol
         prect1(i) = zero
         preci1(i) = zero
       end do
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Get humidity and saturation vapor pressures  big loop1

       do k=1,pver

! find wet bulk temperature and saturation value for provisional t and q without
! condensation

!        call vqsatd_water(t(1,k),p(1,k),es,qs,gammas,ncol)

         do i=1,ncol

#ifdef NEMS_GSM
           esl(i,k) = min(fpvsl(t(i,k)), p(i,k))
           esi(i,k) = min(fpvsi(t(i,k)), p(i,k))
#else
           esl(i,k) = min(polysvp(t(i,k),0), p(i,k))
           esi(i,k) = min(polysvp(t(i,k),1), p(i,k))
#endif

           esloesi  = esl(i,k) / esi(i,k)
! hm fix, make sure when above freezing that esi=esl, not active yet
           if (t(i,k) > tmelt) esi(i,k) = esl(i,k)

           qs(i)       = epsqs*esl(i,k)/(p(i,k)-omeps*esl(i,k))

           relhum(i,k) = min(q(i,k)/qs(i), one) !Anning limiting relhum

!     if (lprint .and. k==29) write(0,*)' esl=',esl(i,k)
!    &,' esi=',esi(i,k),' pres=',p(i,k),' t=',t(i,k)
!    &,' relhum=',relhum(i,k),' q=',q(i,k),' qs=',qs(i)

! get cloud fraction, check for minimum

           cldm(i,k)  = max(cldn(i,k), mincld)
           cldmw(i,k) = max(cldn(i,k), mincld)

           icldm(i,k) = max(icecldf(i,k), mincld)
           lcldm(i,k) = max(liqcldf(i,k), mincld)


           if (qc(i,k) >= qsmall) then
             tx1        = one / lcldm(i,k)
             dum2l(i,k) = (ncold(i,k)+npccn(i,k)*deltat) * tx1
             dum2l(i,k) = max(dum2l(i,k),cdnl*irho(i,k))
             ncmax(i,k) = max(dum2l(i,k)*lcldm(i,k), zero)
             dum2l(i,k) = npccn(i,k)*deltat*tx1
           else
             dum2l(i,k) = zero
           end if



! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)

           nfice(i,k) = zero
           dumfice    = qc(i,k) + qi(i,k)
           if (dumfice > qsmall .and. qi(i,k) > qsmall) then
             nfice(i,k) = qi(i,k)/dumfice
           endif

           if (t(i,k) < tmelt - five) then
             if (liu_in) then

! if aerosols interact with ice set number of activated ice nuclei
               dum2 = naai(i,k)

             else
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
               dum2 = 0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))*1000._r8
! put limit on number of nucleated crystals, set to number at T=-30 C
! cooper (limit to value at -35 C)
               dum2 = min(dum2, 208.9e3_r8)*irho(i,k)
             endif

             dumnnuc = max((dum2*icldm(i,k)-ni(i,k))*deltam, zero)


! get provisional ni and qi after nucleation in order to calculate
! Bergeron process below
             ninew = ni(i,k) + dumnnuc*deltat
             qinew = qi(i,k) + dumnnuc*deltat*mi0
!T>268
           else
             ninew = ni(i,k)
             qinew = qi(i,k)
           end if

! Initialize CME components

           cme(i,k)  = zero
           cmei(i,k) = zero


!-------------------------------------------------------------------
!Bergeron process
! make sure to initialize bergeron process to zero
           berg(i,k) = zero
           prd       = zero

!condensation loop.

! get in-cloud qi and ni after nucleation
           if (icldm(i,k) > zero) then
             tx1 = one / icldm(i,k)
             qiic(i,k) = qinew * tx1
             niic(i,k) = ninew * tx1
           else
             qiic(i,k) = zero
             niic(i,k) = zero
           endif


           if (nicons) then
             niic(i,k) = ninst*irho(i,k)
           end if



!if T < 0 C then bergeron.
           if (t(i,k) < 273.15) then
!if ice exists
             if (qi(i,k) > qsmall) then

               bergtsf = zero

               qvi = epsqs*esi(i,k) / (p(i,k)-omeps*esi(i,k))
               qvl = epsqs*esl(i,k) / (p(i,k)-omeps*esl(i,k))

               dqsidt  = xxls*qvi / (rv*t(i,k)*t(i,k))
               abi     = one + dqsidt*xxls/cpp
               oneoabi = one / abi

! get ice size distribution parameters

               if (qiic(i,k) >= qsmall) then

                 lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**oneodi

!      miu_ice=mui_hemp_l(lami(k))
                 miu_ice(k) = max(min(0.008_r8*(lami(k)*0.01)**0.87_r8,
     &                                      10.0_r8), 0.1_r8) 
                 tx1        = one + miu_ice(k)
                 tx2        = one / gamma(tx1)
                 aux        = (gamma(tx1+di)*tx2) ** oneodi
                 lami(k)    = aux*lami(k)
                 n0i(k)     = niic(i,k) * lami(k)**tx1 * tx2


! check for slope
! adjust vars
                 if (lami(k) < lammini*aux) then
                   miu_ice(k) = zero
                   lami(k)    = lammini
                   n0i(k)     = lami(k)**(di+one)*qiic(i,k)/(ci*cons1)
                   niic(i,k)  = n0i(k)/lami(k)
                 else if (lami(k) > lammaxi*aux) then
                   miu_ice(k) = zero
                   lami(k)    = lammaxi
                   n0i(k)     = lami(k)**(di+one)*qiic(i,k)/(ci*cons1)
                   niic(i,k)  = n0i(k)/lami(k)
                 end if

                 epsi = (miu_ice(k)+two)*(miu_ice(k)+one)*pi*
     &                  niic(i,k)*rho(i,k)*Dv(i,k)/lami(k)

!if liquid exists
                 if (qc(i,k) > qsmall) then

!begin bergeron process
!     do bergeron (vapor deposition with RHw=1)
!     code to find berg (a rate) goes here

! calculate Bergeron process

                   prd = epsi*(qvl-qvi)*oneoabi

                 else
                   prd = zero
                 end if

! multiply by cloud fraction and transfer of existing cloud liquid to ice

                 berg(i,k) = max(zero, prd*min(icldm(i,k),lcldm(i,k)))

               end if                  ! qiic(i,k) >= qsmall


               if (berg(i,k) > zero) then
                 tx1     = qc(i,k) * dti
                 bergtsf = max(zero, tx1/berg(i,k))
                 if(bergtsf < one) berg(i,k) = max(zero, tx1)
               endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               if (bergtsf < one .or. icldm(i,k) >  lcldm(i,k)) then

                 if (qiic(i,k) >= qsmall) then

! first case is for case when liquid water is present, but is completely depleted in time step, i.e., bergrsf > 0 but < 1

                   if (qc(i,k) >= qsmall) then
                     rhin = (one + relhum(i,k)) * half
                     if (rhin*esloesi > one) then
                       prd = epsi*(rhin*qvl-qvi)*oneoabi

! multiply by cloud fraction assuming liquid/ice maximum overlap and add to cmei
                       cmei(i,k) = cmei(i,k)
     &                           + prd * min(icldm(i,k),lcldm(i,k))
     &                           * (one- bergtsf)
!     if (lprint .and. k == 29) write(0,*)' cmei1=',cmei(i,k),
!    &' prd=',prd,' bergtsf=',bergtsf

                     end if
                   end if

!    moved here by Moorthi
! note: for case of no liquid, need to set liquid cloud fraction to zero
! store liquid cloud fraction in 'dum'
                   if (qc(i,k) < qsmall) then
                     dum = zero
                   else
                     dum = lcldm(i,k)
                   end if

! second case is for pure ice cloud, either no liquid, or icldm > lcldm
                   if (qc(i,k) < qsmall) then

! note: for case of no liquid, need to set liquid cloud fraction to zero
! store liquid cloud fraction in 'dum'

!Moorthi             if (qc(i,k) < qsmall) then
!                      dum = 0._r8
!                    else
!                      dum = lcldm(i,k)
!                    end if

! set RH to grid-mean value for pure ice cloud
                     rhin = relhum(i,k)

                     if (rhin*esloesi > one) then

                       prd = epsi*(rhin*qvl-qvi)*oneoabi

! multiply by relevant cloud fraction for pure ice cloud
! assuming maximum overlap of liquid/ice
                       cmei(i,k) = cmei(i,k)
     &                           + prd * max((icldm(i,k)-dum), zero)
!     if (lprint .and. k == 29) write(0,*)' cmei2=',cmei(i,k),
!    &' prd=',prd,' icldm=',icldm(i,k),' dum=',dum

                     end if
                   end if
                 end if
               end if

!     if deposition, it should not reduce grid mean rhi below 1.0
               if(cmei(i,k) > zero .and. relhum(i,k)*esloesi > one)
     &            cmei(i,k) = min(cmei(i,k),(q(i,k)-qs(i)/esloesi)
     &                                              * oneoabi * dti)

!     if (lprint .and. k == 29) write(0,*)' cmei3=',cmei(i,k)

             end if ! if (qi(i,k) > qsmall)


!-------------------------------------------------------------------
           end if
!..............................................................

! evaporation should not exceed available water

           tx1 = qc(i,k)*dti
           if (-berg(i,k) < -tx1) berg(i,k) = max(tx1, zero) ! is this correct?????

!          berg(i,k) = min(berg(i,k), qc(i,k)/deltat) ! Moorthi - ask anning

!sublimation process...

           if (relhum(i,k)*esloesi < one .and.
     &          qiic(i,k) >= qsmall ) then

             qvi = epsqs*esi(i,k)/(p(i,k)-omeps*esi(i,k))
             qvl = epsqs*esl(i,k)/(p(i,k)-omeps*esl(i,k))
             dqsidt  = xxls*qvi/(rv*t(i,k)*t(i,k))
             abi     = one + dqsidt*xxls/cpp
             oneoabi = one / abi

! get ice size distribution parameters

             lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**oneodi

!            miu_ice(k)=mui_hemp_l(lami(k))
             miu_ice(k) = max(min(0.008_r8*(lami(k)*0.01)**0.87_r8,
     &                                        10.0_r8), 0.1_r8)
             tx1        = one + miu_ice(k)
             tx2        = one / gamma(tx1)
             aux        = (gamma(tx1+di) * tx2) ** oneodi
             lami(k)    = aux*lami(k)

             n0i(k)     = niic(i,k) * lami(k)**tx1 * tx2
! check for slope
! adjust vars
             if (lami(k) < lammini*aux) then
               miu_ice(k) = zero
               lami(k) = lammini
               n0i(k)  = lami(k)**(di+one)*qiic(i,k)/(ci*cons1)
             else if (lami(k) > lammaxi*aux) then
               miu_ice(k) = zero
               lami(k) = lammaxi
               n0i(k)  = lami(k)**(di+one)*qiic(i,k)/(ci*cons1)
             end if

             epsi = (miu_ice(k)+two)*(miu_ice(k)+one)*pi*
     &               niic(i,k)*rho(i,k)*Dv(i,k) / lami(k)

! modify for ice fraction below
             prd = epsi*(relhum(i,k)*qvl-qvi)*oneoabi * icldm(i,k)
             cmei(i,k) = min(prd, zero)

!     if (lprint .and. k == 29) write(0,*)' cmei3a=',cmei(i,k)
!    &,' prd=',prd
           endif

! sublimation should not exceed available ice

           cmei(i,k) = max(cmei(i,k), -qi(i,k)*dti)

!     if (lprint .and. k == 29) write(0,*)' cmei3b=',cmei(i,k)
!    &,' qi=',qi(i,k),' deltat=',deltat

! sublimation should not increase grid mean rhi above 1.0
           if(cmei(i,k) < zero .and. relhum(i,k)*esloesi < one)
     &           cmei(i,k) = min(zero, max(cmei(i,k),(q(i,k)-qs(i)
     &                                /esloesi) * oneoabi * dti ))
!     if (lprint .and. k == 29) write(0,*)' cmei3c=',cmei(i,k)
!    &,' q=',q(i,k),' qs=',qs(i),' esi=',esi(i,k)
!    &,' esl=',esl(i,k),' abi=',abi

! limit cmei due for roundoff error

           cmei(i,k) = cmei(i,k)*omsm
!     if (lprint .and. k == 29) write(0,*)' cmei4=',cmei(i,k),
!    &' omsm=',omsm

! conditional for ice nucleation

           if (t(i,k) < (tmelt - five)) then
             if ( liu_in ) then

! using Liu et al. (2007)/Barahona & Nenes (2009) ice nucleation with hooks into simulated aerosol
! ice nucleation rate (dum2) has already been calculated and read in (naai)

               dum2i(i,k) = naai(i,k)
             else
   
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
               dum2i(i,k) = 0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k)))
     &                     * 1000._r8
! put limit on number of nucleated crystals, set to number at T=-30 C
! cooper (limit to value at -35 C)
               dum2i(i,k) = min(dum2i(i,k),208.9e3_r8)*irho(i,k)
             endif
           else
             dum2i(i,k) = zero
           end if

         end do ! end big i loop
       end do   !end big loop1 - k loop

!!
       do i=1,ncol
         rflx1(i,1) = zero     ! initialize sub-step precip flux variables
         sflx1(i,1) = zero

         rflx(i,1)  = zero     ! initialize final precip flux variables.
         sflx(i,1)  = zero
         ltrue(i)   = 0
       end do
       do k=1,pver
         do i=1,ncol
           cldo(i,k)    = cldn(i,k)

           rflx1(i,k+1) = zero ! initialize sub-step precip flux variables
           sflx1(i,k+1) = zero

           rflx(i,k+1)  = zero ! initialize final precip flux variables.
           sflx(i,k+1)  = zero

! skip microphysical calculations if no cloud water

           if ((qc(i,k) >= qsmall .or. qi(i,k) >= qsmall .or.
     &         cmei(i,k) >= qsmall).and.q(i,k)>=qvsmall) ltrue(i) = 1

           rate1ord_cw2pr_st(i,k) = zero
         end do
       end do
!!


! assign number of sub-steps to iter
! use 2 sub-steps, following tests described in MG2008
! Anning Cheng 9/17/2016
       if(fprcp==1) then
         iter = 1
       else
         iter  = 2
       end if

       riter = one / float(iter)
! get sub-step time step
       deltat = deltat * riter
       dti    = one / deltat

! since activation/nucleation processes are fast, need to take into account
! factor mtime = mixing timescale in cloud / model time step
! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk

!        note: mtime for bulk aerosols was set to: mtime=deltat/1200._r8

       mtime = one

!!!! skip calculations if no cloud water

       do i=1,ncol   !big i loop2
         if (ltrue(i) == 0) then
           prect(i) = zero
           preci(i) = zero
           do k=1,pver
             tlat(i,k)    = zero
             qvlat(i,k)   = zero
             qctend(i,k)  = zero
             qitend(i,k)  = zero
             qnitend(i,k) = zero
             qrtend(i,k)  = zero
             nctend(i,k)  = zero
             nitend(i,k)  = zero
             nrtend(i,k)  = zero
             nstend(i,k)  = zero
             qniic(i,k)   = zero
             qric(i,k)    = zero
             nsic(i,k)    = zero
             nric(i,k)    = zero
             rainrt(i,k)  = zero
           enddo

         else

           do k=1,pver
             qcsinksum_rate1ord(k) = zero
             qcsum_rate1ord(k)     = zero
           enddo


!!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.....................................................................................................
           do it=1,iter ! big iter loop

! initialize sub-step microphysical tendencies

             do k=1,pver
               tlat(i,k)    = zero
               qvlat(i,k)   = zero
               qctend(i,k)  = zero
               qitend(i,k)  = zero
               qnitend(i,k) = zero
               qrtend(i,k)  = zero
               nctend(i,k)  = zero
               nitend(i,k)  = zero
               nrtend(i,k)  = zero
               nstend(i,k)  = zero

! initialize diagnostic precipitation to zero

               qniic(i,k)   = zero
               qric(i,k)    = zero
               nsic(i,k)    = zero
               nric(i,k)    = zero

               rainrt(i,k)  = zero

             enddo


! begin new i,k loop, calculate new cldmax after adjustment to cldm above

! initialize vertically-integrated rain and snow tendencies

             qrtot = zero
             nrtot = zero
             qstot = zero
             nstot = zero

! initialize precip at surface

             prect(i) = zero
             preci(i) = zero

             do k=1,pver

               km = k - 1
! set cwml and cwmi to current qc and qi
  
               cwml(i,k) = qc(i,k)
               cwmi(i,k) = qi(i,k)

! initialize precip fallspeeds to zero

               ums(k) = zero
               uns(k) = zero
               umr(k) = zero
               unr(k) = zero

! calculate precip fraction based on maximum overlap assumption

               if (k == 1) then
                 cldmax(i,k) = cldm(i,k)
               else
! if rain or snow mix ratio is smaller than
! threshold, then set cldmax to cloud fraction at current level
                 if (qric(i,km)  >= qsmall .or.
     &               qniic(i,km) >= qsmall) then
                   cldmax(i,k) = max(cldmax(i,km),cldm(i,k))
                 else
                   cldmax(i,k) = cldm(i,k)
                 end if
               end if

               rdz = rho(i,k)  * dz(i,k)
               rdzi = one / rdz

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%

               if (dum2i(i,k) > zero .and. t(i,k) < (tmelt-five) .and.
     &            relhum(i,k)*esl(i,k)/esi(i,k) > rhmini+0.05_r8) then

!if NCAI > 0. then set numice = ncai (as before)
!note: this is gridbox averaged

                 nimax     = dum2i(i,k)*icldm(i,k)
!                nnuccd(k) = (dum2i(i,k)-ni(i,k)/icldm(i,k))/deltat*icldm(i,k)
                 nnuccd(k) = max(zero, (nimax-ni(i,k))*dti)

!Calc mass of new particles using new crystal mass...
!also this will be multiplied by mtime as nnuccd is...

                 mnuccd(k) = nnuccd(k) * mi0

!  add mnuccd to cmei....
                 cmei(i,k)= cmei(i,k) + mnuccd(k) * mtime

!     if (lprint .and. k == 29) write(0,*)' cmei5=',cmei(i,k)

!  limit cmei

                 qvi = epsqs*esi(i,k)/(p(i,k)-omeps*esi(i,k))
                 dqsidt    = xxls*qvi/(rv*t(i,k)*t(i,k))
                 abi       = one + dqsidt*xxls/cpp
                 cmei(i,k) = min(cmei(i,k),(q(i,k)-qvi)/(abi*deltat))

! limit for roundoff error
                 cmei(i,k) = cmei(i,k)*omsm

!     if (lprint .and. k == 29) write(0,*)' cmei6=',cmei(i,k)
               else
                 nnuccd(k) = zero
                 nimax     = zero
                 mnuccd(k) = zero
               end if


!c............................................................................
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
! for microphysical process calculations
! units are kg/kg for mixing ratio, 1/kg for number conc

! limit in-cloud values to 0.005 kg/kg

               tx1       = one / lcldm(i,k)
               tx2       = one / icldm(i,k)
               qcic(i,k) = min(cwml(i,k)*tx1, 5.e-3_r8)
               qiic(i,k) = min(cwmi(i,k)*tx2, 5.e-3_r8)
               ncic(i,k) = max(nc(i,k)*tx1, zero)
               niic(i,k) = max(ni(i,k)*tx2, zero)


               if (nccons) then
                 ncic(i,k) = ncnst*irho(i,k)
               end if

               if (nicons) then
                 niic(i,k) = ninst*irho(i,k)
               end if

               tx1 = qc(i,k) - berg(i,k)*deltat
               if (tx1 < qsmall) then
                 qcic(i,k) = zero
                 ncic(i,k) = zero
                 if (tx1 < zero) then
                   berg(i,k) = qc(i,k)*dti*omsm
                 end if
               end if

               tx1 = qi(i,k) + (cmei(i,k)+berg(i,k))*deltat
               if (tx1 < qsmall) then
                 qiic(i,k) = zero
                 niic(i,k) = zero
                 if (tx1 < zero) then
                   cmei(i,k) = (-qi(i,k)*dti-berg(i,k))*omsm
!     if (lprint .and. k == 29) write(0,*)' cmei7=',cmei(i,k)
                 end if
               end if

! add to cme output

               cmeout(i,k) = cmeout(i,k) + cmei(i,k)

! decrease in number concentration due to sublimation/evap
! divide by cloud fraction to get in-cloud decrease
! don't reduce Nc due to bergeron process

!Moved it here since nsubi since cmei was not limited before(DONIF 03/13/2015)
               if (cmei(i,k) < zero .and. qi(i,k) > qsmall
     &                              .and. cldm(i,k) > mincld) then
                 nsubi(k) = cmei(i,k)*ni(i,k) / (qi(i,k)*cldm(i,k))
               else
                 nsubi(k) = zero
               end if

               nsubc(k) = zero




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! get size distribution parameters based on in-cloud cloud water/ice
! these calculations also ensure consistency between number and mixing ratio
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!......................................................................
! cloud ice

               if (qiic(i,k) >= qsmall) then

! add upper limit to in-cloud number concentration to prevent numerical error
                 niic(i,k) = min(niic(i,k),qiic(i,k)*1.e20_r8)


                 lami(k) = (cons1*ci*niic(i,k)/qiic(i,k))**oneodi

!      miu_ice(k)=mui_hemp_l(lami(k))
                 miu_ice(k) = max(min(0.008_r8*(lami(k)*0.01)**0.87_r8,
     &                                10.0_r8), 0.1_r8)
                 tx1        = one + miu_ice(k)
                 tx2        = one / gamma(tx1)
                 aux        = (gamma(tx1+di) * tx2) **oneodi
                 lami(k)    = aux*lami(k)

                 n0i(k)     = niic(i,k)*lami(k)**tx1 * tx2

! check for slope
! adjust vars

                 if (lami(k) < lammini) then
                   lami(k)   = lammini
                   n0i(k)    = lami(k)**(di+one)*qiic(i,k)/(ci*cons1)
                   niic(i,k) = n0i(k)/lami(k)
                 else if (lami(k) > lammaxi) then
                   lami(k)   = lammaxi
                   n0i(k)    = lami(k)**(di+one)*qiic(i,k)/(ci*cons1)
                   niic(i,k) = n0i(k)/lami(k)
                 end if

               else
                 lami(k) = zero
                 n0i(k)  = zero
               end if

               if (qcic(i,k) >= qsmall) then

! add upper limit to in-cloud number concentration to prevent numerical error
                 ncic(i,k) = min(ncic(i,k),qcic(i,k)*1.e20_r8)
 
                 ncic(i,k) = max(ncic(i,k),cdnl*irho(i,k))

! get pgam from fit to observations of martin et al. 1994

                 pgam(k) = 0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))
     &                   + 0.2714_r8


                 if (.true.) then
                   if ((ncic(i,k) > 1.0e-3) .and.
     &                 (qcic(i,k) > 1.0e-11)) then
                     xs = 0.07_r8*(1000._r8*qcic(i,k)/ncic(i,k))
     &                                               **(-0.14_r8)
                   else
                     xs = 1.2
                   end if

                   xs = max(min(xs, 1.7_r8), 1.1_r8)
                   xs = xs*xs*xs
                   xs = (xs + sqrt(xs+8.0_r8)*sqrt(xs) - four)/8.0_r8
                   pgam(k) = sqrt(xs)

                 end if



                 pgam(k) = one / (pgam(k)*pgam(k)) - one
                 pgam(k) = max(two, min(15._r8, pgam(k)))


! calculate lamc
                 tx1     = pirhow * gamma(pgam(k)+four)
                 tx3     = gamma(pgam(k)+one)
                 tx2     = 6._r8 * qcic(i,k) * tx3

                 lamc(k) = (tx1*ncic(i,k)/tx2) ** oneb3

! lammin, 50 micron diameter max mean size

                 lammin = (pgam(k)+one) / 50.e-6_r8
                 lammax = (pgam(k)+one) / 2.e-6_r8

                 if (lamc(k) < lammin) then
                   lamc(k)   = lammin
                   ncic(i,k) = tx2*lamc(k)*lamc(k)*lamc(k) / tx1
                 else if (lamc(k) > lammax) then
                   lamc(k)   = lammax
                   ncic(i,k) = tx2*lamc(k)*lamc(k)*lamc(k) / tx1
                 end if

! parameter to calculate droplet freezing

                 cdist1(k) = ncic(i,k) / tx3

               else
                 lamc(k)   = zero
                 cdist1(k) = zero
               end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin micropysical process calculations
!.................................................................
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
! minimum qc of 1 x 10^-8 prevents floating point error

               xs = 0.0

               if (qcic(i,k) >= 1.e-8_r8) then

! nprc is increase in rain number conc due to autoconversion
! nprc1 is decrease in cloud droplet conc due to autoconversion

! assume exponential sub-grid distribution of qc, resulting in additional
! factor related to qcvar below


!                prc(k) = cons2/(cons3*cons18)*1350._r8
!    &                  * qcic(i,k)**2.47_r8
!    &                  * (1.e-6_r8*ncic(i,k)*rho(i,k))**(-1.79_r8)
  
                 if (auto_option == 1) then
  
                   tx1    = qcic(i,k)*rho(i, k)/3.0e-4
                   prc(k) = 1.0e-3 * qcic(i,k) * (one-exp(-tx1*tx1))
     &                    * gamma(one+qcvar)/(cons3*qcvar)

                 elseif (auto_option == 2) then

                   tx1    = qcic(i,k)*rho(i, k)/7.5e-4
                   tx1    = tx1 * tx1
                   prc(k) = 1.0e-4 * qcic(i,k) * (one-exp(-tx1*tx1))
     &                    * gamma(one+qcvar)/(cons3*qcvar)

                 elseif (auto_option == 3) then

                   tx1    = qcic(i,k) / 3.0e-4

                   prc(k) = 1.0e-3 * (one-EXP(-tx1*tx1))
     &                    * cons2/(cons3*cons18)
     &                    / (1.e-6_r8/50.0*ncic(i,k)*rho(i,k))**1.79_r8

                 elseif (auto_option  == 4) then

                   xs = one / (one+pgam(k))

                   beta6 = (one+3.0*xs)*(one+4.0*xs)*(one+5.0*xs)
     &                   / ((one+xs)*(one+xs+xs))

                   LW = 1.0e-3_r8 * qcic(i,k) * rho(i,k)
                   NW = ncic(i,k) * rho(i,k)  * 1.e-6_r8

                   xs = min(20.0, 1.03e16*(LW*LW)/(NW*SQRT(NW)))
                   prc(k) = 1.1e10*beta6*LW*LW*LW
     &                    * (one-exp(-(xs**miu_disp))) / NW
                   prc(k) = prc(k)*1.0e3*irho(i,k)
                   prc(k) = prc(k) * gamma(two+qcvar)
     &                    / (gamma(qcvar)*(qcvar*qcvar))

                   prc(k) = prc(k)*dcrit
  
                   xs = 1/xs
                 else
                   prc(k) = cons2/(cons3*cons18)*1350._r8
     &                    * qcic(i,k)**2.47_r8
     &                    * (1.e-6_r8*ncic(i,k)*rho(i,k))**(-1.79_r8)
                 endif

                 nprc(k) = prc(k)/cons22

                 nprc1(k) = ncic(i,k)*prc(k) / (qcic(i,k)*(one+xs))

               else
                 prc(k)   = zero
                 nprc(k)  = zero
                 nprc1(k) = zero
               endif



! add autoconversion to precip from above to get provisional rain mixing ratio
! and number concentration (qric and nric)

! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)

             if(fprcp==1) then
               qric(i,k) = min(qrn(i,k)/lcldm(i,k), 5.e-3_r8)
               nric(i,k) = max(nrn(i,k)/lcldm(i,k), zero)
             else
               dum  = 0.45_r8
               dum1 = 0.45_r8

               if (k == 1) then
                 tx1 = lcldm(i,k)*dz(i,k)/(cldmax(i,k)*dum)
                 qric(i,k) = prc(k)  * tx1
                 nric(i,k) = nprc(k) * tx1
               else
                 if (qric(i,km) >= qsmall) then
!      dum=umr(k-1)
!      dum1=unr(k-1)
!      Anning Cheng find a possible untable case here
                   dum  = max(umr(km),dum)
                   dum1 = max(unr(km),dum1)
                 endif

! no autoconversion of rain number if rain/snow falling from above
! this assumes that new drizzle drops formed by autoconversion are rapidly collected
! by the existing rain/snow particles from above

                 if (qric(i,km)  >= 1.e-9_r8 .or.
     &               qniic(i,km) >= 1.e-9_r8) then
                   nprc(k) = zero
                 endif

                 tx1 = rho(i,km) * cldmax(i,km)
                 tx3 = rho(i,k)  * cldmax(i,k)
                 qric(i,k) = (tx1*umr(km)*qric(i,km)
     &                  + (rdz*((pra(km)+prc(k))*lcldm(i,k)
     &                  + (pre(km)-pracs(km)-mnuccr(km))*cldmax(i,k))))
     &                  / (dum*tx3)


!      nric(i,k) = (rho(i,k-1)*unr(k-1)*nric(i,k-1)*cldmax(i,k-1)+
!    & (rho(i,k)*dz(i,k)*(nprc(k)*lcldm(i,k)+(nsubr(k-1)-npracs(k-1)-
!    &nnuccr(k-1)+nragg(k-1))*cldmax(i,k)))) /(dum1*rho(i,k)*cldmax(i,
!    &k))

!      Anning nsubr never given a value before
                 nric(i,k) = (tx1*unr(km)*nric(i,km)
     &                + (rdz*(nprc(k)*lcldm(i,k)
     &                +(-npracs(km)-nnuccr(km)+nragg(km))*cldmax(i,k))))
     &                / (dum1*tx3)

               endif
            endif !fprcp

!.......................................................................
! Autoconversion of cloud ice to snow (Ice_aut)
! similar to Ferrier (1994)

               if (t(i,k) <= 273.15_r8 .and. qiic(i,k) >= qsmall) then

! note: assumes autoconversion timescale of 180 sec
!vaux = 180.0_r8*10.0_r8

                 if (.false.) then
                   vaux = ts_auto_ice * 10.0_r8
  
                   nprci(k) = (niic(i, k)/vaux)*exp(-lami(k)*dcs)
                   tx1 = one / lami(k)
                   tx2 = tx1 * tx1
                   prci(k)  = pi*irho(i,k)*niic(i,k)*lami(k)
     &                      / (6._r8*vaux)
     &                      * (cons23*tx1+three*cons24*tx2
     &                      + 6._r8*dcs*tx1*tx2+6._r8*tx2*tx2)
     &                      * exp(-lami(k)*dcs)

                 else

!                  miu_ice(k) = mui_hemp_l(lami(k))
                   miu_ice(k) =max(min(0.008_r8*(lami(k)*0.01)**0.87_r8,
     &                                10.0_r8), 0.1_r8)
                   tx1      = lami(k)*dcs
                   nprci(k) = (niic(i,k)/ts_auto_ice)
     &                      * (one - gamma_incomp(miu_ice(k), tx1))

                   prci(k) = (qiic(i,k)/ts_auto_ice)
     &                     * (one - gamma_incomp(miu_ice(k)+three, tx1))

                 end if
               else
                 prci(k)  = zero
                 nprci(k) = zero
               end if


! add autoconversion to flux from level above to get provisional snow mixing ratio
! and number concentration (qniic and nsic)
! Anning Cheng 9/16/2016 forecasting rain and snow, corresponding to MG2
             if(fprcp==1) then
               qniic(i,k) = min(qsnw(i,k)/icldm(i,k), 5.e-3_r8)
               nsic(i,k) = max(nsnw(i,k)/icldm(i,k), 0._r8)
             else

               dum  = (asn(i,k)*cons25)
               dum1 = (asn(i,k)*cons25)

               if (k == 1) then
                 tx1 = icldm(i,k)*dz(i,k)/(cldmax(i,k)*dum)
                 qniic(i,k) = prci(k)  * tx1
                 nsic(i,k)  = nprci(k) * tx1
               else
                 if (qniic(i,km) >= qsmall) then
                   dum  = ums(km)
                   dum1 = uns(km)
                 end if

                 tx1 = rho(i,km) * cldmax(i,km)
                 tx3 = rho(i,k)  * cldmax(i,k)
                 qniic(i,k) = (tx1*ums(km)*qniic(i,km) + (rdz*
     &            ((prci(k)+prai(km)+psacws(km)+bergs(km))*icldm(i,k)
     &            +(prds(km)+pracs(km)+mnuccr(km))*cldmax(i,k))))
     &            / (dum*tx3)

!      nsic(i,k) = (rho(i,k-1)*uns(k-1)*nsic(i,k-1)*cldmax(i,k-1)+
!    & (rho(i,k)*dz(i,k)*(nprci(k)*icldm(i,k)+(nsubs(k-1)+nsagg(k-1)+
!    &nnuccr(k-1))*cldmax(i,k)))) /(dum1*rho(i,k)*cldmax(i,k))

!      nsubs never given a value before
                 nsic(i,k) = (tx1*uns(km)*nsic(i,km)
     &                     + (rdz*(nprci(k)*icldm(i,k)+(nsagg(km)
     &                     + nnuccr(km))*cldmax(i,k)))) /(dum1*tx3)

               end if
             end if !fprcp

! if precip mix ratio is zero so should number concentration

               if (qniic(i,k) < qsmall) then
                qniic(i,k) = zero
                nsic(i,k)  = zero
               end if

               if (qric(i,k) < qsmall) then
                 qric(i,k) = zero
                 nric(i,k) = zero
               end if

! make sure number concentration is a positive number to avoid
! taking root of negative later

               nric(i,k) = max(nric(i,k),zero)
               nsic(i,k) = max(nsic(i,k),zero)

!.......................................................................
! get size distribution parameters for precip
!......................................................................
! rain

               if (qric(i,k) >= qsmall) then
                 lamr(k) = (pirhow*nric(i,k)/qric(i,k))**oneb3
                 n0r(k)  = nric(i,k)*lamr(k)

! check for slope
! adjust vars

                 if (lamr(k) < lamminr) then
                   lamr(k)   = lamminr
                   tx1       = lamminr * lamminr
                   n0r(k)    = tx1 * tx1 * qric(i,k)/pirhow
                   nric(i,k) = n0r(k)/lamr(k)
                 else if (lamr(k) > lammaxr) then
                   lamr(k)   = lammaxr
                   tx1       = lammaxr * lammaxr
                   n0r(k)    = tx1 * tx1 * qric(i,k)/pirhow
                   nric(i,k) = n0r(k)/lamr(k)
                 end if

! provisional rain number and mass weighted mean fallspeed (m/s)

                 tx1    = arn(i,k) / lamr(k) ** br
                 tx2    = 9.1_r8*rhof(i,k)
                 unr(k) = min(tx1*cons4, tx2)
                 umr(k) = min(tx1*(cons5/6._r8), tx2)

               else
                 lamr(k) = zero
                 n0r(k)  = zero
                 umr(k)  = zero
                 unr(k)  = zero
               end if

!......................................................................
! snow

               if (qniic(i,k) >= qsmall) then
                 lams(k) = (cons6*cs*nsic(i,k) / qniic(i,k))**(one/ds)
                 n0s(k)  = nsic(i,k)*lams(k)

! check for slope
! adjust vars

                 if (lams(k) < lammins) then
                   lams(k)   = lammins
                   n0s(k)    = lams(k)**(ds+one)*qniic(i,k)/(cs*cons6)
                   nsic(i,k) = n0s(k)/lams(k)
                 else if (lams(k) > lammaxs) then
                   lams(k)   = lammaxs
                   n0s(k)    = lams(k)**(ds+one)*qniic(i,k)/(cs*cons6)
                   nsic(i,k) = n0s(k)/lams(k)
                 end if

! provisional snow number and mass weighted mean fallspeed (m/s)

                 tx1    = asn(i,k) / lams(k)**bs
                 tx2    = 1.2_r8*rhof(i,k)
                 ums(k) = min(tx1*(cons8/6._r8), tx2)
                 uns(k) = min(tx1*cons7, tx2)

               else
                 lams(k) = zero
                 n0s(k)  = zero
                 ums(k)  = zero
                 uns(k)  = zero
               end if

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
checked up to here moorthi

! heterogeneous freezing of cloud water

               if (qcic(i,k) >= qsmall .and. t(i,k) < 269.15_r8) then

! immersion freezing (Bigg, 1953)

                 tx2       = one / (lamc(k) * lamc(k) * lamc(k))
                 tx3       = min(aimm*(273.15_r8-t(i,k)), 25.0)
                 tx1       = cdist1(k) * tx2 * bimm * exp(tx3)
                 mnuccc(k) = cons9/(cons3*cons19)* pi*pirhow/36._r8
     &                     * gamma(7._r8+pgam(k))*tx1*tx2

                 nnuccc(k) = cons10/(cons3*qcvar)* pi/6._r8
     &                     * gamma(pgam(k)+four) * tx1

                 if (.true.) then
                   nnuccc(k) = nimm(i,k)
                   mnuccc(k) = nimm(i,k)*qcic(i,k)/max(cdist1(k), one)
                 end if

! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
! dust size and number in 4 bins are read in from companion routine

!DONIFFFF this implementation would oveestimate strongly contact IN
! According to Young 1974 the base concentration must be the active contact IN at -4 C
! This implementation is using all the dust at all levels
!Freezing fraction of collected IN is assumed to be 1.

                 tcnt      = (270.16_r8-t(i,k))**1.3_r8
                 viscosity = 1.8e-5_r8*(t(i,k)/298.0_r8)**0.85_r8
                 mfp       = two*viscosity
     &                     / (p(i,k) *sqrt(8.0_r8*28.96e-3_r8/(pi*
     &                                     8.314409_r8*t(i,k))))
                 taux      = t(i,k) - three
                 taux      = max(taux-273.15, -40.0)

                 nsdust    = max(exp(-0.517*taux + 8.934) -3.76e6, 0.0)
                 nssoot    = max(1.0e4*exp((-0.0101*taux-0.8525)*taux
     &                                      +0.7667) -3.77e9, 0.0)

                 do ind_aux = 1, nbincontactdust
                   tx1 = rndst(i,k,ind_aux)
                   nslip(ind_aux) = one + (mfp/tx1) *
     &                    (1.257_r8+(0.4_r8*Exp(-(1.1_r8*tx1/mfp))))
                   ndfaer(ind_aux) = 1.381e-23_r8*t(i,k)* nslip(ind_aux)
     &                             /  (6._r8*pi*viscosity*tx1)

                   ndfaer(ind_aux) = ndfaer(ind_aux)
     &                             * (1.0-exp(-nsdust*12.5664*tx1*tx1))
                 end do

                 nslipsoot  = one + (mfp/rnsoot(i,k)) *
     &                (1.257_r8+(0.4_r8*Exp(-(1.1_r8*rnsoot(i,k)/mfp))))
                 ndfaersoot = 1.381e-23_r8*t(i,k)*nslipsoot
     &                    / (6._r8*pi*viscosity*rnsoot(i,k))
                 ndfaersoot = ndfaersoot *
     &                (1.0-exp(-nssoot*12.5664*rnsoot(i,k)*rnsoot(i,k)))



                 tx1 = sum(ndfaer(:)*nacon(i,k,:))+ndfaersoot*nsoot(i,k)
                 tx1 = tx1 * pi * cdist1(k)
                 tx2 = one / lamc(k)
                 tx3 = tx2 * tx2
                 mnucct(k) = (pi*cons10/(cons3*qcvar*three)) * tx1
     &                     * rhow*gamma(pgam(k)+five)*tx3*tx3


                 nnucct(k) = (tx1+tx1) * gamma(pgam(k)+two) * tx2

! make sure number of droplets frozen does not exceed available ice nuclei concentration
! this prevents 'runaway' droplet freezing
#ifdef CAM
                 if (nnuccc(k)*lcldm(i,k) > nnuccd(k)) then
                   dum = (nnuccd(k)/(nnuccc(k)*lcldm(i,k)))
! scale mixing ratio of droplet freezing with limit
                   mnuccc(k) = mnuccc(k)*dum
                   nnuccc(k) = nnuccd(k)/lcldm(i,k)
                 end if
#endif


                 if ((nnucct(k)+nnuccc(k))*deltat > ncic(i, k)) then
  
                   tx1 = tx2 * tx3
                   nnuccc(k) = ncic(i, k)*dti
                   mnuccc(k) = cons9/(cons3*cons19)* pi*pirhow/36._r8
     &                       * cdist1(k)*gamma(7._r8+pgam(k))*ncic(i,k)
     &                       * tx1 * tx1

                   nnucct(k) = zero
                   mnucct(k) = zero
                 end if


               else
                 mnuccc(k) = zero
                 nnuccc(k) = zero
                 mnucct(k) = zero
                 nnucct(k) = zero
               end if

!.......................................................................
! snow self-aggregation from passarelli, 1978, used by reisner, 1998
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

               if (qniic(i,k) >= qsmall .and. t(i,k) <= 273.15_r8) then
                 tx1 = (two+bs)*oneb3
                 nsagg(k) = -1108._r8*asn(i,k)*Eii* pi**((one-bs)*oneb3)
     &                    * (rho(i,k)*qniic(i,k)/rhosn) ** tx1
     &                    * (nsic(i,k)*rho(i,k))**((four-bs)*oneb3)
     &                    / (four*720._r8*rho(i,k))
               else
                 nsagg(k) = zero
               end if

!.......................................................................
! accretion of cloud droplets onto snow/graupel
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice
! since minimum size ice particle for accretion is 50 - 150 micron

! ignore collision of snow with droplets above freezing

               if (qniic(i,k) >= qsmall .and. t(i,k) <= tmelt
     &                                  .and. qcic(i,k) >= qsmall) then

! put in size dependent collection efficiency
! mean diameter of snow is area-weighted, since
! accretion is function of crystal geometric area
! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

                 dc0 = (pgam(k)+one)/lamc(k)
                 ds0 = one/lams(k)
                 dum = dc0*dc0*uns(k)*rhow/(9._r8*mu(i,k)*ds0)
                 eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

                 eci = min(one, max(eci, zero))


! no impact of sub-grid distribution of qc since psacws
! is linear in qc

                 tx1      = pi/four*asn(i,k)*rho(i,k)*n0s(k)*Eci*cons11
     &                    / lams(k)**(bs+three)
                 psacws(k)  = tx1 * qcic(i,k)
                 npsacws(k) = tx1 * ncic(i,k)
               else
                 psacws(k)  = zero
                 npsacws(k) = zero
               end if

! add secondary ice production due to accretion of droplets by snow
! (Hallet-Mossop process) (from Cotton et al., 1986)

               if((t(i,k) < 270.16_r8) .and. (t(i,k) >= 268.16_r8)) then
                 ni_secp   = 0.5*3.5e8_r8*(270.16_r8-t(i,k))*psacws(k)
                 nsacwi(k) = ni_secp
                 msacwi(k) = min(ni_secp*mi0, psacws(k))
               else if((t(i,k) < 268.16_r8) .and.
     &                 (t(i,k) >= 265.16_r8)) then
                 ni_secp   = oneb3*3.5e8_r8*(t(i,k)-265.16_r8)*psacws(k)
                 nsacwi(k) = ni_secp
                 msacwi(k) = min(ni_secp*mi0, psacws(k))
               else
                 ni_secp   = zero
                 nsacwi(k) = zero
                 msacwi(k) = zero
               endif
               psacws(k) = max(zero, psacws(k)-ni_secp*mi0)

!.......................................................................
! accretion of rain water by snow
! formula from ikawa and saito, 1991, used by reisner et al., 1998

               if (qric(i,k) >= 1.e-8_r8 .and. qniic(i,k) >= 1.e-8_r8
     &          .and. t(i,k) <= 273.15_r8) then
! Anning decrease pracs it can reach 2.3e5 so/1.e6
                 tx1 = 1.2_r8*umr(k) - 0.95_r8*ums(k)
                 tx2 = sqrt(tx1*tx1+0.08_r8*ums(k)*umr(k))
                 tx1 = one / lamr(k)
                 tx3 = one / lams(k)
                 tx4 = tx1 * tx1
                 tx5 = pi * ecr * rho(i,k) *n0r(k) * n0s(k)

                 pracs(k) = pirhow*tx2*tx5 * 
     &          (tx4*tx4*tx3*(five*tx4+tx3*(two*tx1+half*tx3)))


                 tx2 = unr(k) - uns(k)
                 tx2 = sqrt(1.7_r8*tx2*tx2 + 0.3_r8*unr(k)*uns(k)) 

                 npracs(k) = half*tx2*tx5*tx1*tx3*(tx4+tx3*(tx1+tx3))

               else
                 pracs(k)  = zero
                 npracs(k) = zero
               end if

!.......................................................................
! heterogeneous freezing of rain drops
! follows from Bigg (1953)

!      if (t(i,k).lt.269.15_r8 .and. qric(i,k).ge.qsmall) then
!      Anning change to prevent huge value of mnuccr

               if (t(i,k) < 269.15_r8 .and. qric(i,k) >= qsmall
     &                                .and. t(i,k) > 260.0_r8) then

!                tx1 = exp(min(aimm*(273.15_r8-t(i,k))-one, 50.0))
                 tx1 = exp(min(aimm*(273.15_r8-t(i,k))-one, 25.0))
                 tx2 = 1.0 / (lamr(k)*lamr(k)*lamr(k))

                 nnuccr(k) = pi * nric(i,k) * bimm * tx1 * tx2
                 mnuccr(k) = 20._r8 * pirhow * nnuccr(k) * tx2

               else
                 mnuccr(k) = zero
                 nnuccr(k) = zero

               end if

!.......................................................................
! accretion of cloud liquid water by rain
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected

               if (qric(i,k) >= qsmall .and. qcic(i,k) >= qsmall) then

! include sub-grid distribution of cloud water

                 pra(k)  = cons12/(cons3*cons20)
     &                   * 67._r8*(qcic(i,k)*qric(i,k))**1.15_r8
                 npra(k) = pra(k) * (ncic(i,k)/qcic(i,k))

               else
                 pra(k)  = zero
                 npra(k) = zero
               end if

!Not done above
!.......................................................................
! Self-collection of rain drops
! from Beheng(1994)
  
               if (qric(i,k) >= qsmall) then
               nragg(k) = -8._r8*nric(i,k)*qric(i,k)*rho(i,k)
               else
                 nragg(k) = zero
               end if

!.......................................................................
! Accretion of cloud ice by snow
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection

               if (qniic(i,k) >= qsmall .and. qiic(i,k) >= qsmall
     &                                  .and.t(i,k) <= 273.15_r8) then

                 tx1     = (pi/four)*asn(i,k)*rho(i,k)*n0s(k)*Eii*cons11
     &                   / lams(k)**(bs+three)
                 prai(k)  = tx1 * qiic(i,k)
                 nprai(k) = tx1 * niic(i,k)

                 nprai(k)= min(nprai(k), 1.0e10)
  
               else
                 prai(k)  = zero
                 nprai(k) = zero
               end if

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate evaporation/sublimation of rain and snow
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

! initialize evap/sub tendncies
               pre(k)  = zero
               prds(k) = zero

! evaporation of rain
! only calculate if there is some precip fraction > cloud fraction

               if (qcic(i,k)+qiic(i,k) < 1.e-6_r8 .or.
     &             cldmax(i,k) > lcldm(i,k)) then

! set temporary cloud fraction to zero if cloud water + ice is very small
! this will ensure that evaporation/sublimation of precip occurs over
! entire grid cell, since min cloud fraction is specified otherwise
                 if (qcic(i,k)+qiic(i,k) < 1.e-6_r8) then
                   dum = zero
                 else
                   dum = lcldm(i,k)
                 end if

! saturation vapor pressure
!                esn = polysvp(t(i,k),0)
                 esn = min(fpvsl(t(i,k)), p(i,k))
                 qsn = min(epsqs*esn/(p(i,k)-omeps*esn), one)

! recalculate saturation vapor pressure for liquid and ice
                 esl(i,k) = esn
!                esi(i,k) = polysvp(t(i,k),1)
                 esi(i,k) = min(fpvsi(t(i,k)), p(i,k))
! hm fix, make sure when above freezing that esi=esl, not active yet
                 if (t(i,k) > tmelt) esi(i,k) = esl(i,k)

! calculate q for out-of-cloud region

                 qclr = (q(i,k)-dum*qsn) / (one-dum)

                 if (qric(i,k) >= qsmall) then

                   qvs   = epsqs*esl(i,k)/(p(i,k)-omeps*esl(i,k))
                   dqsdt = xxlv*qvs/(rv*t(i,k)*t(i,k))
                   ab    = one + dqsdt*xxlv/cpp
                   epsr = (pi+pi)*n0r(k)*rho(i,k)*Dv(i,k)
     &                  * (f1r/(lamr(k)*lamr(k))
     &                  +  f2r*sqrt(arn(i,k)*rho(i,k)/mu(i,k))
     &                  *  sc(i,k)**oneb3
     &                  *  cons13 / lamr(k)**((five+br)*half))

                   pre(k) = epsr*(qclr-qvs)/ab

! only evaporate in out-of-cloud region
! and distribute across cldmax
                   pre(k) = min(pre(k)*(cldmax(i,k)-dum), zero)
                   pre(k) = pre(k) / cldmax(i,k)

                 end if

! sublimation of snow
                 if (qniic(i,k) >= qsmall) then
                   qvi    = epsqs*esi(i,k)/(p(i,k)-omeps*esi(i,k))
                   dqsidt = xxls*qvi/(rv*t(i,k)*t(i,k))
                   abi    = one + dqsidt*xxls/cpp
                   epss   = (pi+pi)*n0s(k)*rho(i,k)*Dv(i,k)
     &                    * (f1s/(lams(k)*lams(k))
     &                    +  f2s*sqrt((asn(i,k)*rho(i,k)/mu(i,k)))
     &                    *  sc(i,k)**oneb3
     &                    *  cons14/ lams(k)**((five+bs)*half))
                   prds(k) = epss*(qclr-qvi)/abi

! only sublimate in out-of-cloud region and distribute over cldmax
                   prds(k) = min(prds(k)*(cldmax(i,k)-dum), zero)
                   prds(k) = prds(k)/cldmax(i,k)
                 end if


! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
! get updated RH at end of time step based on cloud water/ice condensation/evap

                 tx1  = pre(k)  * cldmax(i,k)
                 tx2  = prds(k) * cldmax(i,k)
                 qtmp = q(i,k) - (cmei(i,k)+tx1+tx2) * deltat
                 ttmp = t(i,k) + (tx1*xxlv + (cmei(i,k)+tx2)*xxls)
     &                         * (deltat/cpp)

!limit range of temperatures!
                 ttmp = max(180._r8,min(ttmp,323._r8))

!                esn = polysvp(ttmp,0)
                 esn = min(fpvsl(ttmp), p(i,k))
                 qsn = min(epsqs*esn/(p(i,k)-omeps*esn), one)

! modify precip evaporation rate if q > qsat
                 if (qtmp > qsn) then
                   if (pre(k)+prds(k) < -1.e-20) then
                     dum1 = pre(k) / (pre(k)+prds(k))
! recalculate q and t after cloud water cond but without precip evap
                     qtmp = q(i,k) - cmei(i,k)*deltat
                     ttmp = t(i,k) + cmei(i,k)*xxls*deltat/cpp
!                    esn = polysvp(ttmp,0)
                     esn = min(fpvsl(ttmp), p(i,k))
                     qsn = min(epsqs*esn/(p(i,k)-omeps*esn), one)
                     tx1 = one / (cpp*rv*ttmp*ttmp)
                     dum = min(zero, (qtmp-qsn)/(one+cons27*qsn*tx1))

! modify rates if needed, divide by cldmax to get local (in-precip) value
                     pre(k) = dum*dum1/(deltat*cldmax(i,k))

! do separately using RHI for prds....
!                    esn = polysvp(ttmp,1)
                     esn = min(fpvsi(ttmp), p(i,k))
                     qsn = min(epsqs*esn/(p(i,k)-omeps*esn), one)
                     dum = min(zero, (qtmp-qsn)/(one+cons28*qsn*tx1))

! modify rates if needed, divide by cldmax to get local (in-precip) value
                     prds(k) = dum*(one-dum1)/(deltat*cldmax(i,k))
                   end if
                 end if

               end if

! bergeron process - evaporation of droplets and deposition onto snow
               if (qniic(i,k) >= qsmall .and. qcic(i,k) >= qsmall
     &                                  .and. t(i,k) < tmelt) then
                 qvi = epsqs*esi(i,k)/(p(i,k)-omeps*esi(i,k))
                 qvs = epsqs*esl(i,k)/(p(i,k)-omeps*esl(i,k))
                 dqsidt = xxls*qvi/(rv*t(i,k)*t(i,k))
                 abi = one + dqsidt*xxls/cpp
                 epss = (pi+pi)*n0s(k)*rho(i,k)*Dv(i,k)
     &                * (f1s/(lams(k)*lams(k))
     &                 + f2s*sqrt(asn(i,k)*rho(i,k)/mu(i,k))
     &                 * sc(i,k)**oneb3
     &                 * cons14 / (lams(k)**((five+bs)*half)))
                 bergs(k) = epss*(qvs-qvi)/abi
               else
                 bergs(k) = zero 
               end if



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! conservation to ensure no negative values of cloud water/precipitation
! in case microphysical process rates are large

! make sure and use end-of-time step values for cloud water, ice, due
! condensation/deposition

! note: for check on conservation, processes are multiplied by omsm
! to prevent problems due to round off error

! include mixing timescale  (mtime)

               qce = max(qc(i,k)-berg(i,k)*deltat, zero)
               nce = nc(i,k) + npccn(i,k)*deltat*mtime
               qie = qi(i,k) + (cmei(i,k)+berg(i,k))*deltat
               nie = ni(i,k) + nnuccd(k)*deltat*mtime

! conservation of qc

               tx1 = lcldm(i,k) * deltat
               dum = (prc(k) + pra(k) + mnuccc(k) + mnucct(k)
     &             +  msacwi(k) + psacws(k) + bergs(k)) * tx1

               if (dum > qce) then

                 ratio = qce/dum * omsm

                 prc(k)    = prc(k)    * ratio
                 pra(k)    = pra(k)    * ratio
                 mnuccc(k) = mnuccc(k) * ratio
                 mnucct(k) = mnucct(k) * ratio
                 msacwi(k) = msacwi(k) * ratio
                 psacws(k) = psacws(k) * ratio
                 bergs(k)  = bergs(k)  * ratio
               end if

! conservation of nc


               dum = (nprc1(k) + npra(k) + nnuccc(k) + nnucct(k)
     &              + npsacws(k) - nsubc(k)) * tx1

               if (dum > nce) then
                 ratio = nce/dum * omsm

                 nprc1(k)   = nprc1(k)   * ratio
                 npra(k)    = npra(k)    * ratio
                 nnuccc(k)  = nnuccc(k)  * ratio
                 nnucct(k)  = nnucct(k)  * ratio
                 npsacws(k) = npsacws(k) * ratio
                 nsubc(k)   = nsubc(k)   * ratio
               end if

! conservation of qi

               dum = ((-mnuccc(k)-mnucct(k)-msacwi(k))*lcldm(i,k)
     &               + (prci(k)+prai(k))*icldm(i,k)) * deltat

               if (dum > qie) then
                 if (prci(k)+prai(k)  > zero) then
                   ratio = (qie*dti
     &                    +(mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(i,k))
     &                   / ((prci(k)+prai(k))*icldm(i,k))*omsm
                 else
                   ratio = zero
                 end if

                 prci(k) = prci(k)*ratio
                 prai(k) = prai(k)*ratio
               end if

! conservation of ni

               dum = ((-nnucct(k)-nsacwi(k))*lcldm(i,k)
     &             +  (nprci(k)+ nprai(k)-nsubi(k))*icldm(i,k))*deltat

               if (dum > nie) then
                 if ( abs(nprci(k)+nprai(k)-nsubi(k)) > zero) then
                   ratio = (nie*dti+(nnucct(k)+nsacwi(k))*lcldm(i,k))
     &                  / ((nprci(k)+nprai(k)-nsubi(k))*icldm(i,k))*omsm
                 else
  
                   ratio = zero
                 end if
                 nprci(k) = nprci(k) * ratio
                 nprai(k) = nprai(k) * ratio
                 nsubi(k) = nsubi(k) * ratio
               end if

! for preciptiation conservation, use logic that vertical integral
! of tendency from current level to top of model (i.e., qrtot) cannot be negative

! conservation of rain mixing rat

               if ( ((prc(k)+pra(k))*lcldm(i,k)
     &             + (-mnuccr(k)+pre(k)-pracs(k))*cldmax(i,k))*rdz
     &             + qrtot < zero) then

                 if (-pre(k)+pracs(k)+mnuccr(k) >= qsmall) then

                   ratio = (qrtot*rdzi + (prc(k)+pra(k))*lcldm(i,k))
     &                / ((-pre(k)+pracs(k)+mnuccr(k))*cldmax(i,k))*omsm

                   pre(k)    = pre(k)    * ratio
                   pracs(k)  = pracs(k)  * ratio
                   mnuccr(k) = mnuccr(k) * ratio
                 end if
               end if

! conservation of nr -  for now neglect evaporation of nr

               nsubr(k) = zero

               if ((nprc(k)*lcldm(i,k)+(-nnuccr(k)+nsubr(k)-npracs(k) +
     &             nragg(k))*cldmax(i,k))*rdz + nrtot < zero) then

                 if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k) >= qsmall)
     &                                                         then
                   ratio = (nrtot*rdzi+nprc(k)*lcldm(i,k))
     &                   / ((-nsubr(k)-nragg(k)+npracs(k)
     &                       +nnuccr(k))*cldmax(i,k))*omsm

                   nsubr(k)  = nsubr(k)  * ratio
                   npracs(k) = npracs(k) * ratio
                   nnuccr(k) = nnuccr(k) * ratio
                   nragg(k)  = nragg(k)  * ratio
                 end if
               end if

! conservation of snow mix ratio

               tx1 = (bergs(k)+psacws(k))*lcldm(i,k)
     &             + (prai(k)+prci(k))*icldm(i,k)
               if ((tx1+(pracs(k)+ mnuccr(k)+prds(k))*cldmax(i,k))
     &                  *rdz+qstot < zero) then

                 if (-prds(k) >= qsmall) then

                   ratio = (qstot*rdzi + tx1
     &                   + (pracs(k)+mnuccr(k))*cldmax(i,k))
     &                   / (-prds(k)*cldmax(i,k))*omsm

                   prds(k) = prds(k) * ratio
                 end if
               end if

! conservation of ns

! calculate loss of number due to sublimation
! for now neglect sublimation of ns
               nsubs(k) = zero

               if ((nprci(k)*icldm(i,k)
     &           +(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax(i,k))
     &           * rdz + nstot < zero) then

                 if (-nsubs(k)-nsagg(k) >= qsmall) then

                   ratio = (nstot*rdzi
     &                    + nprci(k)*icldm(i,k)+ nnuccr(k)*cldmax(i,k))
     &                    / ((-nsubs(k)-nsagg(k))*cldmax(i,k))*omsm

                   nsubs(k) = nsubs(k) * ratio
                   nsagg(k) = nsagg(k) * ratio
                 end if
               end if

! get tendencies due to microphysical conversion processes
! note: tendencies are multiplied by appropaiate cloud/precip
! fraction to get grid-scale values
! note: cmei is already grid-average values

               qvlat(i,k) = qvlat(i,k)
     &                    - (pre(k)+prds(k))*cldmax(i,k) - cmei(i,k)
!     if (lprint .and. k == 29) write(0,*)' qvlata=',qvlat(i,k),
!    &' pre=',pre(k),' prds=',prds(k),' cldmax=',cldmax(i,k),cmei(i,k)
!    &,' it=',it

               tlat(i,k) = tlat(i,k)
     &                   +  pre(k)*cldmax(i,k) * xxlv
     &                   + (prds(k)*cldmax(i,k)+cmei(i,k)) * xxls
     &                   + ((bergs(k)+psacws(k)+mnuccc(k)+
     &                       mnucct(k)+msacwi(k))*lcldm(i,k)
     &                    + (mnuccr(k)+ pracs(k))*cldmax(i,k)
     &                      + berg(i,k)) * xlf


!      if(xlon<0.01.and.xlon>-0.01.and.xlat>1.346        
!    & .and.xlat<1.347.and.k==38)                        
!    & write(*,*)"anning_m0",pre(k),prds(k),cmei(i,k),
!    & bergs(k),psacws(k),
!    & mnuccc(k),mnucct(k),msacwi(k),mnuccr(k),pracs(k),berg(i,k)

               qctend(i,k) = qctend(i,k)
     &                     + (-pra(k)-prc(k)-mnuccc(k)-mnucct(k)
     &                        -msacwi(k)-psacws(k)-bergs(k))*lcldm(i,k)
     &                     - berg(i,k)

               qitend(i,k) = qitend(i,k)
     &                     + (mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(i,k)
     &                     + (-prci(k)-prai(k))*icldm(i,k)
     &                     + cmei(i,k) + berg(i,k)

               qrtend(i,k) = qrtend(i,k)
     &                     + (pra(k)+prc(k))*lcldm(i,k)
     &                     + (pre(k)-pracs(k)-mnuccr(k))*cldmax(i,k)
       
               qnitend(i,k) = qnitend(i,k)
     &                      + (prai(k)+prci(k))*icldm(i,k)
     &                      + (psacws(k)+bergs(k))*lcldm(i,k)
     &                      + (prds(k)+pracs(k)+mnuccr(k))*cldmax(i,k)

!     if (lprint) write(0,*)' k=',k,' qnitend=',qnitend(i,k),
!    & prai(k), prci(k), icldm(i,k),psacws(k),bergs(k),lcldm(i,k)
!    &,prds(k),pracs(k),mnuccr(k),' cldmax=',cldmax(i,k) 

! add output for cmei (accumulate)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               cmeiout(i,k) = cmeiout(i,k) + cmei(i,k)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! assign variables for trop_mozart, these are grid-average
! evaporation/sublimation is stored here as positive term

               evapsnow(i,k) = evapsnow(i,k) + prds(k) * cldmax(i,k)
               nevapr(i,k)   = nevapr(i,k)   + pre(k)  * cldmax(i,k)

! change to make sure prain is positive: do not remove snow from
! prain used for wet deposition
               prain(i,k)    = prain(i,k)
     &                       + (pra(k)+prc(k))*lcldm(i,k)
     &                       + (-pracs(k)-mnuccr(k))*cldmax(i,k)

               prodsnow(i,k) = prodsnow(i,k)
     &                       + (prai(k)+prci(k))*icldm(i,k)
     &                       + (psacws(k)+bergs(k))*lcldm(i,k)
     &                       + (pracs(k)+mnuccr(k))*cldmax(i,k)

! following are used to calculate 1st order conversion rate of cloud water
!    to rain and snow (1/s), for later use in aerosol wet removal routine
! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
!    used to calculate pra, prc, ... in this routine
! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
!                      (no cloud ice or bergeron terms)
! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }

               qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k)
     &                             +(pra(k)+prc(k)+psacws(k))*lcldm(i,k)
               qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc(i,k)

! microphysics output, note this is grid-averaged
               prao(i,k)     = prao(i,k)     + pra(k)    * lcldm(i,k)
               prco(i,k)     = prco(i,k)     + prc(k)    * lcldm(i,k)
               mnuccco(i,k)  = mnuccco(i,k)  + mnuccc(k) * lcldm(i,k)
               mnuccto(i,k)  = mnuccto(i,k)  + mnucct(k) * lcldm(i,k)
               msacwio(i,k)  = msacwio(i,k)  + msacwi(k) * lcldm(i,k)
               psacwso(i,k)  = psacwso(i,k)  + psacws(k) * lcldm(i,k)
               bergso(i,k)   = bergso(i,k)   + bergs(k)  * lcldm(i,k)
               bergo(i,k)    = bergo(i,k)    + berg(i,k)
               prcio(i,k)    = prcio(i,k)    + prci(k)   * icldm(i,k)
               praio(i,k)    = praio(i,k)    + prai(k)   * icldm(i,k)
               mnuccro(i,k)  = mnuccro(i,k)  + mnuccr(k) * cldmax(i,k)
               pracso (i,k)  = pracso (i,k)  + pracs (k) * cldmax(i,k)

               mnuccdo(i,k)  = mnuccdo(i,k)  + mnuccd(k)
               nnuccto(i,k)  = nnuccto(i,k)  + nnucct(k) * lcldm(i,k)

               nnuccdo(i,k)  = nnuccdo(i,k)  + nnuccd(k)
               nnuccco(i,k)  = nnuccco(i,k)  + nnuccc(k) * lcldm(i,k)
               nsacwio(i,k)  = nsacwio(i,k)  + nsacwi(k) * icldm(i,k)
               nsubio(i,k)   = nsubio(i,k)   + nsubi(k)  * icldm(i,k)
               nprcio(i,k)   = nprcio(i,k)   + nprci(k)  * icldm(i,k)
               npraio(i,k)   = npraio(i,k)   + nprai(k)  * icldm(i,k)

               npccno(i, k)  = npccno(i,k)   + npccn(i,k)*lcldm(i,k)
               npsacwso(i,k) = npsacwso(i,k) + npsacws(k)*lcldm(i,k)
               nsubco(i,k)   = nsubco(i,k)   + nsubc(k)  * lcldm(i,k)
               nprao(i,k)    = nprao(i,k)    + npra(k)   * lcldm(i,k)
               nprc1o(i,k)   = nprc1o(i,k)   + nprc1(k)  * lcldm(i,k)


! multiply activation/nucleation by mtime to account for fast timescale

               nctend(i,k) = nctend(i,k) + npccn(i,k)*mtime
     &                     + (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k)
     &                        -npra(k)-nprc1(k))*lcldm(i,k)

               nitend(i,k) = nitend(i,k) + nnuccd(k)*mtime
     &                     + (nnuccc(k)+nnucct(k)+nsacwi(k))*lcldm(i,k)
     &                     + (nsubi(k)-nprci(k)- nprai(k))*icldm(i,k)

               nstend(i,k) = nstend(i,k)
     &                     + (nsubs(k)+ nsagg(k)+nnuccr(k))*cldmax(i,k)
     &                     + nprci(k)*icldm(i,k)

               nrtend(i,k) = nrtend(i,k) + nprc(k)*lcldm(i,k)
     &                     + (nsubr(k)-npracs(k)-nnuccr(k)+nragg(k))
     &                     * cldmax(i,k)

! make sure that nc and ni at advanced time step do not exceed
! maximum (existing N + source terms*dt), which is possible due to
! fast nucleation timescale

               if (nctend(i,k) > zero .and.
     &             nc(i,k)+nctend(i,k)*deltat > ncmax(i,k)) then
                 nctend(i,k) = max(zero, (ncmax(i,k)-nc(i,k))*dti)
               end if
               if (nitend(i,k) > zero .and.
     &             ni(i,k)+nitend(i,k)*deltat > nimax) then
                 nitend(i,k) = max(zero, (nimax-ni(i,k))*dti)
               end if


! get final values for precipitation q and N, based on
! flux of precip from above, source/sink term, and terminal fallspeed
! see eq. 15-16 in MG2008
            if(fprcp==0) then
! rain

               if (qric(i,k) >= qsmall) then
                 if (k == 1) then
                   qric(i,k) = qrtend(i,k)*dz(i,k)/(cldmax(i,k)*umr(k))
                   nric(i,k) = nrtend(i,k)*dz(i,k)/(cldmax(i,k)*unr(k))
                 else
                   tx1       = rho(i,km) * cldmax(i,km)
                   tx3       = rho(i,k)  * cldmax(i,k)

                   qric(i,k) = (tx1*umr(km)*qric(i,km)
     &                       +  rdz*qrtend(i,k)) / (umr(k)*tx3)

                   nric(i,k) = (tx1*unr(km)*nric(i,km)
     &                       +  rdz*nrtend(i,k)) / (unr(k)*tx3)

                 end if
               else
                 qric(i,k) = zero
                 nric(i,k) = zero
               end if

! snow

               if (qniic(i,k) >= qsmall) then
                 if (k == 1) then
                   tx1        = dz(i,k)/cldmax(i,k)
                   qniic(i,k) = qnitend(i,k)*tx1/ums(k)
                   nsic(i,k)  = nstend(i,k)*tx1/uns(k)
                 else
                   tx1        = rho(i,km) * cldmax(i,km)
                   tx3        = rho(i,k)  * cldmax(i,k)

                   qniic(i,k) = (tx1*ums(km)*qniic(i,km)
     &                        +  rdz*qnitend(i,k)) / (ums(k)*tx3)
                   nsic(i,k)  = (tx1*uns(km)*nsic(i,km)
     &                        +  rdz*nstend(i,k))  / (uns(k)*tx3)
                 end if
               else
                 qniic(i,k) = zero
                 nsic(i,k)  = zero
               end if

! calculate precipitation flux at surface
! divide by density of water to get units of m/s

               tx1      = rdz/rhow
               prect(i) = prect(i) + (qrtend(i,k)+ qnitend(i,k))*tx1
               preci(i) = preci(i) +  qnitend(i,k)*tx1

!     if (lprint) write(0,*)' prect=',prect(i),' preci=',preci(i)
!    &,' qrtend=',qrtend(i,k),' qnitend=',qnitend(i,k),' rdz=',rdz
!    &,' k=',k,' it=',it, 'rhow=',rhow,' rho=',rho(i,k),' dz=',dz(i,k)



               rainrt(i,k) = qric(i,k)*rho(i,k)*umr(k)
     &                     / rhow*3600._r8*1000._r8

! vertically-integrated precip source/sink terms (note: grid-averaged)

               qrtot = max(qrtot + qrtend(i,k)*rdz,  zero)
               qstot = max(qstot + qnitend(i,k)*rdz, zero)
               nrtot = max(nrtot + nrtend(i,k)*rdz,  zero)
               nstot = max(nstot + nstend(i,k)*rdz,  zero)

!!!! done up to here - Moorthi
! calculate melting and freezing of precip
! melt snow at +2 C
               taux = 1.0

               if (.true.) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 tx1 = t(i,k) + tlat(i,k)*onebcp*deltat
                 if (tx1 > 275.15_r8) then

                   if (qstot > zero) then

! make sure melting snow doesn't reduce temperature below threshold
                     dum = xlfocp*qstot*rdzi
                     if (tx1-dum*deltat < 273.15_r8) then
                       tx2 = (tx1-275.15_r8) * dti
                       dum = min(one,max(zero,tx2/dum))
                     else
                       dum = one
                     end if

                     qric(i,k)  = qric(i,k) + dum*qniic(i,k)
                     nric(i,k)  = nric(i,k) + dum*nsic(i,k)
                     qniic(i,k) = (one-dum)*qniic(i,k)
                     nsic(i,k)  = (one-dum)*nsic(i,k)
! heating tendency
                     tmp = -xlf*dum*qstot*rdzi
                     meltsdt(i,k) = meltsdt(i,k) + tmp

                     tlat(i,k) = tlat(i,k) + tmp

!      if(xlon<0.01.and.xlon>-0.01.and.xlat>1.346        
!    & .and.xlat<1.347.and.k==38)                        
!    & write(*,*)"anning_m1",tmp

                     qrtot    = qrtot + dum*qstot
                     nrtot    = nrtot + dum*nstot
                     qstot    = (one-dum)*qstot
                     nstot    = (one-dum)*nstot
                     preci(i) = (one-dum)*preci(i)
                   endif
                 endif

               end if
               tlataux(i, k) = tlat(i, k)

! freeze all rain at -5C for Arctic
               if(.true.) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 tx1 = t(i,k) + tlat(i,k)*onebcp*deltat
                 if (tx1 < tmelt-five) then

                   if (qrtot > zero) then

! make sure freezing rain doesn't increase temperature above threshold
                     dum = xlfocp*qrtot*rdzi
                     if (tx1+dum*deltat > tmelt-five) then
                       tx2 = -(tx1 - (tmelt-five)) * dti
                       dum = min(one,max(zero,tx2/dum))
                     else
                       dum = one
                     endif
                     qniic(i,k)  = qniic(i,k) + dum*qric(i,k)
                     nsic(i,k)   = nsic(i,k)  + dum*nric(i,k)
                     qric(i,k)   = (one-dum)*qric(i,k)
                     nric(i,k)   = (one-dum)*nric(i,k)
! heating tendency
                     tmp         = xlf*dum*qrtot*rdzi
                     frzrdt(i,k) = frzrdt(i,k) + tmp

                     tlat(i,k)   = tlat(i,k) + tmp

!      if(xlon<0.01.and.xlon>-0.01.and.xlat>1.346        
!    & .and.xlat<1.347.and.k==38)                        
!    & write(*,*)"anning_m2",tmp
  
                     qstot    = qstot + dum*qrtot
                     qrtot    = (one-dum)*qrtot
                     nstot    = nstot + dum*nrtot
                     nrtot    = (one-dum)*nrtot
                     preci(i) = preci(i) + dum*(prect(i)-preci(i))
                   end if
                 end if
               end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



! if rain/snow mix ratio is zero so should number concentration

               if (qniic(i,k) < qsmall) then
                 qniic(i,k) = zero 
                 nsic(i,k)  = zero
               end if

               if (qric(i,k) < qsmall) then
                 qric(i,k) =  zero
                 nric(i,k) =  zero
               end if

! make sure number concentration is a positive number to avoid
! taking root of negative

               nric(i,k) = max(nric(i,k), zero)
               nsic(i,k) = max(nsic(i,k), zero)

!.......................................................................
! get size distribution parameters for fallspeed calculations
!......................................................................
! rain

               if (qric(i,k) >= qsmall) then
                 lamr(k) = (pirhow*nric(i,k)/qric(i,k))**oneb3
                 n0r(k)  = nric(i,k)*lamr(k)

! check for slope
! change lammax and lammin for rain and snow
! adjust vars

                 if (lamr(k) < lamminr) then
                   lamr(k)   = lamminr
                   tx1       = lamr(k) * lamr(k)
                   n0r(k)    = tx1*tx1*qric(i,k)/pirhow
                   nric(i,k) = n0r(k)/lamr(k)
                 else if (lamr(k) > lammaxr) then
                   lamr(k)   = lammaxr
                   tx1       = lamr(k) * lamr(k)
                   n0r(k)    = tx1*tx1*qric(i,k)/pirhow
                   nric(i,k) = n0r(k)/lamr(k)
                 end if


! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

                 tx1    = arn(i,k) / lamr(k)**br
                 tx2    = 9.1_r8*rhof(i,k)
                 unr(k) = min(tx1*cons4, tx2)
                 umr(k) = min(tx1*(cons5/6._r8),tx2)
  
               else
                 lamr(k) = zero
                 n0r(k)  = zero
                 umr(k)  = zero
                 unr(k)  = zero
               end if

!calculate mean size of combined rain and snow

               if (lamr(k) > zero) then
                 Artmp = n0r(k) * (0.5*pi) / (lamr(k)*lamr(k)*lamr(k))
               else
                 Artmp = zero
               endif

               if (lamc(k) > zero) then
                 Actmp = cdist1(k) * pi * gamma(pgam(k)+three)
     &                 / (four * lamc(k)*lamc(k))
               else
                 Actmp = zero
               endif

               if (Actmp > zero .or.Artmp > zero) then
                 rercld(i,k) = rercld(i,k) + three*(qric(i,k)+qcic(i,k))
     &                                     /(four*rhow*(Actmp+Artmp))
                 arcld(i,k)  = arcld(i,k) + one
               endif

!......................................................................
! snow

               if (qniic(i,k) >= qsmall) then
                 lams(k) = (cons6*cs*nsic(i,k) / qniic(i,k))**(one/ds)
                 n0s(k)  = nsic(i,k)*lams(k)

! check for slope
! adjust vars

                 if (lams(k) < lammins) then
                   lams(k)   = lammins
                   n0s(k)    = lams(k)**(ds+one)*qniic(i,k)/(cs*cons6)
                   nsic(i,k) = n0s(k)/lams(k)

                 else if (lams(k) > lammaxs) then
                   lams(k)   = lammaxs
                   n0s(k)    = lams(k)**(ds+one)*qniic(i,k)/(cs*cons6)
                   nsic(i,k) = n0s(k)/lams(k)
                 end if

! 'final' values of number and mass weighted mean fallspeed for snow (m/s)

                 tx1    = asn(i,k) / lams(k)**bs
                 tx2    = 1.2_r8*rhof(i,k)
                 ums(k) = min(tx1*(cons8/6._r8), tx2)
                 uns(k) = min(tx1*cons7, tx2)

               else
                 lams(k) = zero
                 n0s(k)  = zero
                 ums(k)  = zero
                 uns(k)  = zero
               end if

!c........................................................................
! sum over sub-step for average process rates

! convert rain/snow q and N for output to history, note,
! output is for gridbox average

               qrout(i,k) = qrout(i,k) + qric(i,k)*cldmax(i,k)
               qsout(i,k) = qsout(i,k) + qniic(i,k)*cldmax(i,k)
               tx1        = rho(i,k)*cldmax(i,k)
               nrout(i,k) = nrout(i,k) + nric(i,k)*tx1
               nsout(i,k) = nsout(i,k) + nsic(i,k)*tx1

             end if !fprcp Anning Cheng 9/16/2016
               tlat1(i,k)     = tlat1(i,k)     + tlat(i,k)
               tlat1_aux(i,k) = tlat1_aux(i,k) + tlataux(i,k)

               qvlat1(i,k)  = qvlat1(i,k)  + qvlat(i,k)
               qctend1(i,k) = qctend1(i,k) + qctend(i,k)
               qitend1(i,k) = qitend1(i,k) + qitend(i,k)
               nctend1(i,k) = nctend1(i,k) + nctend(i,k)
               nitend1(i,k) = nitend1(i,k) + nitend(i,k)

               t(i,k)  = t(i,k)  + tlat(i,k)*deltat/cpp
               q(i,k)  = q(i,k)  + qvlat(i,k)*deltat
               qc(i,k) = qc(i,k) + qctend(i,k)*deltat
               qi(i,k) = qi(i,k) + qitend(i,k)*deltat
               nc(i,k) = nc(i,k) + nctend(i,k)*deltat
               ni(i,k) = ni(i,k) + nitend(i,k)*deltat

               rainrt1(i,k) = rainrt1(i,k) + rainrt(i,k)

!divide rain radius over substeps for average
               if (arcld(i,k) > zero) then
                 rercld(i,k) = rercld(i,k) / arcld(i,k)
               end if
!calculate precip fluxes and adding them to summing sub-stepping variables

               rflx(i,1) = zero
               sflx(i,1) = zero


               rflx(i,k+1) = qrout(i,k)*rho(i,k)*umr(k)
               sflx(i,k+1) = qsout(i,k)*rho(i,k)*ums(k)


               rflx1(i,k+1) = rflx1(i,k+1) + rflx(i,k+1)
               sflx1(i,k+1) = sflx1(i,k+1) + sflx(i,k+1)

!c........................................................................

             enddo    ! big k loop

             prect1(i) = prect1(i) + prect(i)
             preci1(i) = preci1(i) + preci(i)

!          if (lprint) write(0,*)' prect1=',prect1(i),' prect=',
!    &prect(i),' iter=',iter,' it=',it

           enddo !end of  big iter loop

           do k = 1, pver
             rate1ord_cw2pr_st(i,k) = qcsinksum_rate1ord(k)

     &                          / max(qcsum_rate1ord(k),1.0e-30_r8)
           enddo

         endif ! end of if (ltrue(i) == 0) then
       enddo   ! end of big i loop2

! convert dt from sub-step back to full time step

       deltat = deltat*real(iter)
       dti    = one / deltat


!.............................................................................

       do i=1,ncol  !big i loop3

         if (ltrue(i) == 0) then  ! skip all calculations if no cloud water

           do k=1,pver
                                  ! assign default values for effective radius
             effc(i,k)    = 10._r8
             effi(i,k)    = 25._r8
             effc_fn(i,k) = 10._r8
             lamcrad(i,k) = zero
             pgamrad(i,k) = zero
             deffi(i,k)   = zero
           end do

         else
!
           nstep = 1 ! initialize nstep for sedimentation sub-steps

! divide precip rate by number of sub-steps to get average over time step

           prect(i) = prect1(i) * riter
           preci(i) = preci1(i) * riter

!     if (lprint) write(0,*)' prect=',prect(i),' prect1=',prect1(i)
!    &,' riter=',riter

           do k=1,pver

! assign variables back to start-of-timestep values before updating after sub-steps

             t(i,k)  = t1(i,k)
             q(i,k)  = q1(i,k)
             qc(i,k) = qc1(i,k)
             qi(i,k) = qi1(i,k)
             nc(i,k) = nc1(i,k)
             ni(i,k) = ni1(i,k)

! divide microphysical tendencies by number of sub-steps to get average over time step

             tlat(i,k)    = tlat1(i,k)   * riter
             qvlat(i,k)   = qvlat1(i,k)  * riter
             qctend(i,k)  = qctend1(i,k) * riter
             qitend(i,k)  = qitend1(i,k) * riter
             nctend(i,k)  = nctend1(i,k) * riter
             nitend(i,k)  = nitend1(i,k) * riter

             tlataux(i,k) = tlat1_aux(i,k) * riter


             rainrt(i,k)  = rainrt1(i,k) * riter

! divide by number of sub-steps to find final values
             rflx(i,k+1)  = rflx1(i,k+1) * riter
             sflx(i,k+1)  = sflx1(i,k+1) * riter

! divide output precip q and N by number of sub-steps to get average over time step

             qrout(i,k)   = qrout(i,k) * riter
             qsout(i,k)   = qsout(i,k) * riter
             nrout(i,k)   = nrout(i,k) * riter
             nsout(i,k)   = nsout(i,k) * riter

! divide trop_mozart variables by number of sub-steps to get average over time step

             nevapr(i,k)   = nevapr(i,k)   * riter
             evapsnow(i,k) = evapsnow(i,k) * riter
             prain(i,k)    = prain(i,k)    * riter
             prodsnow(i,k) = prodsnow(i,k) * riter
             cmeout(i,k)   = cmeout(i,k)   * riter

             cmeiout(i,k) = cmeiout(i,k) * riter
             meltsdt(i,k) = meltsdt(i,k) * riter
             frzrdt (i,k) = frzrdt (i,k) * riter


! microphysics output
             prao(i,k)    = prao(i,k)    * riter
             prco(i,k)    = prco(i,k)    * riter
             mnuccco(i,k) = mnuccco(i,k) * riter
             mnuccto(i,k) = mnuccto(i,k) * riter
             msacwio(i,k) = msacwio(i,k) * riter
             psacwso(i,k) = psacwso(i,k) * riter
             bergso(i,k)  = bergso(i,k)  * riter
             bergo(i,k)   = bergo(i,k)   * riter
             prcio(i,k)   = prcio(i,k)   * riter
             praio(i,k)   = praio(i,k)   * riter

             mnuccdo(i,k) = mnuccdo(i,k) * riter
             mnuccto(i,k) = mnuccto(i,k) * riter

             mnuccro(i,k) = mnuccro(i,k) * riter
             pracso (i,k) = pracso (i,k) * riter

!!!!DONIFFFF==========================
             nnuccdo(i,k) = nnuccdo(i,k) * riter
             nnuccco(i,k) = nnuccco(i,k) * riter
             nsacwio(i,k) = nsacwio(i,k) * riter
             nsubio(i,k)  = nsubio(i,k)  * riter
             nprcio(i,k)  = nprcio(i,k)  * riter
             npraio(i,k)  = npraio(i,k)  * riter

             npccno(i,k)   = npccno(i,k)   * riter
             npsacwso(i,k) = npsacwso(i,k) * riter
             nsubco(i,k)   = nsubco(i,k)   * riter
             nprao(i,k)    = nprao(i,k)    * riter
             nprc1o(i,k)   = nprc1o(i,k)   * riter

!=====================================


! modify to include snow. in prain & evap (diagnostic here: for wet dep)

             prain(i,k) = prain(i,k) + prodsnow(i,k)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate sedimentation for cloud water and ice
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! update in-cloud cloud mixing ratio and number concentration
! with microphysical tendencies to calculate sedimentation, assign to dummy vars
! note: these are in-cloud values***, hence we divide by cloud fraction

             tx1 = one / lcldm(i,k)
             tx2 = one / icldm(i,k)
             dumc(i,k)  = (qc(i,k) + qctend(i,k)*deltat) * tx1
             dumi(i,k)  = (qi(i,k) + qitend(i,k)*deltat) * tx2
             dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)*tx1, zero)
             dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)*tx2, zero)

             if (nccons) then
               dumnc(i,k) = ncnst*irho(i,k)
             end if

             if (nicons) then
               dumni(i,k) = ninst*irho(i,k)
             end if

! obtain new slope parameter to avoid possible singularity

             if (dumi(i,k) >= qsmall) then
! add upper limit to in-cloud number concentration to prevent numerical error
               dumni(i,k) = min(dumni(i,k),dumi(i,k)*1.e20_r8)
               lami(k)    = (cons1*ci* dumni(i,k)/dumi(i,k))**oneodi

!      miu_ice(k)=mui_hemp_l(lami(k)) Anning changed here
               miu_ice(k) = max(min(0.008_r8*(lami(k)*0.01)**0.87_r8,
     &                                             10.0_r8), 0.1_r8)
               tx1        = one + miu_ice(k)
               tx2        = one / gamma(tx1)
               aux        = (gamma(tx1+di)*tx2) ** oneodi
               lami(k)    = aux*lami(k)

               n0i(k)     = niic(i,k)*lami(k)**tx1 * tx2


               if (lami(k) < lammini*aux) then
                 lami(k)    = lammini
                 miu_ice(k) = zero
                 niic(i,k)  = n0i(k)/lami(k)
               end if
               if (lami(k) > lammaxi*aux) then
                 lami(k)    = lammaxi
                 miu_ice(k) = zero
                 niic(i,k)  = n0i(k)/lami(k)
               end if

             else
               lami(k) = zero
             end if

             if (dumc(i,k) >= qsmall) then
! add upper limit to in-cloud number concentration to prevent numerical error
               dumnc(i,k) = min(dumnc(i,k),dumc(i,k)*1.e20_r8)
! add lower limit to in-cloud number concentration
               dumnc(i,k) = max(dumnc(i,k),cdnl*irho(i,k))
               tx1 = 0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
               pgam(k)    = max(two, min(15._r8, one/(tx1*tx1)-one))

               lamc(k) = (cr*dumnc(i,k)*gamma(pgam(k)+four)
     &                 / (dumc(i,k)*gamma(pgam(k)+one)))**oneb3
               lammin  = (pgam(k)+one) / 50.e-6_r8
               lammax  = (pgam(k)+one) / 2.e-6_r8
               lamc(k) = min(lammax, max(lamc(k),lammin))
             else
               lamc(k) = zero
             end if

! calculate number and mass weighted fall velocity for droplets
! include effects of sub-grid distribution of cloud water


             if (dumc(i,k) >= qsmall) then
               tx1 = lamc(k) ** bc
               unc = acn(i,k)*gamma(one+bc+pgam(k))
     &             / (tx1*gamma(pgam(k)+one))
               umc = acn(i,k)*gamma(four+bc+pgam(k))
     &             / (tx1*gamma(pgam(k)+four))
! fallspeed for output
               vtrmc(i,k) = umc
             else
               umc = zero
               unc = zero
             end if



! calculate number and mass weighted fall velocity for cloud ice

             if (dumi(i,k) >= qsmall) then
               cons16 = gamma(one+bi+miu_ice(k))/gamma(one+miu_ice(k))
               cons17 = gamma(four+bi+miu_ice(k))/gamma(four+miu_ice(k))

               tx1    = ain(i,k) / lami(k)**bi
               tx2    = 1.2_r8*rhof(i,k)
               uni    = min(tx1*cons16, tx2)
               umi    = min(tx1*cons17, tx2)
! fallspeed
               vtrmi(i,k) = umi
             else
               umi = zero
               uni = zero
             end if

!DONIFFFF tune up sedimentation vel!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             umi = umi*ui_scale
             uni = uni*ui_scale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             tx1    = g*rho(i,k)
             fi(k)  = tx1*umi
             fni(k) = tx1*uni
             fc(k)  = tx1*umc
             fnc(k) = tx1*unc

! calculate number of split time steps to ensure courant stability criteria
! for sedimentation calculations

             rgvm  = max(fi(k),fc(k),fni(k),fnc(k))
             nstep = max(int(rgvm*deltat/pdel(i,k)+one),nstep)

! redefine dummy variables - sedimentation is calculated over grid-scale
! quantities to ensure conservation

             dumc(i,k)  = qc(i,k) + qctend(i,k)*deltat
             dumi(i,k)  = qi(i,k) + qitend(i,k)*deltat
             dumnc(i,k) = max(nc(i,k) + nctend(i,k)*deltat, zero)
             dumni(i,k) = max(ni(i,k) + nitend(i,k)*deltat, zero)

             if (dumc(i,k) < qsmall) dumnc(i,k) = zero
             if (dumi(i,k) < qsmall) dumni(i,k) = zero

           end do         ! end of k loop

           do n = 1,nstep

             do k = 1,pver
               falouti(k)  = max(fi(k)  * dumi(i,k),  zero)
               faloutni(k) = max(fni(k) * dumni(i,k), zero)
               faloutc(k)  = max(fc(k)  * dumc(i,k),  zero)
               faloutnc(k) = max(fnc(k) * dumnc(i,k), zero)
             end do 

! top of model

             k = 1
             tx1      = one / pdel(i,k)
             faltndi  = falouti(k)  * tx1
             faltndni = faloutni(k) * tx1
             faltndc  = faloutc(k)  * tx1
             faltndnc = faloutnc(k) * tx1

! add fallout terms to microphysical tendencies

             tx2     = one / float(nstep)
             tx3     = deltat*tx2

             qitend(i,k) = qitend(i,k) - faltndi  * tx2
             nitend(i,k) = nitend(i,k) - faltndni * tx2
             qctend(i,k) = qctend(i,k) - faltndc  * tx2
             nctend(i,k) = nctend(i,k) - faltndnc * tx2

! sedimentation tendencies for output
             qcsedten(i,k) = qcsedten(i,k) - faltndc * tx2
             qisedten(i,k) = qisedten(i,k) - faltndi * tx2

             dumi(i,k)  = dumi(i,k)  - faltndi *tx3
             dumni(i,k) = dumni(i,k) - faltndni*tx3
             dumc(i,k)  = dumc(i,k)  - faltndc *tx3
             dumnc(i,k) = dumnc(i,k) - faltndnc*tx3

             do k = 2,pver

! for cloud liquid and ice, if cloud fraction increases with height
! then add flux from above to both vapor and cloud water of current level
! this means that flux entering clear portion of cell from above evaporates
! instantly
               tx4      = tx1
               tx1      = one / pdel(i,k)

               if (lcldm(i,k-1) > zero) then
                 dum  = min(one, lcldm(i,k)/lcldm(i,k-1))
               else
                 dum = min(one, lcldm(i,k))
               endif
               if (icldm(i,k-1) > zero) then
                 dum1 = min(one, icldm(i,k)/icldm(i,k-1))
               else
                 dum1 = min(one, icldm(i,k))
               endif

               faltndqie = (falouti(k)  - falouti(k-1))       * tx1
               faltndi   = (falouti(k)  - dum1*falouti(k-1))  * tx1
               faltndni  = (faloutni(k) - dum1*faloutni(k-1)) * tx1

               faltndqce = (faloutc(k) -     faloutc(k-1))    * tx1
               faltndc   = (faloutc(k) - dum*faloutc(k-1))    * tx1
               faltndnc  = (faloutnc(k)- dum*faloutnc(k-1))   * tx1

! add fallout terms to eulerian tendencies

               qitend(i,k) = qitend(i,k) - faltndi  * tx2
               nitend(i,k) = nitend(i,k) - faltndni * tx2
               qctend(i,k) = qctend(i,k) - faltndc  * tx2
               nctend(i,k) = nctend(i,k) - faltndnc * tx2


! sedimentation tendencies for output
               qcsedten(i,k) = qcsedten(i,k) - faltndc * tx2
               qisedten(i,k) = qisedten(i,k) - faltndi * tx2

! add terms to to evap/sub of cloud water

               qvlat(i,k)   = qvlat(i,k)   - (faltndqie-faltndi) * tx2

!     if (lprint .and. k == 29) write(0,*)' qvlatb=',qvlat(i,k),
!    &' tx2=',tx2,' faltndqie=',faltndqie,' faltndi=',faltndi
! for output
               qisevap(i,k) = qisevap(i,k) + (faltndqie-faltndi) * tx2


               qvlat(i,k)   = qvlat(i,k)   - (faltndqce-faltndc) * tx2
!     if (lprint .and. k == 29) write(0,*)' qvlatc=',qvlat(i,k),
!    &' tx2=',tx2,' faltndqce=',faltndqce,' faltndc=',faltndc

               qcsevap(i,k) = qcsevap(i,k) + (faltndqce-faltndc) * tx2


               tlat(i,k) = tlat(i,k) + ((faltndqie-faltndi)*xxls
     &                               +  (faltndqce-faltndc)*xxlv) * tx2


               dumi(i,k)  = dumi(i,k)  - faltndi  * tx3
               dumni(i,k) = dumni(i,k) - faltndni * tx3
               dumc(i,k)  = dumc(i,k)  - faltndc  * tx3
               dumnc(i,k) = dumnc(i,k) - faltndnc * tx3

               Fni(K) = MAX(Fni(K)*tx1, Fni(K-1)*tx4) * pdel(i,K)
               FI(K)  = MAX(FI(K)*tx1,  FI(K-1)*tx4)  * pdel(i,K)
               fnc(k) = max(fnc(k)*tx1, fnc(k-1)*tx4) * pdel(i,k)
               Fc(K)  = MAX(Fc(K)*tx1,  Fc(K-1)*tx4)  * pdel(i,K)

             end do            ! end of k loop

! units below are m/s
! cloud water/ice sedimentation flux at surface
! is added to precip flux at surface to get total precip (cloud + precip water)
! rate

!     if (lprint) write(0,*)' befend prect=',prect(i),' preci=',preci(i)
             tx3      = tx2 / (g*1000._r8)
             prect(i) = prect(i) + (faloutc(pver)+falouti(pver)) * tx3
             preci(i) = preci(i) + falouti(pver) * tx3

!     if (lprint) write(0,*)' end prect=',prect(i),' preci=',preci(i)

           end do        ! end of n loop


! end sedimentation
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! calculate sedimentation for rain and snow
! Anning Cheng 9/19/2016, forecast rain and snow 
! reuse dummy variable for cloud water and ice
!  iter =1 for fprcp >= 1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! initialize nstep for sedimentation sub-steps
! reuse dumc, dumi, dumnc, and dumni
         if(fprcp==1) then
           nstep = 1

           do k=1,pver

!            tx1 = one / lcldm(i,k)
!            tx2 = one / icldm(i,k)
             dumc(i,k)  = max((qrn(i,k)+qrtend(i,k)*deltat),zero)
             dumi(i,k)  = max((qsnw(i,k)+qnitend(i,k)*deltat),zero)
             dumnc(i,k) = max((nrn(i,k)+nrtend(i,k)*deltat),zero)
             dumni(i,k) = max((nsnw(i,k)+nstend(i,k)*deltat),zero)

! if rain/snow mix ratio is zero so should number concentration

               if (dumi(i,k) < qsmall) then
                 dumi(i,k)  = zero
                 dumni(i,k)  = zero
               endif

               if (dumc(i,k) < qsmall) then
                 dumc(i,k)  = zero
                 dumnc(i,k) = zero
               endif

! make sure number concentration is a positive number to avoid
! taking root of negative

               dumnc(i,k) = max(dumnc(i,k),zero)
               dumni(i,k) = max(dumni(i,k),zero)

!.......................................................................
! get size distribution parameters for fallspeed calculations
!......................................................................
! rain

               if (dumc(i,k) >= qsmall) then
                 lamr(k) = (pi*rhow*dumnc(i,k)/dumc(i,k))**oneb3
                 n0r(k)  = dumnc(i,k)*lamr(k)

! check for slope
! change lammax and lammin for rain and snow
! adjust vars

                 if (lamr(k) < lamminr) then

                   lamr(k)   = lamminr

                   n0r(k)    = lamr(k)**4*dumc(i,k)/(pi*rhow)
                   dumnc(i,k) = n0r(k)/lamr(k)
                 else if (lamr(k) > lammaxr) then
                   lamr(k)   = lammaxr
                   n0r(k)    = lamr(k)**4*dumc(i,k)/(pi*rhow)
                   dumnc(i,k) = n0r(k)/lamr(k)
                 end if


! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

                 tx1    = arn(i,k) / lamr(k)**br
                 tx2    = 9.1_r8*rhof(i,k)
                 unr(k) = min(tx1*cons4, tx2)
                 umr(k) = min(tx1*(cons5/6._r8),tx2)

               else
                 lamr(k) = zero
                 n0r(k)  = zero
                 umr(k)  = zero
                 unr(k)  = zero
               end if

!......................................................................
! snow

               if (dumi(i,k) >= qsmall) then
                 lams(k) = (cons6*cs*dumni(i,k)/ dumi(i,k))**(one/ds)
                 n0s(k)  = dumni(i,k)*lams(k)

! check for slope
! adjust vars

                 if (lams(k) < lammins) then
                   lams(k)   = lammins
                   n0s(k)    = lams(k)**(ds+one)*dumi(i,k)/(cs*cons6)
                   dumni(i,k) = n0s(k)/lams(k)

                 else if (lams(k) > lammaxs) then
                   lams(k)   = lammaxs
                   n0s(k)    = lams(k)**(ds+one)*dumi(i,k)/(cs*cons6)
                   dumni(i,k) = n0s(k)/lams(k)
                 end if

! 'final' values of number and mass weighted mean fallspeed for snow (m/s)

                 tx1    = asn(i,k) / lams(k)**bs
                 tx2    = 1.2_r8*rhof(i,k)
                 ums(k) = min(tx1*(cons8/6._r8), tx2)
                 uns(k) = min(tx1*cons7, tx2)

               else
                 lams(k) = zero
                 n0s(k)  = zero
                 ums(k)  = zero
                 uns(k)  = zero
               end if

             tx1   = g*rho(i,k)
             fi(k)  = tx1*ums(k)
             fni(k) = tx1*uns(k)
             fc(k)  = tx1*umr(k)
             fnc(k) = tx1*unr(k)

! calculate number of split time steps to ensure courant stability criteria
! for sedimentation calculations

             rgvm = max(fi(k),fc(k),fni(k),fnc(k))
             nstep = max(int(rgvm*deltat/pdel(i,k)+one),nstep)

! redefine dummy variables - sedimentation is calculated over grid-scale
! quantities to ensure conservation

             qrn(i,k)  = (qrn(i,k)  + qrtend(i,k)*deltat)
             qsnw(i,k) = (qsnw(i,k) + qnitend(i,k)*deltat)
             nrn(i,k)  = max((nrn(i,k)  + nrtend(i,k)*deltat),zero)
             nsnw(i,k) = max((nsnw(i,k) + nstend(i,k)*deltat),zero)

             if (qrn(i,k)  < qsmall) nrn(i,k)  = zero
             if (qsnw(i,k) < qsmall) nsnw(i,k) = zero

           enddo         ! end of k loop

           tx2     = one / float(nstep)
           tx3     = deltat*tx2
           do n = 1,nstep

             do k = 1,pver
               falouti(k)  = max(fi(k)  * qsnw(i,k),zero)
               faloutni(k) = max(fni(k) * nsnw(i,k),zero)
               faloutc(k)  = max(fc(k)  * qrn(i,k),zero)
               faloutnc(k) = max(fnc(k) * nrn(i,k),zero)
             end do 

! top of model

             k = 1
             tx1      = one / pdel(i,k)
             faltndi  = falouti(k)  * tx1
             faltndni = faloutni(k) * tx1
             faltndc  = faloutc(k)  * tx1
             faltndnc = faloutnc(k) * tx1

! add fallout terms to microphysical tendencies

!            qnitend(i,k) = qnitend(i,k) - faltndi  * tx2
!            nstend(i,k) = nstend(i,k) - faltndni * tx2
!            qrtend(i,k) = qrtend(i,k) - faltndc  * tx2
!            nrtend(i,k) = nrtend(i,k) - faltndnc * tx2

! sedimentation tendencies for output

             qsnw(i,k)  = qsnw(i,k)  - faltndi *tx3
             nsnw(i,k) = nsnw(i,k) - faltndni*tx3
             qrn(i,k)  = qrn(i,k)  - faltndc *tx3
             nrn(i,k) = nrn(i,k) - faltndnc*tx3

             do k = 2,pver

! for rain and snow
               tx1      = one / pdel(i,k)
!              dum  = min(one, lcldm(i,k)/lcldm(i,k-1))
!              dum1 = min(one, icldm(i,k)/icldm(i,k-1))

               faltndc   = (faloutc(k)- faloutc(k-1)) * tx1
               faltndnc  = (faloutnc(k)- faloutnc(k-1)) * tx1
   
               faltndi   = (falouti(k) - falouti(k-1)) * tx1
               faltndni  = (faloutni(k) - faloutni(k-1)) * tx1

               qsnw(i,k)  = qsnw(i,k)  - faltndi  * tx3
               nsnw(i,k) = nsnw(i,k) - faltndni * tx3
               qrn(i,k)  = qrn(i,k)  - faltndc  * tx3
               nrn(i,k) = nrn(i,k) - faltndnc * tx3

             end do


             tx5      = tx2 / (g*1000._r8)
             prect(i) = prect(i) + (faloutc(pver)+falouti(pver)) * tx5
             preci(i) = preci(i) + (falouti(pver)) * tx5

           end do        ! end of n loop
         end if !fprcp ==1
! end sedimentation for rain and snow
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! get new update for variables that includes sedimentation tendency
! note : here dum variables are grid-average, NOT in-cloud

! DONE UP TO HERE
           do k=1,pver
             if (fprcp==1) then
! calculate melting and freezing of precip
! melt snow at +2 C

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               tx1 = t(i,k) + tlat(i,k)*onebcp*deltat
               if (tx1 > 275.15_r8) then
                 if (qsnw(i,k) > zero) then

! make sure melting snow doesn't reduce temperature below threshold
                   dum = xlfocp*qsnw(i,k)

                   if (tx1-dum < 273.15_r8)  then
                     tx2 = tx1-275.15_r8
                     dum = min(one,max(zero,tx2/dum))
                   else
                     dum = one
                   end if

                   qrn(i,k)  = qrn(i,k) + dum*qsnw(i,k)
                   nrn(i,k)  = nrn(i,k) + dum*nsnw(i,k)
                   qsnw(i,k) = (one-dum)*qsnw(i,k)
                   nsnw(i,k)  = (one-dum)*nsnw(i,k)
! heating tendency
                   tmp = -xlf*dum*qsnw(i,k)/deltat

                   tlat(i,k) = tlat(i,k) + tmp

                   preci(i) = (one-dum)*preci(i)
                 end if
               end if

! freeze all rain at -5C for Arctic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               tx1 = t(i,k) + tlat(i,k)*onebcp*deltat
               if (tx1 < tmelt-five) then

                 if (qrn(i,k) > zero) then

! make sure freezing rain doesn't increase temperature above threshold
                   dum = xlfocp*qrn(i,k)
                   if (tx1+dum > tmelt-five) then
                     tx2 = -(tx1-(tmelt-5._r8))
                     dum = min(one,max(zero,tx2/dum))
                   else
                     dum = one
                   end if

                   qsnw(i,k) = qsnw(i,k) + dum*qrn(i,k)
                   nsnw(i,k) = nsnw(i,k)  + dum*nrn(i,k)
                   qrn(i,k)  = (one-dum)*qrn(i,k)
                   nrn(i,k)  = (one-dum)*nrn(i,k)
! heating tendency
                   tmp       = xlf*dum*qrn(i,k)/deltat

                   tlat(i,k) = tlat(i,k) + tmp

                   preci(i)  = preci(i) + dum*(prect(i)-preci(i))
                 end if
               end if

! if rain/snow mix ratio is zero so should number concentration

               if (qsnw(i,k) < qsmall) then
                 nsnw(i,k)  = zero
               end if

               if (qrn(i,k) < qsmall) then
                 nrn(i,k) = zero
               end if

! make sure number concentration is a positive number to avoid
! taking root of negative

               nrn(i,k)  = max(nrn(i,k),zero)
               nsnw(i,k) = max(nsnw(i,k),zero)

!.......................................................................
             end if !fprcp ==1



             dumc(i,k)  = max(qc(i,k)+qctend(i,k)*deltat, zero)
             dumi(i,k)  = max(qi(i,k)+qitend(i,k)*deltat, zero)
             dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat, zero)
             dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat, zero)


             if (nccons) then
               dumnc(i,k) = ncnst * irho(i,k) * lcldm(i,k)
             endif

             if (nicons) then
               dumni(i,k) = ninst * irho(i,k) * icldm(i,k)
             endif

             if (dumc(i,k) < qsmall) dumnc(i,k) = zero
             if (dumi(i,k) < qsmall) dumni(i,k) = zero

! calculate instantaneous processes (melting, homogeneous freezing)

             tx1 = t(i,k) + tlat(i,k)*onebcp*deltat
             if (tx1 > tmelt) then
               if (dumi(i,k) > zero) then

! limit so that melting does not push temperature below freezing
                 dum = dumi(i,k)*xlfocp
                 if (tx1-dum < tmelt) then
                   dum = (tx1-tmelt)*cpoxlf / dumi(i,k)
                   dum = max(zero, min(one, dum))
                 else
                   dum = one
                 end if

                 tx2 = dum*dumi(i,k)*dti
                 qctend(i,k) = qctend(i,k) + tx2
! for output
                 melto(i,k)  = tx2

! assume melting ice produces droplet
! mean volume radius of 8 micron

                 nctend(i,k) = nctend(i,k) + three*tx2
     &                                     / (four*pi*5.12e-16_r8*rhow)

                 qitend(i,k) = ((one-dum)*dumi(i,k)-qi(i,k))  * dti
                 nitend(i,k) = ((one-dum)*dumni(i,k)-ni(i,k)) * dti
                 tlat(i,k)   = tlat(i,k) - xlf*tx2

               endif
             endif

! homogeneously freeze droplets between  -35 C and -40 C

             tx1 = t(i,k) + tlat(i,k)*onebcp*deltat
             if (tx1 < 233.15_r8) then

               if (dumc(i,k) > zero .and. qc(i,k) > zero) then

! limit so that freezing does not push temperature above threshold
                 dum = dumc(i,k)*xlfocp
                 if (tx1+dum > 233.15_r8) then
                   dum = -(tx1-233.15_r8)*cpoxlf / dumc(i,k)
                   dum = max(zero, min(one, dum))
                 else
                   dum = one
                 end if

                 tx2         = dum*dumc(i,k)*dti
                 qitend(i,k) = qitend(i,k) + tx2
! for output
                 homoo(i,k)  = tx2

! assume 25 micron mean volume radius of homogeneously frozen droplets
! consistent with size of detrained ice in stratiform.F90

                 nitend(i,k) = nitend(i,k) + dum*three*dumc(i,k)
     &                     / (four*3.14_r8*1.563e-14_r8* 500._r8) * dti
                 qctend(i,k) = ((one-dum)*dumc(i,k)-qc(i,k))  * dti
                 nctend(i,k) = ((one-dum)*dumnc(i,k)-nc(i,k)) * dti
                 tlat(i,k)   = tlat(i,k) + xlf*tx2

               endif
             endif

! remove any excess over-saturation, which is possible due to non-linearity when adding
! together all microphysical processes
! follow code similar to old CAM scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      qsn = min(epsqs*esn/(p(i,k)-omeps*esn),one)
!      if (qtmp > qsn .and. qsn > 0) then
! expression below is approximate since there may be ice deposition
!      dum = (qtmp-qsn)/(one+cons27*qsn/(cpp*rv*ttmp**2))/deltat
! add to output cme
! now add to tendencies, partition between liquid and ice based on temperature
! now add to tendencies, partition between liquid and ice based on te
!      dum = (qtmp-qsn)/(one+(xxls*dum1+xxlv*(one-dum1))**2 &
!                     *qsn/(cpp*rv*ttmp**2))/deltat
!
!      qctend(i,k)=qctend(i,k)+dum*(one-dum1)
! for output
!        qcreso(i,k)=dum*(one-dum1)
!      qitend(i,k)=qitend(i,k)+dum*dum1
!      qireso(i,k)=dum*dum1
!      qvlat(i,k)=qvlat(i,k)-dum
! for output
!        qvres(i,k)=-dum
!      tlat(i,k)=tlat(i,k)+dum*(one-dum1)*xxlv+dum*dum1*xxls
!      end if
!end if

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!1


!...............................................................................
! calculate effective radius for pass to radiation code
! if no cloud water, default value is 10 micron for droplets,
! 25 micron for cloud ice

! update cloud variables after instantaneous processes to get effective radius
! variables are in-cloud to calculate size dist parameters

             tx1 = one / lcldm(i,k)
             tx2 = one / icldm(i,k)
             dumc(i,k)  = max(qc(i,k)+qctend(i,k)*deltat, zero) * tx1
             dumi(i,k)  = max(qi(i,k)+qitend(i,k)*deltat, zero) * tx2
             dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat, zero) * tx1
             dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat, zero) * tx2

!      if ((dumc(i, k)*1e6 .gt. 1.0) .and. dumnc(i, k) .lt. 1e-20) then!

!        print *, 'dumnc', dumnc(i,k)*1e-6
!        print *, 'dumc',  dumc(i, k)*1e6
!        print *, i, k

!      end if

             if (nccons) then
               dumnc(i,k) = ncnst * irho(i,k)
             end if

             if (nicons) then
               dumni(i,k) = ninst * irho(i,k)
             end if


! limit in-cloud mixing ratio of water and ice to reasonable value of 5 g kg-1

             dumc(i,k) = min(dumc(i,k),5.e-3_r8)
             dumi(i,k) = min(dumi(i,k),5.e-3_r8)
 
!...................
! cloud ice effective radius

             if (dumi(i,k) >= qsmall) then
! add upper limit to in-cloud number concentration to prevent numerical error
               dumni(i,k) = min(dumni(i,k),dumi(i,k)*1.e20_r8)
               lami(k)    = (cons1*ci* dumni(i,k)/dumi(i,k))**oneodi

!              miu_ice(k) = mui_hemp_l(lami(k))
               miu_ice(k) = max(min(0.008_r8*(lami(k)*0.01)**0.87_r8,
     &                                            10.0_r8), 0.1_r8)
               tx1        = one + miu_ice(k)
               tx2        = one / gamma(tx1)
               aux        = (gamma(tx1+di)*tx2) ** oneodi
               lami(k)    = aux*lami(k)
               n0i(k)     = niic(i,k) * lami(k)**tx1 * tx2

               if (lami(k) < lammini*aux) then
                 miu_ice(k)  = zero
                 lami(k)     = lammini
                 n0i(k)      = lami(k)**(di+one)*dumi(i,k)/(ci*cons1)
                 niic(i,k)   = n0i(k)/lami(k)
! adjust number conc if needed to keep mean size in reasonable range
                 nitend(i,k) = (niic(i,k)*icldm(i,k)-ni(i,k))/deltat
               else if (lami(k) > lammaxi*aux) then
                 miu_ice(k)  = zero
                 lami(k)     = lammaxi
                 n0i(k)      = lami(k)**(di+one)*dumi(i,k)/(ci*cons1)
                 niic(i,k)   = n0i(k)/lami(k)
                 nitend(i,k) = (niic(i,k)*icldm(i,k)-ni(i,k))/deltat
               end if

               effi(i,k) = 1.e6_r8*gamma(four+miu_ice(k))
     &                   / (two*lami(k)*gamma(three+miu_ice(k)))
             else
               effi(i,k) = 25._r8
             end if

!...................
! cloud droplet effective radius


             if (dumc(i,k) >= qsmall) then

! add upper limit to in-cloud number concentration to prevent numerical error
               dumnc(i,k) = min(dumnc(i,k),dumc(i,k)*1.e20_r8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set tendency to ensure minimum droplet concentration
! after update by microphysics, except when lambda exceeds bounds on mean drop
! size or if there is no cloud water
               if (dumnc(i,k) < cdnl*irho(i,k)) then
                 nctend(i,k) = (cdnl*irho(i,k)*cldm(i,k)-nc(i,k))*dti
               end if
               dumnc(i,k) = max(dumnc(i,k),cdnl*irho(i,k))

               if (nccons) then
! make sure nc is consistence with the constant N by adjusting tendency, need
! to multiply by cloud fraction
! note that nctend may be further adjusted below if mean droplet size is
! out of bounds

                 nctend(i,k) = (ncnst*irho(i,k)*lcldm(i,k)-nc(i,k))
     &                       * dti
               end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               tx1 = 0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
               pgam(k) = max(two, min(15._r8, one/(tx1*tx1)-one))

               tx1     = gamma(pgam(k)+one)
               tx2     = gamma(pgam(k)+four)
               lamc(k) = (cr*dumnc(i,k)*tx2
     &                 / (dumc(i,k)*tx1))**oneb3
               lammin  = (pgam(k)+one) / 50.e-6_r8
               lammax  = (pgam(k)+one) / 2.e-6_r8
               if (lamc(k) < lammin) then
                 lamc(k)   = lammin
                 ncic(i,k) = 6._r8*lamc(k)*lamc(k)*lamc(k)*dumc(i,k)
     &                     * tx1 / (pirhow*tx2)
! adjust number conc if needed to keep mean size in reasonable range
                nctend(i,k) = (ncic(i,k)*lcldm(i,k)-nc(i,k))*dti

               else if (lamc(k) > lammax) then
                 lamc(k)   = lammax
                 ncic(i,k) = 6._r8*lamc(k)*lamc(k)*lamc(k)*dumc(i,k)
     &                     * tx1 / (pirhow*tx2)
! adjust number conc if needed to keep mean size in reasonable range
                 nctend(i,k) = (ncic(i,k)*lcldm(i,k)-nc(i,k))*dti
               end if

               effc(i,k) = tx2 / (gamma(pgam(k)+three)*lamc(k)*2.e6_r8)
!assign output fields for shape here
               lamcrad(i,k) = lamc(k)
               pgamrad(i,k) = pgam(k)
 
             else
               effc(i,k)    = 10._r8
               lamcrad(i,k) = zero
               pgamrad(i,k) = zero
             end if
 

! ice effective diameter for david mitchell's optics
             deffi(i,k) = effi(i,k) * (rhoi / 917._r8*two)


!!! recalculate effective radius for constant number, in order to separate
! first and second indirect effects
! assume constant number of 10^8 kg-1

             dumnc(i,k) = 1.e8

             if (dumc(i,k) >= qsmall) then
               tx1 = 0.0005714_r8*(ncic(i,k)/1.e6_r8*rho(i,k))+0.2714_r8
               pgam(k) = max(two, min(15._r8, one/(tx1*tx1)-one))
               tx2     = gamma(pgam(k)+four)
               lamc(k) = (cr*dumnc(i,k)*tx2
     &                 / (dumc(i,k)*gamma(pgam(k)+one)))**oneb3
               lammin  = (pgam(k)+one) / 50.e-6_r8
               lammax  = (pgam(k)+one) / 2.e-6_r8
               if (lamc(k) < lammin) then
                 lamc(k) = lammin
               else if (lamc(k) > lammax) then
                 lamc(k) = lammax
               end if
               effc_fn(i,k) = tx2/(gamma(pgam(k)+three)*lamc(k)*2.e6_r8)

             else
               effc_fn(i,k) = 10._r8
             end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

           end do       ! end of k loop after n loop


         endif ! end of if (ltrue(i) == 0) loop

! convert dt from sub-step back to full time step

       deltat = deltat*real(iter)
       dti    = one / deltat


       do k=1,pver
! if updated q (after microphysics) is zero, then ensure updated n is also zero

           if (qc(i,k)+qctend(i,k)*deltat < qsmall)
     &                 nctend(i,k) = -nc(i,k) * dti
           if (qi(i,k)+qitend(i,k)*deltat < qsmall)
     &                 nitend(i,k) = -ni(i,k) * dti
         end do


       end do  !end big i loop3


!print*, 'nctend', 1.0e-6*nctend*deltat

! hm add rain/snow mixing ratio and number concentration as diagnostic

#ifdef  CAM
       call outfld('QRAIN',qrout, pcols, lchnk)
       call outfld('QSNOW',qsout, pcols, lchnk)
       call outfld('NRAIN',nrout, pcols, lchnk)
       call outfld('NSNOW',nsout, pcols, lchnk)
#endif

! add snow output
       do k=1,pver
         do i = 1,ncol
           if (qsout(i,k) > 1.e-7_r8 .and. nsout(i,k) > zero) then
             dsout(i,k) = (pirhosn * nsout(i,k)/qsout(i,k)) **(-oneb3)
           endif
         end do
       end do

#ifdef CAM
       call outfld('DSNOW',dsout, pcols, lchnk)

! add ip fluxes as output fields
       call outfld('MGFLXPRC',rflx, pcols, lchnk)
       call outfld('MGFLXSNW',sflx, pcols, lchnk)
#endif

       do k=1,pver
         do i = 1,ncol
           if ((qc(i,k)+qctend(i,k)*deltat >= qsmall) .and.
     &         (cldmax(i,k) > zero) .and. (lcldm(i, k) > zero) .and.
     &         (nc(i,k)+nctend(i,k)*deltat > 0.0)) then

             tx1 = rho(i,k) / lcldm(i,k)
             tx2 = (qc(i,k)+qctend(i,k)*deltat)*tx1*1000._r8
             dum = tx2 * tx2 * lcldm(i,k)
     &           / (0.109_r8*(nc(i,k)+nctend(i,k)*deltat)
     &                      *tx1*1.e-6_r8*cldmax(i,k))
           else
             dum = zero
           end if
           if (qi(i,k)+qitend(i,k)*deltat >= qsmall) then
             dum1 = ((qi(i,k)+qitend(i,k)*deltat)*rho(i,k)/icldm(i,k)
     &            * 1000._r8/0.1_r8)**(one/0.63_r8)
     &            * icldm(i,k)/cldmax(i,k)
           else
             dum1 = zero
           end if

           if (qsout(i,k) >= qsmall) then
             dum1 = dum1 + (qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)
     &                                          **(one/0.63_r8)
           end if

           refl(i,k) = dum + dum1

           if (rainrt(i,k) >= 0.001) then
             dum = log10(rainrt(i,k)**6._r8) + 16._r8
             dum = 10._r8**(dum/10._r8)
           else
             dum = 0.
           end if

           refl(i,k)   = refl(i,k) + dum

           areflz(i,k) = refl(i,k)

           if (refl(i,k) > minrefl) then
             refl(i,k) = 10._r8*log10(refl(i,k))
           else
             refl(i,k) = -9999._r8
           end if

           if (refl(i,k) > mindbz) then
             arefl(i,k)  = refl(i,k)
             frefl(i,k)  = one
           else
             arefl(i,k)  = zero
             areflz(i,k) = zero
             frefl(i,k)  = zero
           end if

           csrfl(i,k) = min(csmax,refl(i,k))

           if (csrfl(i,k) > csmin) then
             acsrfl(i,k) = refl(i,k)
             fcsrfl(i,k) = one
           else
             acsrfl(i,k) = zero
             fcsrfl(i,k) = zero
           end if

         end do
       end do

#ifdef CAM
       call outfld('REFL',refl, pcols, lchnk)
       call outfld('AREFL',arefl, pcols, lchnk)
       call outfld('AREFLZ',areflz, pcols, lchnk)
       call outfld('FREFL',frefl, pcols, lchnk)
       call outfld('CSRFL',csrfl, pcols, lchnk)
       call outfld('ACSRFL',acsrfl, pcols, lchnk)
       call outfld('FCSRFL',fcsrfl, pcols, lchnk)

       call outfld('RERCLD',rercld, pcols, lchnk)
#endif

! averaging for snow and rain number and diameter

       do k=1,pver
         do i = 1,ncol
           qrout2(i,k) = zero
           nrout2(i,k) = zero
           drout2(i,k) = zero
           freqr(i,k)  = zero

           qsout2(i,k) = zero
           nsout2(i,k) = zero
           dsout2(i,k) = zero
           freqs(i,k)  = zero

           if (qrout(i,k) > 1.e-7_r8 .and. nrout(i,k) > zero) then
             qrout2(i,k) = qrout(i,k)
             nrout2(i,k) = nrout(i,k)
             drout2(i,k) = (pirhow*nrout(i,k) / qrout(i,k))**(-oneb3)
             freqr(i,k)  = one
           endif
           if (qsout(i,k) > 1.e-7_r8 .and. nsout(i,k) > zero) then
             qsout2(i,k) = qsout(i,k)
             nsout2(i,k) = nsout(i,k)
             dsout2(i,k) = (pirhosn*nsout(i,k) / qsout(i,k))**(-oneb3)
             freqs(i,k)  = one
           endif

! output activated liquid and ice (convert from #/kg -> #/m3)
           ncai(i,k) = dum2i(i,k)*rho(i,k)
           ncal(i,k) = dum2l(i,k)*rho(i,k)
         end do
       end do


#ifdef CAM
       call outfld('NCAL',ncal, pcols,lchnk)
       call outfld('NCAI',ncai, pcols,lchnk)

!add averaged output fields.
       call outfld('AQRAIN',qrout2, pcols,lchnk)
       call outfld('AQSNOW',qsout2, pcols,lchnk)
       call outfld('ANRAIN',nrout2, pcols,lchnk)
       call outfld('ANSNOW',nsout2, pcols,lchnk)
       call outfld('ADRAIN',drout2, pcols,lchnk)
       call outfld('ADSNOW',dsout2, pcols,lchnk)
       call outfld('FREQR',freqr, pcols,lchnk)
       call outfld('FREQS',freqs, pcols,lchnk)
#endif

!redefine fice here....
       do k=1,pver
         do i=1,ncol
           nfice(i,k) = zero
           tx1        = qc(i,k) + qctend(i,k)*deltat    ! water
           tx2        = qi(i,k) + qitend(i,k)*deltat    ! ice
           dumfice    = qsout(i,k) + qrout(i,k) + tx1 + tx2
           if (dumfice > zero) then
             nfice(i,k) = (qsout(i,k) + tx2) / dumfice
           endif
         enddo
       enddo

#ifdef CAM
       call outfld('FICE',nfice, pcols, lchnk)
#endif

      return
      end subroutine mmicro_pcond

#ifdef CAM

!##############################################################################

      subroutine findsp1 (lchnk, ncol, q, t, p, tsp, qsp)
!-----------------------------------------------------------------------
!
! Purpose:
!     find the wet bulb temperature for a given t and q
!     in a longitude height section
!     wet bulb temp is the temperature and spec humidity that is
!     just saturated and has the same enthalpy
!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method:
! a Newton method is used
! first guess uses an algorithm provided by John Petch from the UKMO
! we exclude points where the physical situation is unrealistic
! e.g. where the temperature is outside the range of validity for the
!      saturation vapor pressure, or where the water vapor pressure
!      exceeds the ambient pressure, or the saturation specific humidity is
!      unrealistic
!
! Author: P. Rasch
!
!-----------------------------------------------------------------------
!
!     input arguments
!
       integer, intent(in) :: lchnk
       integer, intent(in) :: ncol

       real(r8), intent(in) :: q(pcols,pver)
       real(r8), intent(in) :: t(pcols,pver)
       real(r8), intent(in) :: p(pcols,pver)
!
! output arguments
!
       real(r8), intent(out) :: tsp(pcols,pver)
       real(r8), intent(out) :: qsp(pcols,pver)
!
! local variables
!
       integer i
       integer k
       logical lflg
       integer iter
       integer l
       logical :: error_found

       real(r8) omeps
       real(r8) trinv
       real(r8) es
       real(r8) desdt
!     real(r8) desdp                ! change in sat vap pressure wrt pressure
       real(r8) dqsdt
       real(r8) dgdt
       real(r8) g
       real(r8) weight(pcols)
       real(r8) hlatsb
       real(r8) hlatvp
       real(r8) hltalt(pcols,pver)
       real(r8) tterm
       real(r8) qs
       real(r8) tc

! work variables
       real(r8) t1, q1, dt, dq
       real(r8) dtm, dqm
       real(r8) qvd, a1, tmp
       real(r8) rair
       real(r8) r1b, c1, c2, c3
       real(r8) denom
       real(r8) dttol
       real(r8) dqtol
       integer doit(pcols)
       real(r8) enin(pcols), enout(pcols)
       real(r8) tlim(pcols)

       omeps = one - epsqs
       trinv = one / ttrice
       a1 = 7.5_r8*log(10._r8)
       rair = 287.04_r8
       c3 = rair*a1/cp
       dtm = zero
       dqm = zero
       dttol = 1.e-4_r8
       dqtol = 1.e-4_r8
!  tmin = 173.16 ! the coldest temperature we can deal with
!
! max number of times to iterate the calculation
       iter = 10
!
       do k = 1,pver

!
! first guess on the wet bulb temperature
!
       do i = 1,ncol

#ifdef DEBUG
       if ( (lchnk == lchnklook(nlook) ) .and. (i ==
     & icollook(nlook) ) ) then
       write (iulog,*) ' '
       write (iulog,*) ' level, t, q, p', k, t(i,k), q(i,k), p(i,k)
       endif
#endif

! limit the temperature range to that relevant to the sat vap pres tables
#if ( defined WACCM_PHYS )
       tlim(i) = min(max(t(i,k),128._r8),373._r8)
#else
       tlim(i) = min(max(t(i,k),173._r8),373._r8)
#endif
       es = estblf(tlim(i))
       denom = p(i,k) - omeps*es
       qs = epsqs*es/denom
       doit(i) = 0
       enout(i) = one
! make sure a meaningful calculation is possible
       if (p(i,k) > 5._r8*es .and. qs > zero .and. qs < 0.5_r8) then
!
! Saturation specific humidity
!
       qs = min(epsqs*es/denom,one)
!
! "generalized" analytic expression for t derivative of es
!  accurate to within 1 percent for 173.16 < t < 373.16
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
       tc = tlim(i) - tt0
       lflg = (tc >= -ttrice .and. tc < zero)
       weight(i) = min(-tc*trinv,one)
       hlatsb = hlatv + weight(i)*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       if (tlim(i) < tt0) then
       hltalt(i,k) = hlatsb
       else
       hltalt(i,k) = hlatvp
       end if
       enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
       tmp = q(i,k) - qs
       c1 = hltalt(i,k)*c3
       c2 = (tlim(i) + 36._r8)**2
       r1b = c2/(c2 + c1*qs)
       qvd = r1b*tmp
       tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
#ifdef DEBUG
       if ( (lchnk == lchnklook(nlook) ) .and. (i ==
     & icollook(nlook) ) ) then
       write (iulog,*) ' relative humidity ', q(i,k)/qs
       write (iulog,*) ' first guess ', tsp(i,k)
       endif
#endif
       es = estblf(tsp(i,k))
       qsp(i,k) = min(epsqs*es/(p(i,k) - omeps*es),one)
       else
       doit(i) = 1
       tsp(i,k) = tlim(i)
       qsp(i,k) = q(i,k)
       enin(i) = one
       endif
       end do
!
! now iterate on first guess
!
       do l = 1, iter
       dtm = 0
       dqm = 0
       do i = 1,ncol
       if (doit(i) == 0) then
       es = estblf(tsp(i,k))
!
! Saturation specific humidity
!
       qs = min(epsqs*es/(p(i,k) - omeps*es),one)
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
       tc = tsp(i,k) - tt0
       lflg = (tc >= -ttrice .and. tc < zero)
       weight(i) = min(-tc*trinv,one)
       hlatsb = hlatv + weight(i)*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       if (tsp(i,k) < tt0) then
       hltalt(i,k) = hlatsb
       else
       hltalt(i,k) = hlatvp
       end if
       if (lflg) then
       tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*
     &pcf(5))))
       else
       tterm = zero
       end if
       desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
       dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
!              g = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
       g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
       dgdt = -(cp + hltalt(i,k)*dqsdt)
       t1 = tsp(i,k) - g/dgdt
       dt = abs(t1 - tsp(i,k))/t1
       tsp(i,k) = max(t1,tmin)
       es = estblf(tsp(i,k))
       q1 = min(epsqs*es/(p(i,k) - omeps*es),one)
       dq = abs(q1 - qsp(i,k))/max(q1,1.e-12_r8)
       qsp(i,k) = q1
#ifdef DEBUG
       if ( (lchnk == lchnklook(nlook) ) .and. (i ==
     & icollook(nlook) ) ) then
       write (iulog,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
       endif
#endif
       dtm = max(dtm,dt)
       dqm = max(dqm,dq)
! if converged at this point, exclude it from more iterations
       if (dt < dttol .and. dq < dqtol) then
       doit(i) = 2
       endif
       enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
! bail out if we are too near the end of temp range
#if ( defined WACCM_PHYS )
       if (tsp(i,k) < 130.16_r8) then
#else
       if (tsp(i,k) < 174.16_r8) then
#endif
       doit(i) = 4
       endif
       else
       endif
       end do

       if (dtm < dttol .and. dqm < dqtol) then
       go to 10
       endif

       end do
      10 continue

       error_found = .false.
       if (dtm > dttol .or. dqm > dqtol) then
       do i = 1,ncol
       if (doit(i) == 0) error_found = .true.
       end do
       if (error_found) then
       do i = 1,ncol
       if (doit(i) == 0) then
       write (iulog,*) ' findsp not converging at point i, k ', i, k
       write (iulog,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k),
     & enin(i)
       write (iulog,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k),
     & enout(i)
       call endrun ('FINDSP')
       endif
       end do
       endif
       endif
       do i = 1,ncol
       if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+
     &enout(i))) > 1.e-4_r8) then
       error_found = .true.
       endif
       end do
       if (error_found) then
       do i = 1,ncol
       if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+
     &enout(i))) > 1.e-4_r8) then
       write (iulog,*) ' the enthalpy is not conserved for point ', i,
     & k, enin(i), enout(i)
       write (iulog,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k),
     & enin(i)
       write (iulog,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k),
     & enout(i)
       call endrun ('FINDSP')
       endif
       end do
       endif

       end do

       return
      end subroutine findsp1


      subroutine findsp1_water (lchnk, ncol, q, t, p, tsp, qsp)
!-----------------------------------------------------------------------
!
! Purpose:
!     find the wet bulb temperature for a given t and q
!     in a longitude height section
!     wet bulb temp is the temperature and spec humidity that is
!     just saturated and has the same enthalpy
!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method:
! a Newton method is used
! first guess uses an algorithm provided by John Petch from the UKMO
! we exclude points where the physical situation is unrealistic
! e.g. where the temperature is outside the range of validity for the
!      saturation vapor pressure, or where the water vapor pressure
!      exceeds the ambient pressure, or the saturation specific humidity is
!      unrealistic
!
! Author: P. Rasch
!
!-----------------------------------------------------------------------
!
!     input arguments
!
       integer, intent(in) :: lchnk
       integer, intent(in) :: ncol

       real(r8), intent(in) :: q(pcols,pver)
       real(r8), intent(in) :: t(pcols,pver)
       real(r8), intent(in) :: p(pcols,pver)
!
! output arguments
!
       real(r8), intent(out) :: tsp(pcols,pver)
       real(r8), intent(out) :: qsp(pcols,pver)
!
! local variables
!
       integer i
       integer k
       logical lflg
       integer iter
       integer l
       logical :: error_found

       real(r8) omeps
       real(r8) trinv
       real(r8) es
       real(r8) desdt
!     real(r8) desdp                ! change in sat vap pressure wrt pressure
       real(r8) dqsdt
       real(r8) dgdt
       real(r8) g
       real(r8) weight(pcols)
       real(r8) hlatsb
       real(r8) hlatvp
       real(r8) hltalt(pcols,pver)
       real(r8) tterm
       real(r8) qs
       real(r8) tc

! work variables
       real(r8) t1, q1, dt, dq
       real(r8) dtm, dqm
       real(r8) qvd, a1, tmp
       real(r8) rair
       real(r8) r1b, c1, c2, c3
       real(r8) denom
       real(r8) dttol
       real(r8) dqtol
       integer doit(pcols)
       real(r8) enin(pcols), enout(pcols)
       real(r8) tlim(pcols)

       omeps = one - epsqs
       a1 = 7.5_r8*log(10._r8)
       rair = 287.04_r8
       c3 = rair*a1/cp
       dtm = zero
       dqm = zero
       dttol = 1.e-4_r8
       dqtol = 1.e-4_r8
!
! max number of times to iterate the calculation
       iter = 8
!
       do k = 1,pver

!
! first guess on the wet bulb temperature
!
       do i = 1,ncol

#ifdef DEBUG
       if ( (lchnk == lchnklook(nlook) ) .and. (i ==
     & icollook(nlook) ) ) then
       write (iulog,*) ' '
       write (iulog,*) ' level, t, q, p', k, t(i,k), q(i,k), p(i,k)
       endif
#endif
! limit the temperature range to that relevant to the sat vap pres tables
#if ( defined WACCM_PHYS )
       tlim(i) = min(max(t(i,k),128._r8),373._r8)
#else
       tlim(i) = min(max(t(i,k),173._r8),373._r8)
#endif

       es = polysvp(tlim(i),0)
       denom = p(i,k) - omeps*es
       qs = epsqs*es/denom
       doit(i) = 0
       enout(i) = one
! make sure a meaningful calculation is possible
       if (p(i,k) > 5._r8*es .and. qs > zero .and. qs < 0.5_r8) then
!
! Saturation specific humidity
!
       qs = min(epsqs*es/denom,one)
!
! "generalized" analytic expression for t derivative of es
!  accurate to within 1 percent for 173.16 < t < 373.16

!
! No icephs or water to ice transition
!
       hlatvp = hlatv - 2369.0*(tlim(i)-tt0)
       hlatsb = hlatv
       if (tlim(i) < tt0) then
       hltalt(i,k) = hlatsb
       else
       hltalt(i,k) = hlatvp
       end if
!--xl
       enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
       tmp = q(i,k) - qs
       c1 = hltalt(i,k)*c3
       c2 = (tlim(i) + 36._r8)**2
       r1b = c2/(c2 + c1*qs)
       qvd = r1b*tmp
       tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
#ifdef DEBUG
       if ( (lchnk == lchnklook(nlook) ) .and. (i ==
     & icollook(nlook) ) ) then
       write (iulog,*) ' relative humidity ', q(i,k)/qs
       write (iulog,*) ' first guess ', tsp(i,k)
       endif
#endif
       es = polysvp(tsp(i,k),0)
       qsp(i,k) = min(epsqs*es/(p(i,k) - omeps*es),one)
       else
       doit(i) = 1
       tsp(i,k) = tlim(i)
       qsp(i,k) = q(i,k)
       enin(i) = one
       endif
       end do
!
! now iterate on first guess
!
       do l = 1, iter
       dtm = 0
       dqm = 0
       do i = 1,ncol
       if (doit(i) == 0) then
       es = polysvp(tsp(i,k),0)
!
! Saturation specific humidity
!
       qs = min(epsqs*es/(p(i,k) - omeps*es),one)
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
!
! No icephs or water to ice transition
!
       hlatvp = hlatv - 2369.0*(tsp(i,k)-tt0)
       hlatsb = hlatv
       if (tsp(i,k) < tt0) then
       hltalt(i,k) = hlatsb
       else
       hltalt(i,k) = hlatvp
       end if
       desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k))
       dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
       g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
       dgdt = -(cp + hltalt(i,k)*dqsdt)
       t1 = tsp(i,k) - g/dgdt
       dt = abs(t1 - tsp(i,k))/t1
       tsp(i,k) = max(t1,tmin)

       es = polysvp(tsp(i,k),0)
       q1 = min(epsqs*es/(p(i,k) - omeps*es),one)
       dq = abs(q1 - qsp(i,k))/max(q1,1.e-12_r8)
       qsp(i,k) = q1
#ifdef DEBUG
       if ( (lchnk == lchnklook(nlook) ) .and. (i ==
     & icollook(nlook) ) ) then
       write (iulog,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
       endif
#endif
       dtm = max(dtm,dt)
       dqm = max(dqm,dq)
! if converged at this point, exclude it from more iterations
       if (dt < dttol .and. dq < dqtol) then
       doit(i) = 2
       endif
       enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
! bail out if we are too near the end of temp range
#if ( defined WACCM_PHYS )
       if (tsp(i,k) < 130.16_r8) then
#else
       if (tsp(i,k) < 174.16_r8) then
#endif
       doit(i) = 4
       endif
       else
       endif
       end do

       if (dtm < dttol .and. dqm < dqtol) then
       go to 10
       endif

       end do
      10 continue

       error_found = .false.
       if (dtm > dttol .or. dqm > dqtol) then
       do i = 1,ncol
       if (doit(i) == 0) error_found = .true.
       end do
       if (error_found) then
       do i = 1,ncol
       if (doit(i) == 0) then
       write (iulog,*) ' findsp not converging at point i, k ', i, k
       write (iulog,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k),
     & enin(i)
       write (iulog,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k),
     & enout(i)
       call endrun ('FINDSP')
       endif
       end do
       endif
       endif
       do i = 1,ncol
       if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+
     &enout(i))) > 1.e-4_r8) then
       error_found = .true.
       endif
       end do
       if (error_found) then
       do i = 1,ncol
       if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+
     &enout(i))) > 1.e-4_r8) then
       write (iulog,*) ' the enthalpy is not conserved for point ', i,
     & k, enin(i), enout(i)
       write (iulog,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k),
     & enin(i)
       write (iulog,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k),
     & enout(i)
       call endrun ('FINDSP')
       endif
       end do
       endif

       end do

       return
      end subroutine findsp1_water
#endif
!--jtb : end of GEOS5 exclusion (begins at top of findsp)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! error function in single precision
!
!    Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
!    You may use, copy, modify this code for any purpose and
!    without fee. You may distribute this ORIGINAL package.

       function derf(x)
       implicit real (a - h, o - z)
       real(r8) a,b,x
       dimension a(0 : 64), b(0 : 64)
       integer i,k
       data (a(i), i = 0, 12) / 0.00000000005958930743d0, -
     &0.00000000113739022964d0, 0.00000001466005199839d0, -
     &0.00000016350354461960d0, 0.00000164610044809620d0, -
     &0.00001492559551950604d0, 0.00012055331122299265d0, -
     &0.00085483269811296660d0, 0.00522397762482322257d0, -
     &0.02686617064507733420d0, 0.11283791670954881569d0, -
     &0.37612638903183748117d0, 1.12837916709551257377d0 /
       data (a(i), i = 13, 25) / 0.00000000002372510631d0, -
     &0.00000000045493253732d0, 0.00000000590362766598d0, -
     &0.00000006642090827576d0, 0.00000067595634268133d0, -
     &0.00000621188515924000d0, 0.00005103883009709690d0, -
     &0.00037015410692956173d0, 0.00233307631218880978d0, -
     &0.01254988477182192210d0, 0.05657061146827041994d0, -
     &0.21379664776456006580d0, 0.84270079294971486929d0 /
       data (a(i), i = 26, 38) / 0.00000000000949905026d0, -
     &0.00000000018310229805d0, 0.00000000239463074000d0, -
     &0.00000002721444369609d0, 0.00000028045522331686d0, -
     &0.00000261830022482897d0, 0.00002195455056768781d0, -
     &0.00016358986921372656d0, 0.00107052153564110318d0, -
     &0.00608284718113590151d0, 0.02986978465246258244d0, -
     &0.13055593046562267625d0, 0.67493323603965504676d0 /
       data (a(i), i = 39, 51) / 0.00000000000382722073d0, -
     &0.00000000007421598602d0, 0.00000000097930574080d0, -
     &0.00000001126008898854d0, 0.00000011775134830784d0, -
     &0.00000111992758382650d0, 0.00000962023443095201d0, -
     &0.00007404402135070773d0, 0.00050689993654144881d0, -
     &0.00307553051439272889d0, 0.01668977892553165586d0, -
     &0.08548534594781312114d0, 0.56909076642393639985d0 /
       data (a(i), i = 52, 64) / 0.00000000000155296588d0, -
     &0.00000000003032205868d0, 0.00000000040424830707d0, -
     &0.00000000471135111493d0, 0.00000005011915876293d0, -
     &0.00000048722516178974d0, 0.00000430683284629395d0, -
     &0.00003445026145385764d0, 0.00024879276133931664d0, -
     &0.00162940941748079288d0, 0.00988786373932350462d0, -
     &0.05962426839442303805d0, 0.49766113250947636708d0 /
       data (b(i), i = 0, 12) / -0.00000000029734388465d0,
     & 0.00000000269776334046d0, -0.00000000640788827665d0, -
     &0.00000001667820132100d0, -0.00000021854388148686d0,
     & 0.00000266246030457984d0, 0.00001612722157047886d0, -
     &0.00025616361025506629d0, 0.00015380842432375365d0,
     & 0.00815533022524927908d0, -0.01402283663896319337d0, -
     &0.19746892495383021487d0, 0.71511720328842845913d0 /
       data (b(i), i = 13, 25) / -0.00000000001951073787d0, -
     &0.00000000032302692214d0, 0.00000000522461866919d0,
     & 0.00000000342940918551d0, -0.00000035772874310272d0,
     & 0.00000019999935792654d0, 0.00002687044575042908d0, -
     &0.00011843240273775776d0, -0.00080991728956032271d0,
     & 0.00661062970502241174d0, 0.00909530922354827295d0, -
     &0.20160072778491013140d0, 0.51169696718727644908d0 /
       data (b(i), i = 26, 38) / 0.00000000003147682272d0, -
     &0.00000000048465972408d0, 0.00000000063675740242d0,
     & 0.00000003377623323271d0, -0.00000015451139637086d0, -
     &0.00000203340624738438d0, 0.00001947204525295057d0,
     & 0.00002854147231653228d0, -0.00101565063152200272d0,
     & 0.00271187003520095655d0, 0.02328095035422810727d0, -
     &0.16725021123116877197d0, 0.32490054966649436974d0 /
       data (b(i), i = 39, 51) / 0.00000000002319363370d0, -
     &0.00000000006303206648d0, -0.00000000264888267434d0,
     & 0.00000002050708040581d0, 0.00000011371857327578d0, -
     &0.00000211211337219663d0, 0.00000368797328322935d0,
     & 0.00009823686253424796d0, -0.00065860243990455368d0, -
     &0.00075285814895230877d0, 0.02585434424202960464d0, -
     &0.11637092784486193258d0, 0.18267336775296612024d0 /
       data (b(i), i = 52, 64) / -0.00000000000367789363d0,
     & 0.00000000020876046746d0, -0.00000000193319027226d0, -
     &0.00000000435953392472d0, 0.00000018006992266137d0, -
     &0.00000078441223763969d0, -0.00000675407647949153d0,
     & 0.00008428418334440096d0, -0.00017604388937031815d0, -
     &0.00239729611435071610d0, 0.02064129023876022970d0, -
     &0.06905562880005864105d0, 0.09084526782065478489d0 /
       w = abs(x)
       if (w .lt. 2.2d0) then
       t = w * w
       k = int(t)
       t = t - k
       k = k * 13
       y = ((((((((((((a(k) * t + a(k + 1)) * t + a(k + 2)) * t + a(k +
     & 3)) * t + a(k + 4)) * t + a(k + 5)) * t + a(k + 6)) * t + a(k +
     & 7)) * t + a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + a(k +
     & 11)) * t + a(k + 12)) * w
       else if (w .lt. 6.9d0) then
       k = int(w)
       t = w - k
       k = 13 * (k - 2)
       y = (((((((((((b(k) * t + b(k + 1)) * t + b(k + 2)) * t + b(k +
     & 3)) * t + b(k + 4)) * t + b(k + 5)) * t + b(k + 6)) * t + b(k +
     & 7)) * t + b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + b(k +
     & 11)) * t + b(k + 12)
       y = y * y
       y = y * y
       y = y * y
       y = 1 - y * y
       else
       y = 1
       end if
       if (x .lt. 0) y = -y
       derf = y
       end function derf
!


!cccccccccccccccccccccDONIFccccccccccccccccccccccccccccccccccccccccccccccccc



!**********************************
       FUNCTION MUI_HEMP(T)


       real(r8) :: MUI_HEMP
       REAL(r8), intent(in) :: T
       REAL(r8) :: TC, mui, lambdai
       TC=T-273.15_r8

       TC=MIN(MAX(TC, -70.0), -15.0)

       if (TC  >  -27.0) then
       lambdai = 6.8_r8*exp(-0.096_r8*TC)
       else
       lambdai = 24.8_r8*exp(-0.049_r8*TC)
       end if

       mui=(0.13_r8*(lambdai**0.64_r8))-two
       mui=max(mui, 1.5_r8)
       MUI_HEMP=mui


       END FUNCTION MUI_HEMP


!cccccccccccccccccccccDONIFccccccccccccccccccccccccccccccccccccccccccccccccc



!**********************************
       FUNCTION MUI_HEMP_L(lambda)


       real(r8) :: MUI_HEMP_L
       REAL(r8), intent(in) :: lambda
       REAL(r8) :: TC, mui, lx
       lx = lambda*0.01

       mui=(0.008_r8*(lx**0.87_r8))
       mui=max(min(mui, 10.0_r8), 0.1_r8)
       MUI_HEMP_L=mui  !Anning for multithread to work


       END FUNCTION MUI_HEMP_L




       FUNCTION gamma_incomp(muice, x)



       real(r8)             :: gamma_incomp
       REAL(r8), intent(in) :: muice, x
       REAL(r8)             :: xog, kg, alfa, auxx
       alfa = min(max(muice+1._r8, 1._r8), 20._r8)

       xog  = log(alfa -0.3068_r8)
       kg   = 1.44818*(alfa**0.5357_r8)
       auxx = max(min(kg*(log(x)-xog), 30._r8), -30._r8)
       gamma_incomp= one/(one +exp(-auxx))
       gamma_incomp = max(gamma_incomp, 1.0e-20)

       END FUNCTION gamma_incomp


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       FUNCTION GAMMA(X)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!D    DOUBLE PRECISION FUNCTION DGAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
       INTEGER I,N
       LOGICAL PARITY

       REAL(r8)  :: gamma,CONV,FACT,RES,SUM,X,XDEN,XNUM,Y,Y1,YSQ,Z
       real(r8) ::  C(7),P(8),Q(8)
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
       real(r8), parameter :: ONE=1.0E0_r8,     HALF=0.5E0_r8,          &
     &                        TWELVE=12.0E0_r8, TWO=2.0E0_r8,           &
     &                        ZERO=0.0E0_r8,                            &
     &                        PI=3.1415926535897932384626434E0_r8,      &
     &                        SQRTPI=0.9189385332046727417803297E0_r8

!D    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
!D   1     SQRTPI/0.9189385332046727417803297D0/,
!D   2     PI/3.1415926535897932384626434D0/
!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
       real(r8), parameter :: XBIG=35.040E0_r8, XMININ=1.18E-38_r8,     &
     &                        EPS=1.19E-7_r8,   XINF=3.4E38_r8

!D    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
!D   1     XINF/1.79D308/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
       DATA P/-1.71618513886549492533811E+0_r8,
     &2.47656508055759199108314E+1_r8, -3.79804256470945635097577E+2_r8,
     &6.29331155312818442661052E+2_r8, 8.66966202790413211295064E+2_r8,-
     &3.14512729688483675254357E+4_r8, -3.61444134186911729807069E+4_r8,
     &6.64561438202405440627855E+4_r8/
       DATA Q/-3.08402300119738975254353E+1_r8,
     &3.15350626979604161529144E+2_r8, -1.01515636749021914166146E+3_r8,
     &-3.10777167157231109440444E+3_r8, 2.25381184209801510330112E+4_r8,
     &4.75584627752788110767815E+3_r8, -1.34659959864969306392456E+5_r8,
     &-1.15132259675553483497211E+5_r8/
!D    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
!D   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
!D   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
!D   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
!D    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
!D   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
!D   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
!D   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
       DATA C/-1.910444077728E-03_r8,8.4171387781295E-04_r8, -
     &5.952379913043012E-04_r8,7.93650793500350248E-04_r8, -
     &2.777777777777681622553E-03_r8,8.333333333333333331554247E-02_r8,
     & 5.7083835261E-03_r8/
!D    DATA C/-1.910444077728D-03,8.4171387781295D-04,
!D   1     -5.952379913043012D-04,7.93650793500350248D-04,
!D   2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
!D   3      5.7083835261D-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I,r8)
!D    CONV(I) = DBLE(I)
       PARITY = .FALSE.
       FACT   = ONE
       N = 0
       Y = X
       IF(Y <= ZERO) THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
         Y  = -X
         Y1 = AINT(Y)
         RES = Y - Y1
         IF(RES.NE.ZERO)THEN
           IF(Y1.NE.AINT(Y1*HALF)*TWO) PARITY = .TRUE.
             FACT = -PI/SIN(PI*RES)
             Y    = Y + ONE
         ELSE
           RES=XINF
           GOTO 900
         ENDIF
       ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
       IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
         IF(Y.GE.XMININ)THEN
           RES = ONE/Y
         ELSE
           RES = XINF
           GOTO 900
         ENDIF
       ELSEIF(Y.LT.TWELVE)THEN
         Y1 = Y
         IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
           Z = Y
           Y = Y + ONE
         ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
           N = INT(Y) - 1
           Y = Y - CONV(N)
           Z = Y - ONE
         ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
         XNUM = ZERO
         XDEN = ONE
         DO 260 I=1,8
           XNUM = (XNUM+P(I))*Z
           XDEN = XDEN*Z + Q(I)
260      CONTINUE
         RES = XNUM/XDEN + ONE
         IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
           RES = RES/Y1
         ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
           DO 290 I=1,N
             RES = RES*Y
             Y   = Y + ONE
290        CONTINUE
         ENDIF
       ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
         IF(Y.LE.XBIG)THEN
           YSQ = Y*Y
           SUM = C(7)
           DO 350 I=1,6
             SUM = SUM / YSQ + C(I)
350        CONTINUE
           SUM = SUM / Y - Y + SQRTPI
           SUM = SUM + (Y-HALF)*LOG(Y)
           RES = EXP(SUM)
         ELSE
           RES = XINF
           GOTO 900
         ENDIF
       ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
       IF(PARITY) RES = -RES
       IF(FACT.NE.ONE) RES = FACT/RES
900    GAMMA = RES
!D900 DGAMMA = RES
       RETURN
! ---------- LAST LINE OF GAMMA ----------
       END function gamma



      end module cldwat2m_micro
