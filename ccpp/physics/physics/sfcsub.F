      module sfccyc_module
      implicit none
      save
!
!  grib code for each parameter - used in subroutines sfccycle and setrmsk.
!
      integer kpdtsf,kpdwet,kpdsno,kpdzor,kpdais,kpdtg3,kpdplr,kpdgla,
     &        kpdmxi,kpdscv,kpdsmc,kpdoro,kpdmsk,kpdstc,kpdacn,kpdveg,
     &        kpdvet,kpdsot
     &,       kpdvmn,kpdvmx,kpdslp,kpdabs
     &,       kpdsnd, kpdabs_0, kpdabs_1, kpdalb(4)
      parameter(kpdtsf=11,  kpdwet=86, kpdsno=65,  kpdzor=83,
!    1          kpdalb=84,  kpdais=91, kpdtg3=11,  kpdplr=224,
     1          kpdais=91,  kpdtg3=11, kpdplr=224,
     2          kpdgla=238, kpdmxi=91, kpdscv=238, kpdsmc=144,
     3          kpdoro=8,   kpdmsk=81, kpdstc=11,  kpdacn=91, kpdveg=87,
!cbosu  max snow albedo uses a grib id number of 159, not 255.
     &          kpdvmn=255, kpdvmx=255,kpdslp=236, kpdabs_0=255,    
     &          kpdvet=225, kpdsot=224,kpdabs_1=159,
     &          kpdsnd=66 )
!
      integer, parameter :: kpdalb_0(4)=(/212,215,213,216/)
      integer, parameter :: kpdalb_1(4)=(/189,190,191,192/)
      integer, parameter :: kpdalf(2)=(/214,217/)
!
      integer, parameter :: xdata=5000, ydata=2500, mdata=xdata*ydata
      integer            :: veg_type_landice
      integer            :: soil_type_landice
!
      end module sfccyc_module
      subroutine sfccycle_sub(lugb,len,lsoil,sig1t,deltsfc
     &,                   iy,im,id,ih,fh
     &,                   rla, rlo, slmask,orog,orog_uf,use_ufo,nst_anl
     &,                   sihfcs,sicfcs,sitfcs                 
     &,                   swdfcs,slcfcs      
     &,                   vmnfcs,vmxfcs,slpfcs,absfcs
     &,                   tsffcs,snofcs,zorfcs,albfcs,tg3fcs
     &,                   cnpfcs,smcfcs,stcfcs,slifcs,aisfcs,f10m
     &,                   vegfcs,vetfcs,sotfcs,alffcs
     &,                   cvfcs,cvbfcs,cvtfcs,me,nlunit,ialb
     &,                   isot,ivegsrc)
!
      use machine , only : kind_io8,kind_io4
      use sfccyc_module
      implicit none
      logical use_ufo, nst_anl
      real (kind=kind_io8) sllnd,slsea,aicice,aicsea,tgice,rlapse,
     &                     orolmx,orolmn,oroomx,oroomn,orosmx,
     &                     orosmn,oroimx,oroimn,orojmx,orojmn,
     &                     alblmx,alblmn,albomx,albomn,albsmx,
     &                     albsmn,albimx,albimn,albjmx,albjmn,
     &                     wetlmx,wetlmn,wetomx,wetomn,wetsmx,
     &                     wetsmn,wetimx,wetimn,wetjmx,wetjmn,
     &                     snolmx,snolmn,snoomx,snoomn,snosmx,
     &                     snosmn,snoimx,snoimn,snojmx,snojmn,
     &                     zorlmx,zorlmn,zoromx,zoromn,zorsmx,
     &                     zorsmn,zorimx,zorimn,zorjmx, zorjmn,
     &                     plrlmx,plrlmn,plromx,plromn,plrsmx,
     &                     plrsmn,plrimx,plrimn,plrjmx,plrjmn,
     &                     tsflmx,tsflmn,tsfomx,tsfomn,tsfsmx,
     &                     tsfsmn,tsfimx,tsfimn,tsfjmx,tsfjmn,
     &                     tg3lmx,tg3lmn,tg3omx,tg3omn,tg3smx,
     &                     tg3smn,tg3imx,tg3imn,tg3jmx,tg3jmn,
     &                     stclmx,stclmn,stcomx,stcomn,stcsmx,
     &                     stcsmn,stcimx,stcimn,stcjmx,stcjmn,
     &                     smclmx,smclmn,smcomx,smcomn,smcsmx,
     &                     smcsmn,smcimx,smcimn,smcjmx,smcjmn,
     &                     scvlmx,scvlmn,scvomx,scvomn,scvsmx,
     &                     scvsmn,scvimx,scvimn,scvjmx,scvjmn,
     &                     veglmx,veglmn,vegomx,vegomn,vegsmx,
     &                     vegsmn,vegimx,vegimn,vegjmx,vegjmn,
     &                     vetlmx,vetlmn,vetomx,vetomn,vetsmx,
     &                     vetsmn,vetimx,vetimn,vetjmx,vetjmn,
     &                     sotlmx,sotlmn,sotomx,sotomn,sotsmx,
     &                     sotsmn,sotimx,sotimn,sotjmx,sotjmn,
     &                     alslmx,alslmn,alsomx,alsomn,alssmx,
     &                     alssmn,alsimx,alsimn,alsjmx,alsjmn,
     &                     epstsf,epsalb,epssno,epswet,epszor,
     &                     epsplr,epsoro,epssmc,epsscv,eptsfc,
     &                     epstg3,epsais,epsacn,epsveg,epsvet,
     &                     epssot,epsalf,qctsfs,qcsnos,qctsfi,
     &                     aislim,snwmin,snwmax,cplrl,cplrs,
     &                     cvegl,czors,csnol,csnos,czorl,csots,
     &                     csotl,cvwgs,cvetl,cvets,calfs,
     &                     fcalfl,fcalfs,ccvt,ccnp,ccv,ccvb,
     &                     calbl,calfl,calbs,ctsfs,grboro,
     &                     grbmsk,ctsfl,deltf,caisl,caiss,
     &                     fsalfl,fsalfs,flalfs,falbl,ftsfl,
     &                     ftsfs,fzorl,fzors,fplrl,fsnos,faisl,
     &                     faiss,fsnol,bltmsk,falbs,cvegs,percrit,
     &                     deltsfc,critp2,critp3,blnmsk,critp1,
     &                     fcplrl,fcplrs,fczors,fvets,fsotl,fsots,
     &                     fvetl,fplrs,fvegl,fvegs,fcsnol,fcsnos,
     &                     fczorl,fcalbs,fctsfl,fctsfs,fcalbl,
     &                     falfs,falfl,fh,crit,zsca,ztsfc,tem1,tem2
     &,                    fsihl,fsihs,fsicl,fsics,
     &                     csihl,csihs,csicl,csics,epssih,epssic
     &,                    fvmnl,fvmns,fvmxl,fvmxs,fslpl,fslps,
     &                     fabsl,fabss,cvmnl,cvmns,cvmxl,cvmxs,
     &                     cslpl,cslps,cabsl,cabss,epsvmn,epsvmx,
     &                     epsslp,epsabs
     &,                    sihlmx,sihlmn,sihomx,sihomn,sihsmx,
     &                     sihsmn,sihimx,sihimn,sihjmx,sihjmn,
     &                     siclmx,siclmn,sicomx,sicomn,sicsmx,
     &                     sicsmn,sicimx,sicimn,sicjmx,sicjmn
     &,                    glacir_hice
     &,                    vmnlmx,vmnlmn,vmnomx,vmnomn,vmnsmx,
     &                     vmnsmn,vmnimx,vmnimn,vmnjmx,vmnjmn,
     &                     vmxlmx,vmxlmn,vmxomx,vmxomn,vmxsmx,
     &                     vmxsmn,vmximx,vmximn,vmxjmx,vmxjmn,
     &                     slplmx,slplmn,slpomx,slpomn,slpsmx,
     &                     slpsmn,slpimx,slpimn,slpjmx,slpjmn,
     &                     abslmx,abslmn,absomx,absomn,abssmx,
     &                     abssmn,absimx,absimn,absjmx,absjmn
     &,                    sihnew

      integer imsk,jmsk,ifp,irtscv,irtacn,irtais,irtsno,irtzor,
     &        irtalb,irtsot,irtalf,j,irtvet,irtsmc,irtstc,irtveg,
     &        irtwet,k,iprnt,kk,irttsf,iret,i,igrdbg,iy,im,id,
     &        icalbl,icalbs,icalfl,ictsfs,lugb,len,lsoil,ih,
     &        ictsfl,iczors,icplrl,icplrs,iczorl,icalfs,icsnol,
     &        icsnos,irttg3,me,kqcm, nlunit,ialb
     &,       irtvmn, irtvmx, irtslp, irtabs, isot, ivegsrc
      logical gausm, deads, qcmsk, znlst, monclm, monanl,
     &        monfcs, monmer, mondif, landice

      integer num_parthds
!
!  this is a limited point version of surface program.
!
!  this program runs in two different modes:
!
!  1.  analysis mode (fh=0.)
!
!      this program merges climatology, analysis and forecast guess to create
!      new surface fields.  if analysis file is given, the program
!      uses it if date of the analysis matches with iy,im,id,ih (see note
!      below).
!
!  2.  forecast mode (fh.gt.0.)
!
!      this program interpolates climatology to the date corresponding to the
!      forecast hour.  if surface analysis file is given, for the corresponding
!      dates, the program will use it.
!
!   note:
!
!      if the date of the analysis does not match given iy,im,id,ih, (and fh),
!      the program searches an old analysis by going back 6 hours, then 12 hours,
!      then one day upto nrepmx days (parameter statement in the subrotine fixrd.
!      now defined as 8).  this allows the user to provide non-daily analysis to
!      be used.  if matching field is not found, the forecast guess will be used.
!
!      use of a combined earlier surface analyses and current analysis is
!      not allowed (as was done in the old version for snow analysis in which
!      old snow analysis is used in combination with initial guess), except
!      for sea surface temperature.  for sst anolmaly interpolation, you need to
!      set lanom=.true. and must provide sst analysis at initial time.
!
!      if you want to do complex merging of past and present surface field analysis,
!      you need to create a separate file that contains daily surface field.
!
!      for a dead start, do not supply fnbgsi or set fnbgsi='        '
!
!  lugb           is the unit number used in this subprogram
!  len ...        number of points on which sfccyc operates
!  lsoil .. 	  number of soil layers (2 as of april, 1994)
!  iy,im,id,ih .. year, month, day, and hour of initial state.
!  fh ..          forecast hour
!  rla, rlo --    latitude and longitudes of the len points
!  sig1t .. sigma level 1 temperature for dead start.  should be on gaussian
!           grid.  if not dead start, no need for dimension but set to zero
!           as in the example below.
!
!  variable naming conventions:
!
!     oro .. orography
!     alb .. albedo
!     wet .. soil wetness as defined for bucket model
!     sno .. snow depth
!     zor .. surface roughness length
!     vet .. vegetation type
!     plr .. plant evaporation resistance
!     tsf .. surface skin temperature.  sea surface temp. over ocean.
!     tg3 .. deep soil temperature (at 500cm)
!     stc .. soil temperature (lsoil layrs)
!     smc .. soil moisture (lsoil layrs)
!     scv .. snow cover (not snow depth)
!     ais .. sea ice mask (0 or 1)
!     acn .. sea ice concentration (fraction)
!     gla .. glacier (permanent snow) mask (0 or 1)
!     mxi .. maximum sea ice extent (0 or 1)
!     msk .. land ocean mask (0=ocean 1=land)
!     cnp .. canopy water content
!     cv  .. convective cloud cover
!     cvb .. convective cloud base
!     cvt .. convective cloud top
!     sli .. land/sea/sea-ice mask. (1/0/2 respectively)
!     veg .. vegetation cover
!     sot .. soil type
!cwu [+2l] add sih & sic
!     sih .. sea ice thickness
!     sic .. sea ice concentration
!clu [+6l] add swd,slc,vmn,vmx,slp,abs
!     swd .. actual snow depth
!     slc .. liquid soil moisture (lsoil layers)
!     vmn .. vegetation cover minimum
!     vmx .. vegetation cover maximum
!     slp .. slope type
!     abs .. maximum snow albedo

!
!  definition of land/sea mask. sllnd for land and slsea for sea.
!  definition of sea/ice mask. aicice for ice, aicsea for sea.
!  tgice=max ice temperature
!  rlapse=lapse rate for sst correction due to surface angulation
!
      parameter(sllnd =1.0,slsea =0.0)
      parameter(aicice=1.0,aicsea=0.0)
      parameter(tgice=271.2)
      parameter(rlapse=0.65e-2)
!
!  max/min of fields for check and replace.
!
!     ???lmx .. max over bare land
!     ???lmn .. min over bare land
!     ???omx .. max over open ocean
!     ???omn .. min over open ocean
!     ???smx .. max over snow surface (land and sea-ice)
!     ???smn .. min over snow surface (land and sea-ice)
!     ???imx .. max over bare sea ice
!     ???imn .. min over bare sea ice
!     ???jmx .. max over snow covered sea ice
!     ???jmn .. min over snow covered sea ice
!
      parameter(orolmx=8000.,orolmn=-1000.,oroomx=3000.,oroomn=-1000.,
     &          orosmx=8000.,orosmn=-1000.,oroimx=3000.,oroimn=-1000.,
     &          orojmx=3000.,orojmn=-1000.)
!     parameter(alblmx=0.80,alblmn=0.06,albomx=0.06,albomn=0.06,
!    &          albsmx=0.80,albsmn=0.06,albimx=0.80,albimn=0.80,
!    &          albjmx=0.80,albjmn=0.80)
!cwu [-3l/+9l] change min/max for alb; add min/max for sih & sic
!     parameter(alblmx=0.80,alblmn=0.01,albomx=0.01,albomn=0.01,
!    &          albsmx=0.80,albsmn=0.01,albimx=0.01,albimn=0.01,
!    &          albjmx=0.01,albjmn=0.01)
!  note: the range values for bare land and snow covered land
!        (alblmx, alblmn, albsmx, albsmn) are set below
!        based on whether the old or new radiation is selected
      parameter(albomx=0.06,albomn=0.06,
     &          albimx=0.80,albimn=0.06,
     &          albjmx=0.80,albjmn=0.06)
      parameter(sihlmx=0.0,sihlmn=0.0,sihomx=5.0,sihomn=0.0,
     &          sihsmx=5.0,sihsmn=0.0,sihimx=5.0,sihimn=0.10,
     &          sihjmx=5.0,sihjmn=0.10,glacir_hice=3.0)
!cwu change sicimn & sicjmn Jan 2015
!     parameter(siclmx=0.0,siclmn=0.0,sicomx=1.0,sicomn=0.0,
!    &          sicsmx=1.0,sicsmn=0.0,sicimx=1.0,sicimn=0.50,
!    &          sicjmx=1.0,sicjmn=0.50)
!
!     parameter(sihlmx=0.0,sihlmn=0.0,sihomx=8.0,sihomn=0.0,
!    &          sihsmx=8.0,sihsmn=0.0,sihimx=8.0,sihimn=0.10,
!    &          sihjmx=8.0,sihjmn=0.10,glacir_hice=3.0)
      parameter(siclmx=0.0,siclmn=0.0,sicomx=1.0,sicomn=0.0,
     &          sicsmx=1.0,sicsmn=0.0,sicimx=1.0,sicimn=0.15,
     &          sicjmx=1.0,sicjmn=0.15)

      parameter(wetlmx=0.15,wetlmn=0.00,wetomx=0.15,wetomn=0.15,
     &          wetsmx=0.15,wetsmn=0.15,wetimx=0.15,wetimn=0.15,
     &          wetjmx=0.15,wetjmn=0.15)
      parameter(snolmx=0.0,snolmn=0.0,snoomx=0.0,snoomn=0.0,
     &          snosmx=55000.,snosmn=0.001,snoimx=0.,snoimn=0.0,
     &          snojmx=10000.,snojmn=0.01)
      parameter(zorlmx=300.,zorlmn=1.0,zoromx=1.0,zoromn=1.e-05,
     &          zorsmx=300.,zorsmn=1.0,zorimx=1.0,zorimn=1.0,
     &          zorjmx=1.0,zorjmn=1.0)
      parameter(plrlmx=1000.,plrlmn=0.0,plromx=1000.0,plromn=0.0,
     &          plrsmx=1000.,plrsmn=0.0,plrimx=1000.,plrimn=0.0,
     &          plrjmx=1000.,plrjmn=0.0)
!clu [-1l/+1l] relax tsfsmx (for noah lsm)
      parameter(tsflmx=353.,tsflmn=173.0,tsfomx=313.0,tsfomn=271.2,
     &          tsfsmx=305.0,tsfsmn=173.0,tsfimx=271.2,tsfimn=173.0,
     &          tsfjmx=273.16,tsfjmn=173.0)
!     parameter(tsflmx=353.,tsflmn=173.0,tsfomx=313.0,tsfomn=271.21,
!*   &          tsfsmx=273.16,tsfsmn=173.0,tsfimx=271.21,tsfimn=173.0,
!    &          tsfsmx=305.0,tsfsmn=173.0,tsfimx=271.21,tsfimn=173.0,
      parameter(tg3lmx=310.,tg3lmn=200.0,tg3omx=310.0,tg3omn=200.0,
     &          tg3smx=310.,tg3smn=200.0,tg3imx=310.0,tg3imn=200.0,
     &          tg3jmx=310.,tg3jmn=200.0)
      parameter(stclmx=353.,stclmn=173.0,stcomx=313.0,stcomn=200.0,
     &          stcsmx=310.,stcsmn=200.0,stcimx=310.0,stcimn=200.0,
     &          stcjmx=310.,stcjmn=200.0)
!landice mods   force a flag value of soil moisture of 1.0
!               at non-land points
      parameter(smclmx=0.55,smclmn=0.0,smcomx=1.0,smcomn=1.0,
     &          smcsmx=0.55,smcsmn=0.0,smcimx=1.0,smcimn=1.0,
     &          smcjmx=1.0,smcjmn=1.0)
      parameter(scvlmx=0.0,scvlmn=0.0,scvomx=0.0,scvomn=0.0,
     &          scvsmx=1.0,scvsmn=1.0,scvimx=0.0,scvimn=0.0,
     &          scvjmx=1.0,scvjmn=1.0)
      parameter(veglmx=1.0,veglmn=0.0,vegomx=0.0,vegomn=0.0,
     &          vegsmx=1.0,vegsmn=0.0,vegimx=0.0,vegimn=0.0,
     &          vegjmx=0.0,vegjmn=0.0)
      parameter(vmnlmx=1.0,vmnlmn=0.0,vmnomx=0.0,vmnomn=0.0,
     &          vmnsmx=1.0,vmnsmn=0.0,vmnimx=0.0,vmnimn=0.0,
     &          vmnjmx=0.0,vmnjmn=0.0)   
      parameter(vmxlmx=1.0,vmxlmn=0.0,vmxomx=0.0,vmxomn=0.0,
     &          vmxsmx=1.0,vmxsmn=0.0,vmximx=0.0,vmximn=0.0,
     &          vmxjmx=0.0,vmxjmn=0.0)  
      parameter(slplmx=9.0,slplmn=1.0,slpomx=0.0,slpomn=0.0,
     &          slpsmx=9.0,slpsmn=1.0,slpimx=0.,slpimn=0.,
     &          slpjmx=0.,slpjmn=0.) 
!  note: the range values for bare land and snow covered land
!        (alblmx, alblmn, albsmx, albsmn) are set below
!        based on whether the old or new radiation is selected
      parameter(absomx=0.0,absomn=0.0,
     &          absimx=0.0,absimn=0.0,
     &          absjmx=0.0,absjmn=0.0)    
!  vegetation type
      parameter(vetlmx=20.,vetlmn=1.0,vetomx=0.0,vetomn=0.0,
     &          vetsmx=20.,vetsmn=1.0,vetimx=0.,vetimn=0.,
     &          vetjmx=0.,vetjmn=0.)
!  soil type
      parameter(sotlmx=16.,sotlmn=1.0,sotomx=0.0,sotomn=0.0,
     &          sotsmx=16.,sotsmn=1.0,sotimx=0.,sotimn=0.,
     &          sotjmx=0.,sotjmn=0.)
!  fraction of vegetation for strongly and weakly zeneith angle dependent
!  albedo
      parameter(alslmx=1.0,alslmn=0.0,alsomx=0.0,alsomn=0.0,
     &          alssmx=1.0,alssmn=0.0,alsimx=0.0,alsimn=0.0,
     &          alsjmx=0.0,alsjmn=0.0)
!
!  criteria used for monitoring
!
      parameter(epstsf=0.01,epsalb=0.001,epssno=0.01,
     &          epswet=0.01,epszor=0.0000001,epsplr=1.,epsoro=0.,
     &          epssmc=0.0001,epsscv=0.,eptsfc=0.01,epstg3=0.01,
     &          epsais=0.,epsacn=0.01,epsveg=0.01,
     &          epssih=0.001,epssic=0.001,
     &          epsvmn=0.01,epsvmx=0.01,epsabs=0.001,epsslp=0.01,
     &          epsvet=.01,epssot=.01,epsalf=.001)
!
!  quality control of analysis snow and sea ice
!
!   qctsfs .. surface temperature above which no snow allowed
!   qcsnos .. snow depth above which snow must exist
!   qctsfi .. sst above which sea-ice is not allowed
!
!clu relax qctsfs (for noah lsm)
!*    parameter(qctsfs=283.16,qcsnos=100.,qctsfi=280.16)
!*    parameter(qctsfs=288.16,qcsnos=100.,qctsfi=280.16)
      parameter(qctsfs=293.16,qcsnos=100.,qctsfi=280.16)
!
!cwu [-2l]
!* ice concentration for ice limit (55 percent)
!
!*    parameter(aislim=0.55)
!
!  parameters to obtain snow depth from snow cover and temperature
!
!     parameter(snwmin=25.,snwmax=100.)
      parameter(snwmin=5.0,snwmax=100.)
      real (kind=kind_io8), parameter :: ten=10.0, one=1.0
!
!  coeeficients of blending forecast and interpolated clim
!  (or analyzed) fields over sea or land(l) (not for clouds)
!  1.0 = use of forecast
!  0.0 = replace with interpolated analysis
!
!    these values are set for analysis mode.
!
!   variables                  land                 sea
!   ---------------------------------------------------------
!   surface temperature        forecast             analysis
!   surface temperature        forecast             forecast (over sea ice)
!   albedo                     analysis             analysis
!   sea-ice                    analysis             analysis
!   snow                       analysis             forecast (over sea ice)
!   roughness                  analysis             forecast
!   plant resistance           analysis             analysis
!   soil wetness (layer)       weighted average     analysis
!   soil temperature           forecast             analysis
!   canopy waver content       forecast             forecast
!   convective cloud cover     forecast             forecast
!   convective cloud bottm     forecast             forecast
!   convective cloud top       forecast             forecast
!   vegetation cover           analysis             analysis
!   vegetation type            analysis             analysis
!   soil type                  analysis             analysis
!   sea-ice thickness          forecast             forecast
!   sea-ice concentration      analysis             analysis
!   vegetation cover min       analysis             analysis
!   vegetation cover max       analysis             analysis
!   max snow albedo            analysis             analysis
!   slope type                 analysis             analysis
!   liquid soil wetness        analysis-weighted    analysis
!   actual snow depth          analysis-weighted    analysis
!
!  note: if analysis file is not given, then time interpolated climatology
!        is used.  if analyiss file is given, it will be used as far as the
!        date and time matches.  if they do not match, it uses forecast.
!
!  critical percentage value for aborting bad points when lgchek=.true.
!
      logical lgchek
      data lgchek/.true./
      data critp1,critp2,critp3/80.,80.,25./
!
!     integer kpdalb(4), kpdalf(2)
!     data kpdalb/212,215,213,216/, kpdalf/214,217/
!     save kpdalb, kpdalf
!
!  mask orography and variance on gaussian grid
!
      real (kind=kind_io8) slmask(len),orog(len), orog_uf(len)
     &,                    orogd(len)
      real (kind=kind_io8) rla(len), rlo(len)
!
!  permanent/extremes
!
      character*500 fnglac,fnmxic
      real (kind=kind_io8), allocatable :: glacir(:),amxice(:),tsfcl0(:)
!
!     tsfcl0 is the climatological tsf at fh=0
!
!  climatology surface fields (last character 'c' or 'clm' indicate climatology)
!
      character*500 fntsfc,fnwetc,fnsnoc,fnzorc,fnalbc,fnaisc,
     &              fnplrc,fntg3c,fnscvc,fnsmcc,fnstcc,fnacnc,
     &              fnvegc,fnvetc,fnsotc
     &,             fnvmnc,fnvmxc,fnslpc,fnabsc, fnalbc2 
      real (kind=kind_io8) tsfclm(len), wetclm(len),   snoclm(len),
     &     zorclm(len), albclm(len,4), aisclm(len),
     &     tg3clm(len), acnclm(len),   cnpclm(len),
     &     cvclm (len), cvbclm(len),   cvtclm(len),
     &     scvclm(len), tsfcl2(len),   vegclm(len),
     &     vetclm(len), sotclm(len),   alfclm(len,2), sliclm(len),
     &     smcclm(len,lsoil), stcclm(len,lsoil)
     &,    sihclm(len), sicclm(len)
     &,    vmnclm(len), vmxclm(len), slpclm(len), absclm(len)
!
!  analyzed surface fields (last character 'a' or 'anl' indicate analysis)
!
      character*500 fntsfa,fnweta,fnsnoa,fnzora,fnalba,fnaisa,
     &             fnplra,fntg3a,fnscva,fnsmca,fnstca,fnacna,
     &             fnvega,fnveta,fnsota
     &,            fnvmna,fnvmxa,fnslpa,fnabsa       
!
      real (kind=kind_io8) tsfanl(len), wetanl(len),   snoanl(len),
     &     zoranl(len), albanl(len,4), aisanl(len),
     &     tg3anl(len), acnanl(len),   cnpanl(len),
     &     cvanl (len), cvbanl(len),   cvtanl(len),
     &     scvanl(len), tsfan2(len),   veganl(len),
     &     vetanl(len), sotanl(len),   alfanl(len,2), slianl(len),
     &     smcanl(len,lsoil), stcanl(len,lsoil)
     &,    sihanl(len), sicanl(len)
     &,    vmnanl(len), vmxanl(len), slpanl(len), absanl(len)
!
      real (kind=kind_io8) tsfan0(len) !  sea surface temperature analysis at ft=0.
!
!  predicted surface fields (last characters 'fcs' indicates forecast)
!
      real (kind=kind_io8) tsffcs(len), wetfcs(len),   snofcs(len),
     &     zorfcs(len), albfcs(len,4), aisfcs(len),
     &     tg3fcs(len), acnfcs(len),   cnpfcs(len),
     &     cvfcs (len), cvbfcs(len),   cvtfcs(len),
     &     slifcs(len), vegfcs(len),
     &     vetfcs(len), sotfcs(len),   alffcs(len,2),
     &     smcfcs(len,lsoil), stcfcs(len,lsoil)
     &,    sihfcs(len), sicfcs(len), sitfcs(len)
     &,    vmnfcs(len), vmxfcs(len), slpfcs(len), absfcs(len)
     &,    swdfcs(len), slcfcs(len,lsoil)
!
! ratio of sigma level 1 wind and 10m wind (diagnozed by model and not touched
! in this program).
!
      real (kind=kind_io8) f10m  (len)
      real (kind=kind_io8) fsmcl(25),fsmcs(25),fstcl(25),fstcs(25)
      real (kind=kind_io8) fcsmcl(25),fcsmcs(25),fcstcl(25),fcstcs(25)

!clu [+1l] add swratio (soil moisture liquid-to-total ratio)
      real (kind=kind_io8) swratio(len,lsoil)
!clu [+1l] add fixratio (option to adjust slc from smc)
      logical fixratio(lsoil)
!
      integer icsmcl(25), icsmcs(25), icstcl(25), icstcs(25)
!
      real (kind=kind_io8) csmcl(25), csmcs(25)
      real (kind=kind_io8) cstcl(25), cstcs(25)
!
      real (kind=kind_io8) slmskh(mdata)
      character*500 fnmskh
      integer kpd7, kpd9
!
      logical icefl1(len), icefl2(len)
!
!  input and output surface fields (bges) file names
!
!
!  sigma level 1 temperature for dead start
!
      real (kind=kind_io8) sig1t(len)
!
      character*32 label
!
!  = 1 ==> forecast is used
!  = 0 ==> analysis (or climatology) is used
!
!     output file  ... primary surface file for radiation and forecast
!
!       rec.  1    label
!       rec.  2    date record
!       rec.  3    tsf
!       rec.  4    soilm(two layers)              ----> 4 layers
!       rec.  5    snow
!       rec.  6    soilt(two layers)              ----> 4 layers
!       rec.  7    tg3
!       rec.  8    zor
!       rec.  9    cv
!       rec. 10    cvb
!       rec. 11    cvt
!       rec. 12    albedo (four types)
!       rec. 13    slimsk
!       rec. 14    vegetation cover
!       rec. 14    plantr                         -----> skip this record
!       rec. 15    f10m                           -----> canopy
!       rec. 16    canopy water content (cnpanl)  -----> f10m
!       rec. 17    vegetation type
!       rec. 18    soil type
!       rec. 19    zeneith angle dependent vegetation fraction (two types)
!       rec. 20    uustar
!       rec. 21    ffmm
!       rec. 22    ffhh
!cwu add sih & sic
!       rec. 23    sih(one category only)
!       rec. 24    sic
!clu [+8l] add prcp, flag, swd, slc, vmn, vmx, slp, abs
!       rec. 25    tprcp
!       rec. 26    srflag
!       rec. 27    swd
!       rec. 28    slc (4 layers)
!       rec. 29    vmn
!       rec. 30    vmx
!       rec. 31    slp
!       rec. 32    abs

!
!  debug only
!   ldebug=.true. creates bges files for climatology and analysis
!   lqcbgs=.true. quality controls input bges file before merging (should have been
!              qced in the forecast program)
!
      logical ldebug,lqcbgs
      logical lprnt
!
!  debug only
!
      character*500 fndclm,fndanl
!
      logical lanom

!
      namelist/namsfc/fnglac,fnmxic,
     &                fntsfc,fnwetc,fnsnoc,fnzorc,fnalbc,fnaisc,
     &                fnplrc,fntg3c,fnscvc,fnsmcc,fnstcc,fnacnc,
     &                fnvegc,fnvetc,fnsotc,fnalbc2,
     &                fnvmnc,fnvmxc,fnslpc,fnabsc,
     &                fntsfa,fnweta,fnsnoa,fnzora,fnalba,fnaisa,
     &                fnplra,fntg3a,fnscva,fnsmca,fnstca,fnacna,
     &                fnvega,fnveta,fnsota,
     &                fnvmna,fnvmxa,fnslpa,fnabsa,
     &                fnmskh,
     &                ldebug,lgchek,lqcbgs,critp1,critp2,critp3,
     &                fndclm,fndanl,
     &                lanom,
     &                ftsfl,ftsfs,falbl,falbs,faisl,faiss,fsnol,fsnos,
     &                fzorl,fzors,fplrl,fplrs,fsmcl,fsmcs,
     &                fstcl,fstcs,fvegl,fvegs,fvetl,fvets,fsotl,fsots,
     &                fctsfl,fctsfs,fcalbl,fcalbs,fcsnol,fcsnos,
     &                fczorl,fczors,fcplrl,fcplrs,fcsmcl,fcsmcs,
     &                fcstcl,fcstcs,fsalfl,fsalfs,fcalfl,flalfs,
     &                fsihl,fsicl,fsihs,fsics,aislim,sihnew,
     &                fvmnl,fvmns,fvmxl,fvmxs,fslpl,fslps,
     &                fabsl,fabss,
     &                ictsfl,ictsfs,icalbl,icalbs,icsnol,icsnos,
     &                iczorl,iczors,icplrl,icplrs,icsmcl,icsmcs,
     &                icstcl,icstcs,icalfl,icalfs,
     &                gausm,  deads, qcmsk, znlst,
     &                monclm, monanl, monfcs, monmer, mondif, igrdbg,
     &                blnmsk, bltmsk, landice
!
      data gausm/.true./,  deads/.false./, blnmsk/0.0/, bltmsk/90.0/
     &,    qcmsk/.false./, znlst/.false./, igrdbg/-1/
     &,    monclm/.false./, monanl/.false./, monfcs/.false./
     &,    monmer/.false./,  mondif/.false./,  landice/.true./
!
!  defaults file names
!
      data fnmskh/'global_slmask.t126.grb'/
      data fnalbc/'global_albedo4.1x1.grb'/
      data fnalbc2/'global_albedo4.1x1.grb'/
      data fntsfc/'global_sstclim.2x2.grb'/
      data fnsotc/'global_soiltype.1x1.grb'/
      data fnvegc/'global_vegfrac.1x1.grb'/
      data fnvetc/'global_vegtype.1x1.grb'/
      data fnglac/'global_glacier.2x2.grb'/
      data fnmxic/'global_maxice.2x2.grb'/
      data fnsnoc/'global_snoclim.1.875.grb'/
      data fnzorc/'global_zorclim.1x1.grb'/
      data fnaisc/'global_iceclim.2x2.grb'/
      data fntg3c/'global_tg3clim.2.6x1.5.grb'/
      data fnsmcc/'global_soilmcpc.1x1.grb'/
!clu [+4l] add fn()c for vmn, vmx, abs, slp
      data fnvmnc/'global_shdmin.0.144x0.144.grb'/
      data fnvmxc/'global_shdmax.0.144x0.144.grb'/
      data fnslpc/'global_slope.1x1.grb'/
      data fnabsc/'global_snoalb.1x1.grb'/
!
      data fnwetc/'        '/
      data fnplrc/'        '/
      data fnstcc/'        '/
      data fnscvc/'        '/
      data fnacnc/'        '/
!
      data fntsfa/'        '/
      data fnweta/'        '/
      data fnsnoa/'        '/
      data fnzora/'        '/
      data fnalba/'        '/
      data fnaisa/'        '/
      data fnplra/'        '/
      data fntg3a/'        '/
      data fnsmca/'        '/
      data fnstca/'        '/
      data fnscva/'        '/
      data fnacna/'        '/
      data fnvega/'        '/
      data fnveta/'        '/
      data fnsota/'        '/
!clu [+4l] add fn()a for vmn, vmx, abs, slp
      data fnvmna/'        '/
      data fnvmxa/'        '/
      data fnslpa/'        '/
      data fnabsa/'        '/
!
      data ldebug/.false./, lqcbgs/.true./
      data fndclm/'        '/
      data fndanl/'        '/
      data lanom/.false./
!
!  default relaxation time in hours to analysis or climatology
      data ftsfl/99999.0/,  ftsfs/0.0/
      data falbl/0.0/,      falbs/0.0/
      data falfl/0.0/,      falfs/0.0/
      data faisl/0.0/,      faiss/0.0/
      data fsnol/0.0/,      fsnos/99999.0/
      data fzorl/0.0/,      fzors/99999.0/
      data fplrl/0.0/,      fplrs/0.0/
      data fvetl/0.0/,      fvets/99999.0/
      data fsotl/0.0/,      fsots/99999.0/
      data fvegl/0.0/,      fvegs/99999.0/
!cwu [+4l] add f()l and f()s for sih, sic and aislim, sihlim
      data fsihl/99999.0/,  fsihs/99999.0/
!     data fsicl/99999.0/,  fsics/99999.0/
      data fsicl/0.0/,      fsics/0.0/
!  default ice concentration limit (50%), new ice thickness (20cm)
!cwu change ice concentration limit (15%) Jan 2015
!     data aislim/0.50/,    sihnew/0.2/
      data aislim/0.15/,    sihnew/0.2/
!clu [+4l] add f()l and f()s for vmn, vmx, abs, slp
      data fvmnl/0.0/,      fvmns/99999.0/
      data fvmxl/0.0/,      fvmxs/99999.0/
      data fslpl/0.0/,      fslps/99999.0/
      data fabsl/0.0/,      fabss/99999.0/
!  default relaxation time in hours to climatology if analysis missing
      data fctsfl/99999.0/, fctsfs/99999.0/
      data fcalbl/99999.0/, fcalbs/99999.0/
      data fcsnol/99999.0/, fcsnos/99999.0/
      data fczorl/99999.0/, fczors/99999.0/
      data fcplrl/99999.0/, fcplrs/99999.0/
!  default flag to apply climatological annual cycle
      data ictsfl/0/, ictsfs/1/
      data icalbl/1/, icalbs/1/
      data icalfl/1/, icalfs/1/
      data icsnol/0/, icsnos/0/
      data iczorl/1/, iczors/0/
      data icplrl/1/, icplrs/0/
!
      data ccnp/1.0/
      data ccv/1.0/,   ccvb/1.0/, ccvt/1.0/
!
      data ifp/0/
!
      save ifp,fnglac,fnmxic,
     &     fntsfc,fnwetc,fnsnoc,fnzorc,fnalbc,fnaisc,
     &     fnplrc,fntg3c,fnscvc,fnsmcc,fnstcc,fnacnc,fnvegc,
     &     fntsfa,fnweta,fnsnoa,fnzora,fnalba,fnaisa,
     &     fnplra,fntg3a,fnscva,fnsmca,fnstca,fnacna,fnvega,
     &     fnvetc,fnveta,
     &     fnsotc,fnsota,
!clu [+2l] add fn()c and fn()a for vmn, vmx, slp, abs
     &     fnvmnc,fnvmxc,fnabsc,fnslpc,
     &     fnvmna,fnvmxa,fnabsa,fnslpa,
     &     ldebug,lgchek,lqcbgs,critp1,critp2,critp3,
     &     fndclm,fndanl,
     &     lanom,
     &     ftsfl,ftsfs,falbl,falbs,faisl,faiss,fsnol,fsnos,
     &     fzorl,fzors,fplrl,fplrs,fsmcl,fsmcs,falfl,falfs,
     &     fstcl,fstcs,fvegl,fvegs,fvetl,fvets,fsotl,fsots,
     &     fctsfl,fctsfs,fcalbl,fcalbs,fcsnol,fcsnos,
     &     fczorl,fczors,fcplrl,fcplrs,fcsmcl,fcsmcs,
     &     fcstcl,fcstcs,fcalfl,fcalfs,
!cwu [+1l] add f()l and f()s for sih, sic and aislim, sihnew
     &     fsihl,fsihs,fsicl,fsics,aislim,sihnew,
!clu [+2l] add f()l and f()s for vmn, vmx, slp, abs
     &     fvmnl,fvmns,fvmxl,fvmxs,fslpl,fslps,
     &     fabsl,fabss,
     &     ictsfl,ictsfs,icalbl,icalbs,icsnol,icsnos,
     &     iczorl,iczors,icplrl,icplrs,icsmcl,icsmcs,
     &     icstcl,icstcs,icalfl,icalfs,
     &     gausm, deads, qcmsk,
     &     monclm, monanl, monfcs, monmer, mondif, igrdbg,
     &     grboro, grbmsk,
!
     &     ctsfl,  ctsfs,  calbl, calfl, calbs, calfs, csmcs,
     &     csnol,  csnos,  czorl, czors, cplrl, cplrs, cstcl,
     &     cstcs,  cvegl,  cvwgs, cvetl, cvets, csotl, csots,
     &     csmcl
!cwu [+1l] add c()l and c()s for sih, sic
     &,    csihl,  csihs,  csicl, csics
!clu [+2l] add c()l and c()s for vmn, vmx, slp, abs
     &,    cvmnl,  cvmns,  cvmxl, cvmxs, cslpl, cslps,
     &     cabsl,  cabss
     &,    imsk, jmsk, slmskh, blnmsk, bltmsk
     &,    glacir, amxice, tsfcl0
     &,    caisl, caiss, cvegs
!
      lprnt = .false.
      iprnt = 1
!     do i=1,len
!       if (ifp .eq. 0 .and. rla(i) .gt. 80.0) print *,' rla=',rla(i)
!    *,' rlo=',rlo(i)
!       tem1 = abs(rla(i) - 48.75)
!       tem2 = abs(rlo(i) - (-68.50))
!       if(tem1 .lt. 0.25 .and. tem2 .lt. 0.50) then
!         lprnt = .true.
!         iprnt = i
!         print *,' lprnt=',lprnt,' iprnt=',iprnt
!         print *,' rla(i)=',rla(i),' rlo(i)=',rlo(i)
!       endif
!     enddo
      if (ialb == 1) then
        kpdabs = kpdabs_1
        kpdalb = kpdalb_1
        alblmx = .99
        albsmx = .99
        alblmn = .01
        albsmn = .01
        abslmx = 1.0
        abssmx = 1.0
        abssmn = .01
        abslmn = .01
      else
        kpdabs = kpdabs_0
        kpdalb = kpdalb_0
        alblmx = .80
        albsmx = .80
        alblmn = .06
        albsmn = .06
        abslmx = .80
        abssmx = .80
        abslmn = .01
        abssmn = .01
      endif
      if(ifp.eq.0) then
        ifp = 1
        do k=1,lsoil
          fsmcl(k) = 99999.
          fsmcs(k) = 0.
          fstcl(k) = 99999.
          fstcs(k) = 0.
        enddo
!       print *,' in sfcsub nlunit=',nlunit,' me=',me,' ialb=',ialb
        rewind(nlunit)
        read (nlunit,namsfc)
!       write(6,namsfc)
!
        if (me .eq. 0) then
          print *,'ftsfl,falbl,faisl,fsnol,fzorl=',
     &    ftsfl,falbl,faisl,fsnol,fzorl
          print *,'fsmcl=',fsmcl(1:lsoil)
          print *,'fstcl=',fstcl(1:lsoil)
          print *,'ftsfs,falbs,faiss,fsnos,fzors=',
     &    ftsfs,falbs,faiss,fsnos,fzors
          print *,'fsmcs=',fsmcs(1:lsoil)
          print *,'fstcs=',fstcs(1:lsoil)
          print *,' aislim=',aislim,' sihnew=',sihnew
          print *,' isot=', isot,' ivegsrc=',ivegsrc
        endif

        if (ivegsrc == 2) then   ! sib
          veg_type_landice=13
        else
          veg_type_landice=15
        endif
        if (isot == 0) then
          soil_type_landice=9
        else
          soil_type_landice=16
        endif
!
        deltf = deltsfc / 24.0
!
        ctsfl=0.                       !...  tsfc over land
        if(ftsfl.ge.99999.) ctsfl=1.
        if((ftsfl.gt.0.).and.(ftsfl.lt.99999))  ctsfl=exp(-deltf/ftsfl)
!
        ctsfs=0.                       !...  tsfc over sea
        if(ftsfs.ge.99999.) ctsfs=1.
        if((ftsfs.gt.0.).and.(ftsfs.lt.99999))  ctsfs=exp(-deltf/ftsfs)
!
        do k=1,lsoil
          csmcl(k)=0.                  !...  soilm over land
          if(fsmcl(k).ge.99999.) csmcl(k)=1.
          if((fsmcl(k).gt.0.).and.(fsmcl(k).lt.99999))
     &                           csmcl(k)=exp(-deltf/fsmcl(k))
          csmcs(k)=0.                  !...  soilm over sea
          if(fsmcs(k).ge.99999.) csmcs(k)=1.
          if((fsmcs(k).gt.0.).and.(fsmcs(k).lt.99999))
     &                           csmcs(k)=exp(-deltf/fsmcs(k))
        enddo
!
        calbl=0.                       !...  albedo over land
        if(falbl.ge.99999.) calbl=1.
        if((falbl.gt.0.).and.(falbl.lt.99999))  calbl=exp(-deltf/falbl)
!
        calfl=0.                       !...  fraction field for albedo over land
        if(falfl.ge.99999.) calfl=1.
        if((falfl.gt.0.).and.(falfl.lt.99999))  calfl=exp(-deltf/falfl)
!
        calbs=0.                       !...  albedo over sea
        if(falbs.ge.99999.) calbs=1.
        if((falbs.gt.0.).and.(falbs.lt.99999))  calbs=exp(-deltf/falbs)
!
        calfs=0.                       !...  fraction field for albedo over sea
        if(falfs.ge.99999.) calfs=1.
        if((falfs.gt.0.).and.(falfs.lt.99999))  calfs=exp(-deltf/falfs)
!
        caisl=0.                       !...  sea ice over land
        if(faisl.ge.99999.) caisl=1.
        if((faisl.gt.0.).and.(faisl.lt.99999))  caisl=1.
!
        caiss=0.                       !...  sea ice over sea
        if(faiss.ge.99999.) caiss=1.
        if((faiss.gt.0.).and.(faiss.lt.99999))  caiss=1.
!
        csnol=0.                       !...  snow over land
        if(fsnol.ge.99999.) csnol=1.
        if((fsnol.gt.0.).and.(fsnol.lt.99999))  csnol=exp(-deltf/fsnol)
!       using the same way to bending snow as narr when fsnol is the negative value
!       the magnitude of fsnol is the thread to determine the lower and upper bound
!       of final swe
        if(fsnol.lt.0.)csnol=fsnol
!
        csnos=0.                       !...  snow over sea
        if(fsnos.ge.99999.) csnos=1.
        if((fsnos.gt.0.).and.(fsnos.lt.99999))  csnos=exp(-deltf/fsnos)
!
        czorl=0.                       !...  roughness length over land
        if(fzorl.ge.99999.) czorl=1.
        if((fzorl.gt.0.).and.(fzorl.lt.99999))  czorl=exp(-deltf/fzorl)
!
        czors=0.                       !...  roughness length over sea
        if(fzors.ge.99999.) czors=1.
        if((fzors.gt.0.).and.(fzors.lt.99999))  czors=exp(-deltf/fzors)
!
!       cplrl=0.                       !...  plant resistance over land
!       if(fplrl.ge.99999.) cplrl=1.
!       if((fplrl.gt.0.).and.(fplrl.lt.99999))  cplrl=exp(-deltf/fplrl)
!
!       cplrs=0.                       !...  plant resistance over sea
!       if(fplrs.ge.99999.) cplrs=1.
!       if((fplrs.gt.0.).and.(fplrs.lt.99999))  cplrs=exp(-deltf/fplrs)
!
        do k=1,lsoil
           cstcl(k)=0.                 !...  soilt over land
           if(fstcl(k).ge.99999.) cstcl(k)=1.
           if((fstcl(k).gt.0.).and.(fstcl(k).lt.99999))
     &                            cstcl(k)=exp(-deltf/fstcl(k))
          cstcs(k)=0.                  !...  soilt over sea
          if(fstcs(k).ge.99999.) cstcs(k)=1.
          if((fstcs(k).gt.0.).and.(fstcs(k).lt.99999))
     &                           cstcs(k)=exp(-deltf/fstcs(k))
        enddo
!
        cvegl=0.                       !...  vegetation fraction over land
        if(fvegl.ge.99999.) cvegl=1.
        if((fvegl.gt.0.).and.(fvegl.lt.99999))  cvegl=exp(-deltf/fvegl)
!
        cvegs=0.                       !...  vegetation fraction over sea
        if(fvegs.ge.99999.) cvegs=1.
        if((fvegs.gt.0.).and.(fvegs.lt.99999))  cvegs=exp(-deltf/fvegs)
!
        cvetl=0.                       !...  vegetation type over land
        if(fvetl.ge.99999.) cvetl=1.
        if((fvetl.gt.0.).and.(fvetl.lt.99999))  cvetl=exp(-deltf/fvetl)
!
        cvets=0.                       !...  vegetation type over sea
        if(fvets.ge.99999.) cvets=1.
        if((fvets.gt.0.).and.(fvets.lt.99999))  cvets=exp(-deltf/fvets)
!
        csotl=0.                       !...  soil type over land
        if(fsotl.ge.99999.) csotl=1.
        if((fsotl.gt.0.).and.(fsotl.lt.99999))  csotl=exp(-deltf/fsotl)
!
        csots=0.                       !...  soil type over sea
        if(fsots.ge.99999.) csots=1.
        if((fsots.gt.0.).and.(fsots.lt.99999))  csots=exp(-deltf/fsots)

!cwu [+16l]---------------------------------------------------------------
!
        csihl=0.                       !...  sea ice thickness over land
        if(fsihl.ge.99999.) csihl=1.
        if((fsihl.gt.0.).and.(fsihl.lt.99999))  csihl=exp(-deltf/fsihl)
!
        csihs=0.                       !...  sea ice thickness over sea
        if(fsihs.ge.99999.) csihs=1.
        if((fsihs.gt.0.).and.(fsihs.lt.99999))  csihs=exp(-deltf/fsihs)
!
        csicl=0.                       !...  sea ice concentration over land
        if(fsicl.ge.99999.) csicl=1.
        if((fsicl.gt.0.).and.(fsicl.lt.99999))  csicl=exp(-deltf/fsicl)
!
        csics=0.                       !...  sea ice concentration over sea
        if(fsics.ge.99999.) csics=1.
        if((fsics.gt.0.).and.(fsics.lt.99999))  csics=exp(-deltf/fsics)

!clu [+32l]---------------------------------------------------------------
!
        cvmnl=0.                       !...  min veg cover over land
        if(fvmnl.ge.99999.) cvmnl=1.
        if((fvmnl.gt.0.).and.(fvmnl.lt.99999))  cvmnl=exp(-deltf/fvmnl)
!
        cvmns=0.                       !...  min veg cover over sea
        if(fvmns.ge.99999.) cvmns=1.
        if((fvmns.gt.0.).and.(fvmns.lt.99999))  cvmns=exp(-deltf/fvmns)
!
        cvmxl=0.                       !...  max veg cover over land
        if(fvmxl.ge.99999.) cvmxl=1.
        if((fvmxl.gt.0.).and.(fvmxl.lt.99999))  cvmxl=exp(-deltf/fvmxl)
!
        cvmxs=0.                       !...  max veg cover over sea
        if(fvmxs.ge.99999.) cvmxs=1.
        if((fvmxs.gt.0.).and.(fvmxs.lt.99999))  cvmxs=exp(-deltf/fvmxs)
!
        cslpl=0.                       !... slope type over land
        if(fslpl.ge.99999.) cslpl=1.
        if((fslpl.gt.0.).and.(fslpl.lt.99999))  cslpl=exp(-deltf/fslpl)
!
        cslps=0.                       !...  slope type over sea
        if(fslps.ge.99999.) cslps=1.
        if((fslps.gt.0.).and.(fslps.lt.99999))  cslps=exp(-deltf/fslps)
!
        cabsl=0.                       !... snow albedo over land
        if(fabsl.ge.99999.) cabsl=1.
        if((fabsl.gt.0.).and.(fabsl.lt.99999))  cabsl=exp(-deltf/fabsl)
!
        cabss=0.                       !... snow albedo over sea
        if(fabss.ge.99999.) cabss=1.
        if((fabss.gt.0.).and.(fabss.lt.99999))  cabss=exp(-deltf/fabss)
!clu ----------------------------------------------------------------------
!
!     read a high resolution mask field for use in grib interpolation
!
        call hmskrd(lugb,imsk,jmsk,fnmskh,
     &              kpdmsk,slmskh,gausm,blnmsk,bltmsk,me)
!       if (qcmsk) call qcmask(slmskh,sllnd,slsea,imsk,jmsk,rla,rlo)
!
        if (me .eq. 0) then
          write(6,*) ' '
          write(6,*) ' lugb=',lugb,' len=',len, ' lsoil=',lsoil
          write(6,*) 'iy=',iy,' im=',im,' id=',id,' ih=',ih,' fh=',fh
     &,            ' sig1t(1)=',sig1t(1)
     &,            ' gausm=',gausm,' blnmsk=',blnmsk,' bltmsk=',bltmsk
          write(6,*) ' '
        endif
!
!  reading permanent/extreme features (glacier points and maximum ice extent)
!
        allocate (tsfcl0(len))
        allocate (glacir(len))
        allocate (amxice(len))
!
!  read glacier
!
        kpd9 = -1
        kpd7 = -1
        call fixrdc(lugb,fnglac,kpdgla,kpd7,kpd9,slmask,
     &              glacir,len,iret
     &,             imsk, jmsk, slmskh, gausm, blnmsk, bltmsk
     &,             rla, rlo, me)
!     znnt=1.
!     call nntprt(glacir,len,znnt)
!
!  read maximum ice extent
!
        kpd7 = -1
        call fixrdc(lugb,fnmxic,kpdmxi,kpd7,kpd9,slmask,
     &              amxice,len,iret
     &,             imsk, jmsk, slmskh, gausm, blnmsk, bltmsk
     &,             rla, rlo, me)
!     znnt=1.
!     call nntprt(amxice,len,znnt)
!
        crit=0.5
        call rof01(glacir,len,'ge',crit)
        call rof01(amxice,len,'ge',crit)
!
!  quality control max ice limit based on glacier points
!
        call qcmxice(glacir,amxice,len,me)
!
      endif                       ! first time loop finished
!
      do i=1,len
        sliclm(i) = 1.
        snoclm(i) = 0.
        icefl1(i) = .true.
      enddo
!     if(lprnt) print *,' tsffcsin=',tsffcs(iprnt)
!
!  read climatology fields
!
      if (me .eq. 0) then
        write(6,*) '=============='
        write(6,*) 'climatology'
        write(6,*) '=============='
      endif
!
      percrit=critp1
!
      call clima(lugb,iy,im,id,ih,fh,len,lsoil,slmask,
     &           fntsfc,fnwetc,fnsnoc,fnzorc,fnalbc,fnaisc,
     &           fntg3c,fnscvc,fnsmcc,fnstcc,fnacnc,fnvegc,
     &           fnvetc,fnsotc,
     &           fnvmnc,fnvmxc,fnslpc,fnabsc,
     &           tsfclm,tsfcl2,wetclm,snoclm,zorclm,albclm,aisclm,
     &           tg3clm,cvclm ,cvbclm,cvtclm,
     &           cnpclm,smcclm,stcclm,sliclm,scvclm,acnclm,vegclm,
     &           vetclm,sotclm,alfclm,
     &           vmnclm,vmxclm,slpclm,absclm,
     &           kpdtsf,kpdwet,kpdsno,kpdzor,kpdalb,kpdais,
     &           kpdtg3,kpdscv,kpdacn,kpdsmc,kpdstc,kpdveg,
     &           kpdvet,kpdsot,kpdalf,tsfcl0,
     &           kpdvmn,kpdvmx,kpdslp,kpdabs,
     &           deltsfc, lanom
     &,          imsk, jmsk, slmskh, rla, rlo, gausm, blnmsk, bltmsk,me
     &,          lprnt, iprnt, fnalbc2, ialb)
!     if(lprnt) print *,'tsfclm=',tsfclm(iprnt),' tsfcl2=',tsfcl2(iprnt)
!
!  scale surface roughness and albedo to model required units
!
      zsca=100.
      call scale(zorclm,len,zsca)
      zsca=0.01
      call scale(albclm,len,zsca)
      call scale(albclm(1,2),len,zsca)
      call scale(albclm(1,3),len,zsca)
      call scale(albclm(1,4),len,zsca)
      call scale(alfclm,len,zsca)
      call scale(alfclm(1,2),len,zsca)
!clu [+4l] scale vmn, vmx, abs from percent to fraction
      zsca=0.01
      call scale(vmnclm,len,zsca)
      call scale(vmxclm,len,zsca)
      call scale(absclm,len,zsca)

!
!  set albedo over ocean to albomx
!
      call albocn(albclm,slmask,albomx,len)
!
!  make sure vegetation type and soil type are non zero over land
!
      call landtyp(vetclm,sotclm,slpclm,slmask,len)
!
!cwu [-1l/+1l]
!* ice concentration or ice mask (only ice mask used in the model now)
!  ice concentration and ice mask (both are used in the model now)
!
      if(fnaisc(1:8).ne.'        ') then
!cwu [+5l/-1l] update sihclm, sicclm
        do i=1,len
         sihclm(i) = 3.0*aisclm(i)
         sicclm(i) = aisclm(i)
          if(slmask(i).eq.0..and.glacir(i).eq.1..and.
     &      sicclm(i).ne.1.) then
            sicclm(i) = sicimx
            sihfcs(i) = glacir_hice
          endif
        enddo
        crit=aislim
!*      crit=0.5
        call rof01(aisclm,len,'ge',crit)
      elseif(fnacnc(1:8).ne.'        ') then
!cwu [+4l] update sihclm, sicclm
        do i=1,len
         sihclm(i) = 3.0*acnclm(i)
         sicclm(i) = acnclm(i)
          if(slmask(i).eq.0..and.glacir(i).eq.1..and.
     &      sicclm(i).ne.1.) then
            sicclm(i) = sicimx
            sihfcs(i) = glacir_hice
          endif
        enddo
        call rof01(acnclm,len,'ge',aislim)
        do i=1,len
         aisclm(i) = acnclm(i)
        enddo
      endif
!
!  quality control of sea ice mask
!
      call qcsice(aisclm,glacir,amxice,aicice,aicsea,sllnd,slmask,
     &            rla,rlo,len,me)
!
!  set ocean/land/sea-ice mask
!
      call setlsi(slmask,aisclm,len,aicice,sliclm)
!     if(lprnt) print *,' aisclm=',aisclm(iprnt),' sliclm='
!    *,sliclm(iprnt),' slmask=',slmask(iprnt)
!
!     write(6,*) 'sliclm'
!     znnt=1.
!     call nntprt(sliclm,len,znnt)
!
!  quality control of snow
!
      call qcsnow(snoclm,slmask,aisclm,glacir,len,snosmx,landice,me)
!
      call setzro(snoclm,epssno,len)
!
!  snow cover handling (we assume climatological snow depth is available)
!  quality control of snow depth (note that snow should be corrected first
!  because it influences tsf
!
      kqcm=1
      call qcmxmn('snow    ',snoclm,sliclm,snoclm,icefl1,
     &            snolmx,snolmn,snoomx,snoomn,snoimx,snoimn,
     &            snojmx,snojmn,snosmx,snosmn,epssno,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!     write(6,*) 'snoclm'
!     znnt=1.
!     call nntprt(snoclm,len,znnt)
!
!  get snow cover from snow depth array
!
      if(fnscvc(1:8).eq.'        ') then
        call getscv(snoclm,scvclm,len)
      endif
!
!  set tsfc over snow to tsfsmx if greater
!
      call snosfc(snoclm,tsfclm,tsfsmx,len,me)
!     call snosfc(snoclm,tsfcl2,tsfsmx,len)

!
!  quality control
!
      do i=1,len
        icefl2(i) = sicclm(i) .gt. 0.99999
      enddo
      kqcm=1
      call qcmxmn('tsfc    ',tsfclm,sliclm,snoclm,icefl2,
     &            tsflmx,tsflmn,tsfomx,tsfomn,tsfimx,tsfimn,
     &            tsfjmx,tsfjmn,tsfsmx,tsfsmn,epstsf,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('tsf2    ',tsfcl2,sliclm,snoclm,icefl2,
     &            tsflmx,tsflmn,tsfomx,tsfomn,tsfimx,tsfimn,
     &            tsfjmx,tsfjmn,tsfsmx,tsfsmn,epstsf,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      do kk = 1, 4
      call qcmxmn('albc    ',albclm(1,kk),sliclm,snoclm,icefl1,
     &            alblmx,alblmn,albomx,albomn,albimx,albimn,
     &            albjmx,albjmn,albsmx,albsmn,epsalb,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      enddo
      if(fnwetc(1:8).ne.'        ') then
        call qcmxmn('wetc    ',wetclm,sliclm,snoclm,icefl1,
     &              wetlmx,wetlmn,wetomx,wetomn,wetimx,wetimn,
     &              wetjmx,wetjmn,wetsmx,wetsmn,epswet,
     &              rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
      call qcmxmn('zorc    ',zorclm,sliclm,snoclm,icefl1,
     &            zorlmx,zorlmn,zoromx,zoromn,zorimx,zorimn,
     &            zorjmx,zorjmn,zorsmx,zorsmn,epszor,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!     if(fnplrc(1:8).ne.'        ') then
!     call qcmxmn('plntc   ',plrclm,sliclm,snoclm,icefl1,
!    &            plrlmx,plrlmn,plromx,plromn,plrimx,plrimn,
!    &            plrjmx,plrjmn,plrsmx,plrsmn,epsplr,
!    &            rla,rlo,len,kqcm,percrit,lgchek,me)
!     endif
      call qcmxmn('tg3c    ',tg3clm,sliclm,snoclm,icefl1,
     &            tg3lmx,tg3lmn,tg3omx,tg3omn,tg3imx,tg3imn,
     &            tg3jmx,tg3jmn,tg3smx,tg3smn,epstg3,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!
!  get soil temp and moisture (after all the qcs are completed)
!
      if(fnsmcc(1:8).eq.'        ') then
        call getsmc(wetclm,len,lsoil,smcclm,me)
      endif
      call qcmxmn('smc1c   ',smcclm(1,1),sliclm,snoclm,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('smc2c   ',smcclm(1,2),sliclm,snoclm,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+8l] add smcclm(3:4)
      if(lsoil.gt.2) then
      call qcmxmn('smc3c   ',smcclm(1,3),sliclm,snoclm,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('smc4c   ',smcclm(1,4),sliclm,snoclm,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
      if(fnstcc(1:8).eq.'        ') then
        call getstc(tsfclm,tg3clm,sliclm,len,lsoil,stcclm,tsfimx)
      endif
      call qcmxmn('stc1c   ',stcclm(1,1),sliclm,snoclm,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('stc2c   ',stcclm(1,2),sliclm,snoclm,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+8l] add stcclm(3:4)
      if(lsoil.gt.2) then
      call qcmxmn('stc3c   ',stcclm(1,3),sliclm,snoclm,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('stc4c   ',stcclm(1,4),sliclm,snoclm,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
      call qcmxmn('vegc    ',vegclm,sliclm,snoclm,icefl1,
     &            veglmx,veglmn,vegomx,vegomn,vegimx,vegimn,
     &            vegjmx,vegjmn,vegsmx,vegsmn,epsveg,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('vetc    ',vetclm,sliclm,snoclm,icefl1,
     &            vetlmx,vetlmn,vetomx,vetomn,vetimx,vetimn,
     &            vetjmx,vetjmn,vetsmx,vetsmn,epsvet,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('sotc    ',sotclm,sliclm,snoclm,icefl1,
     &            sotlmx,sotlmn,sotomx,sotomn,sotimx,sotimn,
     &            sotjmx,sotjmn,sotsmx,sotsmn,epssot,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!cwu [+8l] ---------------------------------------------------------------
      call qcmxmn('sihc    ',sihclm,sliclm,snoclm,icefl1,
     &            sihlmx,sihlmn,sihomx,sihomn,sihimx,sihimn,
     &            sihjmx,sihjmn,sihsmx,sihsmn,epssih,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('sicc    ',sicclm,sliclm,snoclm,icefl1,
     &            siclmx,siclmn,sicomx,sicomn,sicimx,sicimn,
     &            sicjmx,sicjmn,sicsmx,sicsmn,epssic,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+16l] ---------------------------------------------------------------
      call qcmxmn('vmnc    ',vmnclm,sliclm,snoclm,icefl1,
     &            vmnlmx,vmnlmn,vmnomx,vmnomn,vmnimx,vmnimn,
     &            vmnjmx,vmnjmn,vmnsmx,vmnsmn,epsvmn,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('vmxc    ',vmxclm,sliclm,snoclm,icefl1,
     &            vmxlmx,vmxlmn,vmxomx,vmxomn,vmximx,vmximn,
     &            vmxjmx,vmxjmn,vmxsmx,vmxsmn,epsvmx,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('slpc    ',slpclm,sliclm,snoclm,icefl1,
     &            slplmx,slplmn,slpomx,slpomn,slpimx,slpimn,
     &            slpjmx,slpjmn,slpsmx,slpsmn,epsslp,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('absc    ',absclm,sliclm,snoclm,icefl1,
     &            abslmx,abslmn,absomx,absomn,absimx,absimn,
     &            absjmx,absjmn,abssmx,abssmn,epsabs,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu ----------------------------------------------------------------------
!
!  monitoring prints
!
      if (monclm) then
       if (me .eq. 0) then
        print *,' '
        print *,'monitor of time and space interpolated climatology'
        print *,' '
!       call count(sliclm,snoclm,len)
        print *,' '
        call monitr('tsfclm',tsfclm,sliclm,snoclm,len)
        call monitr('albclm',albclm(1,1),sliclm,snoclm,len)
        call monitr('albclm',albclm(1,2),sliclm,snoclm,len)
        call monitr('albclm',albclm(1,3),sliclm,snoclm,len)
        call monitr('albclm',albclm(1,4),sliclm,snoclm,len)
        call monitr('aisclm',aisclm,sliclm,snoclm,len)
        call monitr('snoclm',snoclm,sliclm,snoclm,len)
        call monitr('scvclm',scvclm,sliclm,snoclm,len)
        call monitr('smcclm1',smcclm(1,1),sliclm,snoclm,len)
        call monitr('smcclm2',smcclm(1,2),sliclm,snoclm,len)
        call monitr('stcclm1',stcclm(1,1),sliclm,snoclm,len)
        call monitr('stcclm2',stcclm(1,2),sliclm,snoclm,len)
!clu [+4l] add smcclm(3:4) and stcclm(3:4)
        if(lsoil.gt.2) then
        call monitr('smcclm3',smcclm(1,3),sliclm,snoclm,len)
        call monitr('smcclm4',smcclm(1,4),sliclm,snoclm,len)
        call monitr('stcclm3',stcclm(1,3),sliclm,snoclm,len)
        call monitr('stcclm4',stcclm(1,4),sliclm,snoclm,len)
        endif
        call monitr('tg3clm',tg3clm,sliclm,snoclm,len)
        call monitr('zorclm',zorclm,sliclm,snoclm,len)
!       if (gaus) then
          call monitr('cvaclm',cvclm ,sliclm,snoclm,len)
          call monitr('cvbclm',cvbclm,sliclm,snoclm,len)
          call monitr('cvtclm',cvtclm,sliclm,snoclm,len)
!       endif
        call monitr('sliclm',sliclm,sliclm,snoclm,len)
!       call monitr('plrclm',plrclm,sliclm,snoclm,len)
        call monitr('orog  ',orog  ,sliclm,snoclm,len)
        call monitr('vegclm',vegclm,sliclm,snoclm,len)
        call monitr('vetclm',vetclm,sliclm,snoclm,len)
        call monitr('sotclm',sotclm,sliclm,snoclm,len)
!cwu [+2l] add sih, sic
        call monitr('sihclm',sihclm,sliclm,snoclm,len)
        call monitr('sicclm',sicclm,sliclm,snoclm,len)
!clu [+4l] add vmn, vmx, slp, abs
        call monitr('vmnclm',vmnclm,sliclm,snoclm,len)
        call monitr('vmxclm',vmxclm,sliclm,snoclm,len)
        call monitr('slpclm',slpclm,sliclm,snoclm,len)
        call monitr('absclm',absclm,sliclm,snoclm,len)
       endif
      endif
!
!
      if (me .eq. 0) then
        write(6,*) '=============='
        write(6,*) '   analysis'
        write(6,*) '=============='
      endif
!
!  fill in analysis array with climatology before reading analysis.
!
      call filanl(tsfanl,tsfan2,wetanl,snoanl,zoranl,albanl,aisanl,
     &            tg3anl,cvanl ,cvbanl,cvtanl,
     &            cnpanl,smcanl,stcanl,slianl,scvanl,veganl,
     &            vetanl,sotanl,alfanl,
     &            sihanl,sicanl,
     &            vmnanl,vmxanl,slpanl,absanl, 
     &            tsfclm,tsfcl2,wetclm,snoclm,zorclm,albclm,aisclm,
     &            tg3clm,cvclm ,cvbclm,cvtclm,
     &            cnpclm,smcclm,stcclm,sliclm,scvclm,vegclm,
     &            vetclm,sotclm,alfclm,
     &            sihclm,sicclm,
     &            vmnclm,vmxclm,slpclm,absclm,      
     &            len,lsoil)
!
!  reverse scaling to match with grib analysis input
!
      zsca=0.01
      call scale(zoranl,len, zsca)
      zsca=100.
      call scale(albanl,len,zsca)
      call scale(albanl(1,2),len,zsca)
      call scale(albanl(1,3),len,zsca)
      call scale(albanl(1,4),len,zsca)
      call scale(alfanl,len,zsca)
      call scale(alfanl(1,2),len,zsca)
!clu [+4l] reverse scale for vmn, vmx, abs
      zsca=100.
      call scale(vmnanl,len,zsca)
      call scale(vmxanl,len,zsca)
      call scale(absanl,len,zsca)
!
      percrit=critp2
!
!  read analysis fields
!
      call analy(lugb,iy,im,id,ih,fh,len,lsoil,slmask,
     &           fntsfa,fnweta,fnsnoa,fnzora,fnalba,fnaisa,
     &           fntg3a,fnscva,fnsmca,fnstca,fnacna,fnvega,
     &           fnveta,fnsota,
     &           fnvmna,fnvmxa,fnslpa,fnabsa,      
     &           tsfanl,wetanl,snoanl,zoranl,albanl,aisanl,
     &           tg3anl,cvanl ,cvbanl,cvtanl,
     &           smcanl,stcanl,slianl,scvanl,acnanl,veganl,
     &           vetanl,sotanl,alfanl,tsfan0,
     &           vmnanl,vmxanl,slpanl,absanl,      
     &           kpdtsf,kpdwet,kpdsno,kpdsnd,kpdzor,kpdalb,kpdais,
     &           kpdtg3,kpdscv,kpdacn,kpdsmc,kpdstc,kpdveg,
     &           kpdvet,kpdsot,kpdalf,
     &           kpdvmn,kpdvmx,kpdslp,kpdabs,      
     &           irttsf,irtwet,irtsno,irtzor,irtalb,irtais,
     &           irttg3,irtscv,irtacn,irtsmc,irtstc,irtveg,
     &           irtvet,irtsot,irtalf
     &,          irtvmn,irtvmx,irtslp,irtabs, 
     &           imsk, jmsk, slmskh, rla, rlo, gausm, blnmsk, bltmsk,
     &           me)
!     if(lprnt) print *,' tsfanl=',tsfanl(iprnt)
!
!  scale zor and alb to match forecast model units
!
      zsca=100.
      call scale(zoranl,len, zsca)
      zsca=0.01
      call scale(albanl,len,zsca)
      call scale(albanl(1,2),len,zsca)
      call scale(albanl(1,3),len,zsca)
      call scale(albanl(1,4),len,zsca)
      call scale(alfanl,len,zsca)
      call scale(alfanl(1,2),len,zsca)
!clu [+4] scale vmn, vmx, abs from percent to fraction
      zsca=0.01
      call scale(vmnanl,len,zsca)
      call scale(vmxanl,len,zsca)
      call scale(absanl,len,zsca)
!
!  interpolate climatology but fixing initial anomaly
!
      if(fh.gt.0.0.and.fntsfa(1:8).ne.'        '.and.lanom) then
        call anomint(tsfan0,tsfclm,tsfcl0,tsfanl,len)
      endif
!
!    if the tsfanl is at sea level, then bring it to the surface using
!    unfiltered orography (for lakes).  if the analysis is at lake surface
!    as in the nst model, then this call should be removed - moorthi 09/23/2011
!
        if (use_ufo .and. .not. nst_anl) then
          ztsfc = 0.0
          call tsfcor(tsfanl,orog_uf,slmask,ztsfc,len,rlapse)
        endif
!
!  ice concentration or ice mask (only ice mask used in the model now)
!
      if(fnaisa(1:8).ne.'        ') then
!cwu [+5l/-1l] update sihanl, sicanl
        do i=1,len
         sihanl(i) = 3.0*aisanl(i)
         sicanl(i) = aisanl(i)
          if(slmask(i).eq.0..and.glacir(i).eq.1..and.
     &      sicanl(i).ne.1.) then
            sicanl(i) = sicimx
            sihfcs(i) = glacir_hice
          endif
        enddo
        crit=aislim
!*      crit=0.5
        call rof01(aisanl,len,'ge',crit)
      elseif(fnacna(1:8).ne.'        ') then
!cwu [+17l] update sihanl, sicanl
        do i=1,len
          sihanl(i) = 3.0*acnanl(i)
          sicanl(i) = acnanl(i)
          if(slmask(i).eq.0..and.glacir(i).eq.1..and.
     &     sicanl(i).ne.1.) then
            sicanl(i) = sicimx
            sihfcs(i) = glacir_hice
          endif
        enddo
        crit=aislim
        do i=1,len
          if((slianl(i).eq.0.).and.(sicanl(i).ge.crit)) then
            slianl(i)=2.
!           print *,'cycle - new ice form: fice=',sicanl(i)
          else if((slianl(i).ge.2.).and.(sicanl(i).lt.crit)) then
            slianl(i)=0.
!           print *,'cycle - ice free: fice=',sicanl(i)
          else if((slianl(i).eq.1.).and.(sicanl(i).ge.sicimn)) then
!           print *,'cycle - land covered by sea-ice: fice=',sicanl(i)
            sicanl(i)=0.
          endif
        enddo
!       znnt=10.
!       call nntprt(acnanl,len,znnt)
!     if(lprnt) print *,' acnanl=',acnanl(iprnt)
!       do i=1,len
!         if (acnanl(i) .gt. 0.3 .and. aisclm(i) .eq. 1.0
!    &     .and. aisfcs(i) .ge. 0.75)   acnanl(i) = aislim
!       enddo
!     if(lprnt) print *,' acnanl=',acnanl(iprnt)
        call rof01(acnanl,len,'ge',aislim)
        do i=1,len
          aisanl(i)=acnanl(i)
        enddo
      endif
!     if(lprnt) print *,' aisanl1=',aisanl(iprnt),' glacir='
!    &,glacir(iprnt),' slmask=',slmask(iprnt)
!
      call qcsice(aisanl,glacir,amxice,aicice,aicsea,sllnd,slmask,
     &            rla,rlo,len,me)
!
!  set ocean/land/sea-ice mask
!
      call setlsi(slmask,aisanl,len,aicice,slianl)
!     if(lprnt) print *,' aisanl=',aisanl(iprnt),' slianl='
!    *,slianl(iprnt),' slmask=',slmask(iprnt)
!
!
      do k=1,lsoil
        do i=1,len
          if (slianl(i) .eq. 0) then
            smcanl(i,k) = smcomx
            stcanl(i,k) = tsfanl(i)
          endif
        enddo
      enddo

!     write(6,*) 'slianl'
!     znnt=1.
!     call nntprt(slianl,len,znnt)
!cwu [+8l]----------------------------------------------------------------------
      call qcmxmn('siha    ',sihanl,slianl,snoanl,icefl1,
     &            sihlmx,sihlmn,sihomx,sihomn,sihimx,sihimn,
     &            sihjmx,sihjmn,sihsmx,sihsmn,epssih,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('sica    ',sicanl,slianl,snoanl,icefl1,
     &            siclmx,siclmn,sicomx,sicomn,sicimx,sicimn,
     &            sicjmx,sicjmn,sicsmx,sicsmn,epssic,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!
!  set albedo over ocean to albomx
!
      call albocn(albanl,slmask,albomx,len)
!
!  quality control of snow and sea-ice
!    process snow depth or snow cover
!
      if(fnsnoa(1:8).ne.'        ') then
        call setzro(snoanl,epssno,len)
        call qcsnow(snoanl,slmask,aisanl,glacir,len,ten,landice,me)
        if (.not.landice) then
          call snodpth2(glacir,snosmx,snoanl, len, me)
        endif
        kqcm=1
        call snosfc(snoanl,tsfanl,tsfsmx,len,me)
        call qcmxmn('snoa    ',snoanl,slianl,snoanl,icefl1,
     &              snolmx,snolmn,snoomx,snoomn,snoimx,snoimn,
     &              snojmx,snojmn,snosmx,snosmn,epssno,
     &              rla,rlo,len,kqcm,percrit,lgchek,me)
        call getscv(snoanl,scvanl,len)
        call qcmxmn('sncva   ',scvanl,slianl,snoanl,icefl1,
     &              scvlmx,scvlmn,scvomx,scvomn,scvimx,scvimn,
     &              scvjmx,scvjmn,scvsmx,scvsmn,epsscv,
     &              rla,rlo,len,kqcm,percrit,lgchek,me)
      else
        crit=0.5
        call rof01(scvanl,len,'ge',crit)
        call qcsnow(scvanl,slmask,aisanl,glacir,len,one,landice,me)
        call qcmxmn('sncva   ',scvanl,slianl,scvanl,icefl1,
     &              scvlmx,scvlmn,scvomx,scvomn,scvimx,scvimn,
     &              scvjmx,scvjmn,scvsmx,scvsmn,epsscv,
     &              rla,rlo,len,kqcm,percrit,lgchek,me)
        call snodpth(scvanl,slianl,tsfanl,snoclm,
     &               glacir,snwmax,snwmin,landice,len,snoanl,me)
        call qcsnow(scvanl,slmask,aisanl,glacir,len,snosmx,landice,me)
        call snosfc(snoanl,tsfanl,tsfsmx,len,me)
        call qcmxmn('snowa   ',snoanl,slianl,snoanl,icefl1,
     &              snolmx,snolmn,snoomx,snoomn,snoimx,snoimn,
     &              snojmx,snojmn,snosmx,snosmn,epssno,
     &              rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
!
      do i=1,len
        icefl2(i) = sicanl(i) .gt. 0.99999
      enddo
      call qcmxmn('tsfa    ',tsfanl,slianl,snoanl,icefl2,
     &            tsflmx,tsflmn,tsfomx,tsfomn,tsfimx,tsfimn,
     &            tsfjmx,tsfjmn,tsfsmx,tsfsmn,epstsf,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      do kk = 1, 4
      call qcmxmn('alba    ',albanl(1,kk),slianl,snoanl,icefl1,
     &            alblmx,alblmn,albomx,albomn,albimx,albimn,
     &            albjmx,albjmn,albsmx,albsmn,epsalb,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      enddo
      if(fnwetc(1:8).ne.'        ' .or. fnweta(1:8).ne.'        ' ) then
      call qcmxmn('weta    ',wetanl,slianl,snoanl,icefl1,
     &            wetlmx,wetlmn,wetomx,wetomn,wetimx,wetimn,
     &            wetjmx,wetjmn,wetsmx,wetsmn,epswet,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
      call qcmxmn('zora    ',zoranl,slianl,snoanl,icefl1,
     &            zorlmx,zorlmn,zoromx,zoromn,zorimx,zorimn,
     &            zorjmx,zorjmn,zorsmx,zorsmn,epszor,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!     if(fnplrc(1:8).ne.'        ' .or. fnplra(1:8).ne.'        ' ) then
!     call qcmxmn('plna    ',plranl,slianl,snoanl,icefl1,
!    &            plrlmx,plrlmn,plromx,plromn,plrimx,plrimn,
!    &            plrjmx,plrjmn,plrsmx,plrsmn,epsplr,
!    &            rla,rlo,len,kqcm,percrit,lgchek,me)
!     endif
      call qcmxmn('tg3a    ',tg3anl,slianl,snoanl,icefl1,
     &            tg3lmx,tg3lmn,tg3omx,tg3omn,tg3imx,tg3imn,
     &            tg3jmx,tg3jmn,tg3smx,tg3smn,epstg3,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!
!  get soil temp and moisture
!
      if(fnsmca(1:8).eq.'        ' .and. fnsmcc(1:8).eq.'        ') then
        call getsmc(wetanl,len,lsoil,smcanl,me)
      endif
      call qcmxmn('smc1a   ',smcanl(1,1),slianl,snoanl,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('smc2a   ',smcanl(1,2),slianl,snoanl,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+8l] add smcanl(3:4)
      if(lsoil.gt.2) then
      call qcmxmn('smc3a   ',smcanl(1,3),slianl,snoanl,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('smc4a   ',smcanl(1,4),slianl,snoanl,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
      if(fnstca(1:8).eq.'        ') then
        call getstc(tsfanl,tg3anl,slianl,len,lsoil,stcanl,tsfimx)
      endif
      call qcmxmn('stc1a   ',stcanl(1,1),slianl,snoanl,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('stc2a   ',stcanl(1,2),slianl,snoanl,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+8l] add stcanl(3:4)
      if(lsoil.gt.2) then
      call qcmxmn('stc3a   ',stcanl(1,3),slianl,snoanl,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('stc4a   ',stcanl(1,4),slianl,snoanl,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
      call qcmxmn('vega    ',veganl,slianl,snoanl,icefl1,
     &            veglmx,veglmn,vegomx,vegomn,vegimx,vegimn,
     &            vegjmx,vegjmn,vegsmx,vegsmn,epsveg,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('veta    ',vetanl,slianl,snoanl,icefl1,
     &            vetlmx,vetlmn,vetomx,vetomn,vetimx,vetimn,
     &            vetjmx,vetjmn,vetsmx,vetsmn,epsvet,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('sota    ',sotanl,slianl,snoanl,icefl1,
     &            sotlmx,sotlmn,sotomx,sotomn,sotimx,sotimn,
     &            sotjmx,sotjmn,sotsmx,sotsmn,epssot,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+16l]----------------------------------------------------------------------
      call qcmxmn('vmna    ',vmnanl,slianl,snoanl,icefl1,
     &            vmnlmx,vmnlmn,vmnomx,vmnomn,vmnimx,vmnimn,
     &            vmnjmx,vmnjmn,vmnsmx,vmnsmn,epsvmn,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('vmxa    ',vmxanl,slianl,snoanl,icefl1,
     &            vmxlmx,vmxlmn,vmxomx,vmxomn,vmximx,vmximn,
     &            vmxjmx,vmxjmn,vmxsmx,vmxsmn,epsvmx,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('slpa    ',slpanl,slianl,snoanl,icefl1,
     &            slplmx,slplmn,slpomx,slpomn,slpimx,slpimn,
     &            slpjmx,slpjmn,slpsmx,slpsmn,epsslp,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('absa    ',absanl,slianl,snoanl,icefl1,
     &            abslmx,abslmn,absomx,absomn,absimx,absimn,
     &            absjmx,absjmn,abssmx,abssmn,epsabs,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu ----------------------------------------------------------------------------
!
!  monitoring prints
!
      if (monanl) then
       if (me .eq. 0) then
        print *,' '
        print *,'monitor of time and space interpolated analysis'
        print *,' '
!       call count(slianl,snoanl,len)
        print *,' '
        call monitr('tsfanl',tsfanl,slianl,snoanl,len)
        call monitr('albanl',albanl,slianl,snoanl,len)
        call monitr('aisanl',aisanl,slianl,snoanl,len)
        call monitr('snoanl',snoanl,slianl,snoanl,len)
        call monitr('scvanl',scvanl,slianl,snoanl,len)
        call monitr('smcanl1',smcanl(1,1),slianl,snoanl,len)
        call monitr('smcanl2',smcanl(1,2),slianl,snoanl,len)
        call monitr('stcanl1',stcanl(1,1),slianl,snoanl,len)
        call monitr('stcanl2',stcanl(1,2),slianl,snoanl,len)
!clu [+4l] add smcanl(3:4) and stcanl(3:4)
        if(lsoil.gt.2) then
        call monitr('smcanl3',smcanl(1,3),slianl,snoanl,len)
        call monitr('smcanl4',smcanl(1,4),slianl,snoanl,len)
        call monitr('stcanl3',stcanl(1,3),slianl,snoanl,len)
        call monitr('stcanl4',stcanl(1,4),slianl,snoanl,len)
        endif
        call monitr('tg3anl',tg3anl,slianl,snoanl,len)
        call monitr('zoranl',zoranl,slianl,snoanl,len)
!       if (gaus) then
          call monitr('cvaanl',cvanl ,slianl,snoanl,len)
          call monitr('cvbanl',cvbanl,slianl,snoanl,len)
          call monitr('cvtanl',cvtanl,slianl,snoanl,len)
!       endif
        call monitr('slianl',slianl,slianl,snoanl,len)
!       call monitr('plranl',plranl,slianl,snoanl,len)
        call monitr('orog  ',orog  ,slianl,snoanl,len)
        call monitr('veganl',veganl,slianl,snoanl,len)
        call monitr('vetanl',vetanl,slianl,snoanl,len)
        call monitr('sotanl',sotanl,slianl,snoanl,len)
!cwu [+2l] add sih, sic
        call monitr('sihanl',sihanl,slianl,snoanl,len)
        call monitr('sicanl',sicanl,slianl,snoanl,len)
!clu [+4l] add vmn, vmx, slp, abs
        call monitr('vmnanl',vmnanl,slianl,snoanl,len)
        call monitr('vmxanl',vmxanl,slianl,snoanl,len)
        call monitr('slpanl',slpanl,slianl,snoanl,len)
        call monitr('absanl',absanl,slianl,snoanl,len)
       endif

      endif
!
!  read in forecast fields if needed
!
      if (me .eq. 0) then
        write(6,*) '=============='
        write(6,*) '  fcst guess'
        write(6,*) '=============='
      endif
!
        percrit=critp2
!
      if(deads) then
!
!  fill in guess array with analysis if dead start.
!
        percrit=critp3
        if (me .eq. 0) write(6,*) 'this run is dead start run'
        call filfcs(tsffcs,wetfcs,snofcs,zorfcs,albfcs,
     &              tg3fcs,cvfcs ,cvbfcs,cvtfcs,
     &              cnpfcs,smcfcs,stcfcs,slifcs,aisfcs,
     &              vegfcs,vetfcs,sotfcs,alffcs,
!cwu [+1l] add ()fcs for sih, sic
     &              sihfcs,sicfcs,
!clu [+1l] add ()fcs for vmn, vmx, slp, abs
     &              vmnfcs,vmxfcs,slpfcs,absfcs,
     &              tsfanl,wetanl,snoanl,zoranl,albanl,
     &              tg3anl,cvanl ,cvbanl,cvtanl,
     &              cnpanl,smcanl,stcanl,slianl,aisanl,
     &              veganl,vetanl,sotanl,alfanl,
!cwu [+1l] add ()anl for sih, sic
     &              sihanl,sicanl,
!clu [+1l] add ()anl for vmn, vmx, slp, abs
     &              vmnanl,vmxanl,slpanl,absanl,     
     &              len,lsoil)
        if(sig1t(1).ne.0.) then
          call usesgt(sig1t,slianl,tg3anl,len,lsoil,tsffcs,stcfcs,
     &                tsfimx)
         do i=1,len
            icefl2(i) = sicfcs(i) .gt. 0.99999
          enddo
          kqcm=1
          call qcmxmn('tsff    ',tsffcs,slifcs,snofcs,icefl2,
     &                tsflmx,tsflmn,tsfomx,tsfomn,tsfimx,tsfimn,
     &                tsfjmx,tsfjmn,tsfsmx,tsfsmn,epstsf,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('stc1f   ',stcfcs(1,1),slifcs,snofcs,icefl1,
     &                stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &                stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('stc2f   ',stcfcs(1,2),slifcs,snofcs,icefl1,
     &                stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &                stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
        endif
      else
        percrit=critp2
!
!  make reverse angulation correction to tsf
!  make reverse orography correction to tg3
!
        if (use_ufo) then
          ztsfc = 1.0
          orogd = orog - orog_uf
          call tsfcor(tg3fcs,orogd,slmask,ztsfc,len,-rlapse)
          ztsfc = 0.
          call tsfcor(tsffcs,orogd,slmask,ztsfc,len,-rlapse)
        else
          ztsfc = 0.
          call tsfcor(tsffcs,orog,slmask,ztsfc,len,-rlapse)
        endif

!clu [+12l]  --------------------------------------------------------------
!
!  compute soil moisture liquid-to-total ratio over land
!
        do j=1, lsoil
        do i=1, len
         if(smcfcs(i,j) .ne. 0.)  then
            swratio(i,j) = slcfcs(i,j)/smcfcs(i,j)
           else
            swratio(i,j) = -999.
         endif
        enddo
        enddo
!clu -----------------------------------------------------------------------
!
        if(lqcbgs .and. irtacn .eq. 0) then
          call qcsli(slianl,slifcs,len,me)
          call albocn(albfcs,slmask,albomx,len)
         do i=1,len
            icefl2(i) = sicfcs(i) .gt. 0.99999
          enddo
          kqcm=1
          call qcmxmn('snof    ',snofcs,slifcs,snofcs,icefl1,
     &                snolmx,snolmn,snoomx,snoomn,snoimx,snoimn,
     &                snojmx,snojmn,snosmx,snosmn,epssno,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('tsff    ',tsffcs,slifcs,snofcs,icefl2,
     &                tsflmx,tsflmn,tsfomx,tsfomn,tsfimx,tsfimn,
     &                tsfjmx,tsfjmn,tsfsmx,tsfsmn,epstsf,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          do kk = 1, 4
          call qcmxmn('albf    ',albfcs(1,kk),slifcs,snofcs,icefl1,
     &                alblmx,alblmn,albomx,albomn,albimx,albimn,
     &                albjmx,albjmn,albsmx,albsmn,epsalb,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          enddo
        if(fnwetc(1:8).ne.'        ' .or. fnweta(1:8).ne.'        ' )
     &                                                          then
          call qcmxmn('wetf    ',wetfcs,slifcs,snofcs,icefl1,
     &                wetlmx,wetlmn,wetomx,wetomn,wetimx,wetimn,
     &                wetjmx,wetjmn,wetsmx,wetsmn,epswet,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
        endif
          call qcmxmn('zorf    ',zorfcs,slifcs,snofcs,icefl1,
     &                zorlmx,zorlmn,zoromx,zoromn,zorimx,zorimn,
     &                zorjmx,zorjmn,zorsmx,zorsmn,epszor,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
!       if(fnplrc(1:8).ne.'        ' .or. fnplra(1:8).ne.'        ' )
!         call qcmxmn('plnf    ',plrfcs,slifcs,snofcs,icefl1,
!    &                plrlmx,plrlmn,plromx,plromn,plrimx,plrimn,
!    &                plrjmx,plrjmn,plrsmx,plrsmn,epsplr,
!    &                rla,rlo,len,kqcm,percrit,lgchek,me)
!       endif
          call qcmxmn('tg3f    ',tg3fcs,slifcs,snofcs,icefl1,
     &                tg3lmx,tg3lmn,tg3omx,tg3omn,tg3imx,tg3imn,
     &                tg3jmx,tg3jmn,tg3smx,tg3smn,epstg3,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
!cwu [+8l] ---------------------------------------------------------------
          call qcmxmn('sihf    ',sihfcs,slifcs,snofcs,icefl1,
     &                sihlmx,sihlmn,sihomx,sihomn,sihimx,sihimn,
     &                sihjmx,sihjmn,sihsmx,sihsmn,epssih,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('sicf    ',sicfcs,slifcs,snofcs,icefl1,
     &                siclmx,siclmn,sicomx,sicomn,sicimx,sicimn,
     &                sicjmx,sicjmn,sicsmx,sicsmn,epssic,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('smc1f    ',smcfcs(1,1),slifcs,snofcs,icefl1,
     &                smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &                smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('smc2f   ',smcfcs(1,2),slifcs,snofcs,icefl1,
     &                smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &                smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+8l] add smcfcs(3:4)
          if(lsoil.gt.2) then
          call qcmxmn('smc3f    ',smcfcs(1,3),slifcs,snofcs,icefl1,
     &                smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &                smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('smc4f   ',smcfcs(1,4),slifcs,snofcs,icefl1,
     &                smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &                smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          endif
          call qcmxmn('stc1f   ',stcfcs(1,1),slifcs,snofcs,icefl1,
     &                stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &                stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('stc2f   ',stcfcs(1,2),slifcs,snofcs,icefl1,
     &                stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &                stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+8l] add stcfcs(3:4)
         if(lsoil.gt.2) then
          call qcmxmn('stc3f   ',stcfcs(1,3),slifcs,snofcs,icefl1,
     &                stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &                stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('stc4f   ',stcfcs(1,4),slifcs,snofcs,icefl1,
     &                stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &                stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
         endif
          call qcmxmn('vegf    ',vegfcs,slifcs,snofcs,icefl1,
     &                veglmx,veglmn,vegomx,vegomn,vegimx,vegimn,
     &                vegjmx,vegjmn,vegsmx,vegsmn,epsveg,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('vetf    ',vetfcs,slifcs,snofcs,icefl1,
     &                vetlmx,vetlmn,vetomx,vetomn,vetimx,vetimn,
     &                vetjmx,vetjmn,vetsmx,vetsmn,epsvet,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('sotf    ',sotfcs,slifcs,snofcs,icefl1,
     &                sotlmx,sotlmn,sotomx,sotomn,sotimx,sotimn,
     &                sotjmx,sotjmn,sotsmx,sotsmn,epssot,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)

!clu [+16l] ---------------------------------------------------------------
          call qcmxmn('vmnf    ',vmnfcs,slifcs,snofcs,icefl1,
     &                vmnlmx,vmnlmn,vmnomx,vmnomn,vmnimx,vmnimn,
     &                vmnjmx,vmnjmn,vmnsmx,vmnsmn,epsvmn,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('vmxf    ',vmxfcs,slifcs,snofcs,icefl1,
     &                vmxlmx,vmxlmn,vmxomx,vmxomn,vmximx,vmximn,
     &                vmxjmx,vmxjmn,vmxsmx,vmxsmn,epsvmx,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('slpf    ',slpfcs,slifcs,snofcs,icefl1,
     &                slplmx,slplmn,slpomx,slpomn,slpimx,slpimn,
     &                slpjmx,slpjmn,slpsmx,slpsmn,epsslp,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
          call qcmxmn('absf    ',absfcs,slifcs,snofcs,icefl1,
     &                abslmx,abslmn,absomx,absomn,absimx,absimn,
     &                absjmx,absjmn,abssmx,abssmn,epsabs,
     &                rla,rlo,len,kqcm,percrit,lgchek,me)
!clu -----------------------------------------------------------------------
        endif
      endif
!
      if (monfcs) then
       if (me .eq. 0) then
        print *,' '
        print *,'monitor of guess'
        print *,' '
!       call count(slifcs,snofcs,len)
        print *,' '
        call monitr('tsffcs',tsffcs,slifcs,snofcs,len)
        call monitr('albfcs',albfcs,slifcs,snofcs,len)
        call monitr('aisfcs',aisfcs,slifcs,snofcs,len)
        call monitr('snofcs',snofcs,slifcs,snofcs,len)
        call monitr('smcfcs1',smcfcs(1,1),slifcs,snofcs,len)
        call monitr('smcfcs2',smcfcs(1,2),slifcs,snofcs,len)
        call monitr('stcfcs1',stcfcs(1,1),slifcs,snofcs,len)
        call monitr('stcfcs2',stcfcs(1,2),slifcs,snofcs,len)
!clu [+4l] add smcfcs(3:4) and stcfcs(3:4)
        if(lsoil.gt.2) then
        call monitr('smcfcs3',smcfcs(1,3),slifcs,snofcs,len)
        call monitr('smcfcs4',smcfcs(1,4),slifcs,snofcs,len)
        call monitr('stcfcs3',stcfcs(1,3),slifcs,snofcs,len)
        call monitr('stcfcs4',stcfcs(1,4),slifcs,snofcs,len)
        endif
        call monitr('tg3fcs',tg3fcs,slifcs,snofcs,len)
        call monitr('zorfcs',zorfcs,slifcs,snofcs,len)
!       if (gaus) then
          call monitr('cvafcs',cvfcs ,slifcs,snofcs,len)
          call monitr('cvbfcs',cvbfcs,slifcs,snofcs,len)
          call monitr('cvtfcs',cvtfcs,slifcs,snofcs,len)
!       endif
        call monitr('slifcs',slifcs,slifcs,snofcs,len)
!       call monitr('plrfcs',plrfcs,slifcs,snofcs,len)
        call monitr('orog  ',orog  ,slifcs,snofcs,len)
        call monitr('vegfcs',vegfcs,slifcs,snofcs,len)
        call monitr('vetfcs',vetfcs,slifcs,snofcs,len)
        call monitr('sotfcs',sotfcs,slifcs,snofcs,len)
!cwu [+2l] add sih, sic
        call monitr('sihfcs',sihfcs,slifcs,snofcs,len)
        call monitr('sicfcs',sicfcs,slifcs,snofcs,len)
!clu [+4l] add vmn, vmx, slp, abs
        call monitr('vmnfcs',vmnfcs,slifcs,snofcs,len)
        call monitr('vmxfcs',vmxfcs,slifcs,snofcs,len)
        call monitr('slpfcs',slpfcs,slifcs,snofcs,len)
        call monitr('absfcs',absfcs,slifcs,snofcs,len)
       endif
      endif
!
!...   update annual cycle in the sst guess..
!
!     if(lprnt) print *,'tsfclm=',tsfclm(iprnt),' tsfcl2=',tsfcl2(iprnt)
!    *,' tsffcs=',tsffcs(iprnt),' slianl=',slianl(iprnt)

      if (fh > 0.0) then
        do i=1,len
          if(slianl(i) == 0.0) then
            tsffcs(i) = tsffcs(i) + (tsfclm(i) - tsfcl2(i))
          endif
        enddo
      endif
!
!  quality control analysis using forecast guess
!
      call qcbyfc(tsffcs,snofcs,qctsfs,qcsnos,qctsfi,len,lsoil,
     &            snoanl,aisanl,slianl,tsfanl,albanl,
     &            zoranl,smcanl,
     &            smcclm,tsfsmx,albomx,zoromx,me)
!
!  blend climatology and predicted fields
!
      if(me .eq. 0) then
        write(6,*) '=============='
        write(6,*) '   merging'
        write(6,*) '=============='
      endif
!     if(lprnt) print *,' tsffcs=',tsffcs(iprnt)
!
      percrit=critp3
!
!  merge analysis and forecast.  note tg3, ais are not merged
!
      call merge(len,lsoil,iy,im,id,ih,fh,deltsfc,
     &           sihfcs,sicfcs,
     &           vmnfcs,vmxfcs,slpfcs,absfcs, 
     &           tsffcs,wetfcs,snofcs,zorfcs,albfcs,aisfcs,
     &           cvfcs ,cvbfcs,cvtfcs,
     &           cnpfcs,smcfcs,stcfcs,slifcs,vegfcs,
     &           vetfcs,sotfcs,alffcs,
     &           sihanl,sicanl,                
     &           vmnanl,vmxanl,slpanl,absanl,       
     &           tsfanl,tsfan2,wetanl,snoanl,zoranl,albanl,aisanl,
     &           cvanl ,cvbanl,cvtanl,
     &           cnpanl,smcanl,stcanl,slianl,veganl,
     &           vetanl,sotanl,alfanl,
     &           ctsfl,calbl,caisl,csnol,csmcl,czorl,cstcl,cvegl,
     &           ctsfs,calbs,caiss,csnos,csmcs,czors,cstcs,cvegs,
     &           ccv,ccvb,ccvt,ccnp,cvetl,cvets,csotl,csots,
     &           calfl,calfs,
     &           csihl,csihs,csicl,csics,
     &           cvmnl,cvmns,cvmxl,cvmxs,cslpl,cslps,cabsl,cabss, 
     &           irttsf,irtwet,irtsno,irtzor,irtalb,irtais,
     &           irttg3,irtscv,irtacn,irtsmc,irtstc,irtveg,
     &           irtvmn,irtvmx,irtslp,irtabs,        
     &           irtvet,irtsot,irtalf,landice,me)

      call setzro(snoanl,epssno,len)

!     if(lprnt) print *,' tanlm=',tsfanl(iprnt),' tfcsm=',tsffcs(iprnt)
!     if(lprnt) print *,' sliam=',slianl(iprnt),' slifm=',slifcs(iprnt)

!
!  new ice/melted ice
!
      call newice(slianl,slifcs,tsfanl,tsffcs,len,lsoil,
!cwu [+1l] add sihnew, aislim, sihanl & sicanl
     &            sihnew,aislim,sihanl,sicanl,      
     &            albanl,snoanl,zoranl,smcanl,stcanl,
     &            albomx,snoomx,zoromx,smcomx,smcimx,
!cwu [-1l/+1l] change albimx to albimn - note albimx & albimn have been modified
!    &            tsfomn,tsfimx,albimx,zorimx,tgice,
     &            tsfomn,tsfimx,albimn,zorimx,tgice,
     &            rla,rlo,me)

!     if(lprnt) print *,'tsfanl=',tsfanl(iprnt),' tsffcs=',tsffcs(iprnt)
!     if(lprnt) print *,' slian=',slianl(iprnt),' slifn=',slifcs(iprnt)
!
!  set tsfc to tsnow over snow
!
      call snosfc(snoanl,tsfanl,tsfsmx,len,me)
!
      do i=1,len
        icefl2(i) = sicanl(i) .gt. 0.99999
      enddo
      kqcm=0
      call qcmxmn('snowm   ',snoanl,slianl,snoanl,icefl1,
     &            snolmx,snolmn,snoomx,snoomn,snoimx,snoimn,
     &            snojmx,snojmn,snosmx,snosmn,epssno,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('tsfm    ',tsfanl,slianl,snoanl,icefl2,
     &            tsflmx,tsflmn,tsfomx,tsfomn,tsfimx,tsfimn,
     &            tsfjmx,tsfjmn,tsfsmx,tsfsmn,epstsf,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      do kk = 1, 4
      call qcmxmn('albm    ',albanl(1,kk),slianl,snoanl,icefl1,
     &            alblmx,alblmn,albomx,albomn,albimx,albimn,
     &            albjmx,albjmn,albsmx,albsmn,epsalb,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      enddo
      if(fnwetc(1:8).ne.'        ' .or. fnweta(1:8).ne.'        ' )
     &                                                 then
      call qcmxmn('wetm    ',wetanl,slianl,snoanl,icefl1,
     &            wetlmx,wetlmn,wetomx,wetomn,wetimx,wetimn,
     &            wetjmx,wetjmn,wetsmx,wetsmn,epswet,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
      call qcmxmn('zorm    ',zoranl,slianl,snoanl,icefl1,
     &            zorlmx,zorlmn,zoromx,zoromn,zorimx,zorimn,
     &            zorjmx,zorjmn,zorsmx,zorsmn,epszor,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!     if(fnplrc(1:8).ne.'        ' .or. fnplra(1:8).ne.'        ' )
!    &                                                 then
!     call qcmxmn('plntm   ',plranl,slianl,snoanl,icefl1,
!    &            plrlmx,plrlmn,plromx,plromn,plrimx,plrimn,
!    &            plrjmx,plrjmn,plrsmx,plrsmn,epsplr,
!    &            rla,rlo,len,kqcm,percrit,lgchek,me)
!     endif
      call qcmxmn('stc1m   ',stcanl(1,1),slianl,snoanl,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('stc2m   ',stcanl(1,2),slianl,snoanl,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+8l] add stcanl(3:4)
       if(lsoil.gt.2) then
      call qcmxmn('stc3m   ',stcanl(1,3),slianl,snoanl,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('stc4m   ',stcanl(1,4),slianl,snoanl,icefl1,
     &            stclmx,stclmn,stcomx,stcomn,stcimx,stcimn,
     &            stcjmx,stcjmn,stcsmx,stcsmn,eptsfc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
       endif
      call qcmxmn('smc1m   ',smcanl(1,1),slianl,snoanl,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('smc2m   ',smcanl(1,2),slianl,snoanl,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+8l] add smcanl(3:4)
       if(lsoil.gt.2) then
      call qcmxmn('smc3m   ',smcanl(1,3),slianl,snoanl,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('smc4m   ',smcanl(1,4),slianl,snoanl,icefl1,
     &            smclmx,smclmn,smcomx,smcomn,smcimx,smcimn,
     &            smcjmx,smcjmn,smcsmx,smcsmn,epssmc,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      endif
      kqcm=1
      call qcmxmn('vegm    ',veganl,slianl,snoanl,icefl1,
     &            veglmx,veglmn,vegomx,vegomn,vegimx,vegimn,
     &            vegjmx,vegjmn,vegsmx,vegsmn,epsveg,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('vetm    ',vetanl,slianl,snoanl,icefl1,
     &            vetlmx,vetlmn,vetomx,vetomn,vetimx,vetimn,
     &            vetjmx,vetjmn,vetsmx,vetsmn,epsvet,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('sotm    ',sotanl,slianl,snoanl,icefl1,
     &            sotlmx,sotlmn,sotomx,sotomn,sotimx,sotimn,
     &            sotjmx,sotjmn,sotsmx,sotsmn,epssot,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!cwu [+8l] add sih, sic,
      call qcmxmn('sihm    ',sihanl,slianl,snoanl,icefl1,
     &            sihlmx,sihlmn,sihomx,sihomn,sihimx,sihimn,
     &            sihjmx,sihjmn,sihsmx,sihsmn,epssih,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('sicm    ',sicanl,slianl,snoanl,icefl1,
     &            siclmx,siclmn,sicomx,sicomn,sicimx,sicimn,
     &            sicjmx,sicjmn,sicsmx,sicsmn,epssic,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
!clu [+16l] add vmn, vmx, slp, abs
      call qcmxmn('vmnm    ',vmnanl,slianl,snoanl,icefl1,
     &            vmnlmx,vmnlmn,vmnomx,vmnomn,vmnimx,vmnimn,
     &            vmnjmx,vmnjmn,vmnsmx,vmnsmn,epsvmn,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('vmxm    ',vmxanl,slianl,snoanl,icefl1,
     &            vmxlmx,vmxlmn,vmxomx,vmxomn,vmximx,vmximn,
     &            vmxjmx,vmxjmn,vmxsmx,vmxsmn,epsvmx,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('slpm    ',slpanl,slianl,snoanl,icefl1,
     &            slplmx,slplmn,slpomx,slpomn,slpimx,slpimn,
     &            slpjmx,slpjmn,slpsmx,slpsmn,epsslp,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)
      call qcmxmn('absm    ',absanl,slianl,snoanl,icefl1,
     &            abslmx,abslmn,absomx,absomn,absimx,absimn,
     &            absjmx,absjmn,abssmx,abssmn,epsabs,
     &            rla,rlo,len,kqcm,percrit,lgchek,me)

!
      if(me .eq. 0) then
        write(6,*) '=============='
        write(6,*) 'final results'
        write(6,*) '=============='
      endif
!
!  foreward correction to tg3 and tsf at the last stage
!
!     if(lprnt) print *,' tsfbc=',tsfanl(iprnt)
      if (use_ufo) then
        ztsfc = 1.
        call tsfcor(tg3anl,orogd,slmask,ztsfc,len,rlapse)
        ztsfc = 0.
        call tsfcor(tsfanl,orogd,slmask,ztsfc,len,rlapse)
      else
        ztsfc = 0.
        call tsfcor(tsfanl,orog,slmask,ztsfc,len,rlapse)
      endif
!     if(lprnt) print *,' tsfaf=',tsfanl(iprnt)
!
!  check the final merged product
!
      if (monmer) then
       if(me .eq. 0) then
        print *,' '
        print *,'monitor of updated surface fields'
        print *,'   (includes angulation correction)'
        print *,' '
!       call count(slianl,snoanl,len)
        print *,' '
        call monitr('tsfanl',tsfanl,slianl,snoanl,len)
        call monitr('albanl',albanl,slianl,snoanl,len)
        call monitr('aisanl',aisanl,slianl,snoanl,len)
        call monitr('snoanl',snoanl,slianl,snoanl,len)
        call monitr('smcanl1',smcanl(1,1),slianl,snoanl,len)
        call monitr('smcanl2',smcanl(1,2),slianl,snoanl,len)
        call monitr('stcanl1',stcanl(1,1),slianl,snoanl,len)
        call monitr('stcanl2',stcanl(1,2),slianl,snoanl,len)
!clu [+4l] add smcanl(3:4) and stcanl(3:4)
        if(lsoil.gt.2) then
        call monitr('smcanl3',smcanl(1,3),slianl,snoanl,len)
        call monitr('smcanl4',smcanl(1,4),slianl,snoanl,len)
        call monitr('stcanl3',stcanl(1,3),slianl,snoanl,len)
        call monitr('stcanl4',stcanl(1,4),slianl,snoanl,len)
        call monitr('tg3anl',tg3anl,slianl,snoanl,len)
        call monitr('zoranl',zoranl,slianl,snoanl,len)
        endif
!       if (gaus) then
          call monitr('cvaanl',cvanl ,slianl,snoanl,len)
          call monitr('cvbanl',cvbanl,slianl,snoanl,len)
          call monitr('cvtanl',cvtanl,slianl,snoanl,len)
!       endif
        call monitr('slianl',slianl,slianl,snoanl,len)
!       call monitr('plranl',plranl,slianl,snoanl,len)
        call monitr('orog  ',orog  ,slianl,snoanl,len)
        call monitr('cnpanl',cnpanl,slianl,snoanl,len)
        call monitr('veganl',veganl,slianl,snoanl,len)
        call monitr('vetanl',vetanl,slianl,snoanl,len)
        call monitr('sotanl',sotanl,slianl,snoanl,len)
!cwu [+2l] add sih, sic,
        call monitr('sihanl',sihanl,slianl,snoanl,len)
        call monitr('sicanl',sicanl,slianl,snoanl,len)
!clu [+4l] add vmn, vmx, slp, abs
        call monitr('vmnanl',vmnanl,slianl,snoanl,len)
        call monitr('vmxanl',vmxanl,slianl,snoanl,len)
        call monitr('slpanl',slpanl,slianl,snoanl,len)
        call monitr('absanl',absanl,slianl,snoanl,len)
       endif
      endif
!
      if (mondif) then
        do i=1,len
          tsffcs(i) = tsfanl(i) - tsffcs(i)
          snofcs(i) = snoanl(i) - snofcs(i)
          tg3fcs(i) = tg3anl(i) - tg3fcs(i)
          zorfcs(i) = zoranl(i) - zorfcs(i)
!         plrfcs(i) = plranl(i) - plrfcs(i)
!         albfcs(i) = albanl(i) - albfcs(i)
          slifcs(i) = slianl(i) - slifcs(i)
          aisfcs(i) = aisanl(i) - aisfcs(i)
          cnpfcs(i) = cnpanl(i) - cnpfcs(i)
          vegfcs(i) = veganl(i) - vegfcs(i)
          vetfcs(i) = vetanl(i) - vetfcs(i)
          sotfcs(i) = sotanl(i) - sotfcs(i)
!clu [+2l] add sih, sic
          sihfcs(i) = sihanl(i) - sihfcs(i)
          sicfcs(i) = sicanl(i) - sicfcs(i)
!clu [+4l] add vmn, vmx, slp, abs
          vmnfcs(i) = vmnanl(i) - vmnfcs(i)
          vmxfcs(i) = vmxanl(i) - vmxfcs(i)
          slpfcs(i) = slpanl(i) - slpfcs(i)
          absfcs(i) = absanl(i) - absfcs(i)
        enddo
        do j = 1,lsoil
          do i = 1,len
            smcfcs(i,j) = smcanl(i,j) - smcfcs(i,j)
            stcfcs(i,j) = stcanl(i,j) - stcfcs(i,j)
          enddo
        enddo
        do j = 1,4
          do i = 1,len
            albfcs(i,j) = albanl(i,j) - albfcs(i,j)
          enddo
        enddo
!
!  monitoring prints
!
       if(me .eq. 0) then
        print *,' '
        print *,'monitor of difference'
        print *,'   (includes angulation correction)'
        print *,' '
        call monitr('tsfdif',tsffcs,slianl,snoanl,len)
        call monitr('albdif',albfcs,slianl,snoanl,len)
        call monitr('albdif1',albfcs,slianl,snoanl,len)
        call monitr('albdif2',albfcs(1,2),slianl,snoanl,len)
        call monitr('albdif3',albfcs(1,3),slianl,snoanl,len)
        call monitr('albdif4',albfcs(1,4),slianl,snoanl,len)
        call monitr('aisdif',aisfcs,slianl,snoanl,len)
        call monitr('snodif',snofcs,slianl,snoanl,len)
        call monitr('smcanl1',smcfcs(1,1),slianl,snoanl,len)
        call monitr('smcanl2',smcfcs(1,2),slianl,snoanl,len)
        call monitr('stcanl1',stcfcs(1,1),slianl,snoanl,len)
        call monitr('stcanl2',stcfcs(1,2),slianl,snoanl,len)
!clu [+4l] add smcfcs(3:4) and stc(3:4)
        if(lsoil.gt.2) then
        call monitr('smcanl3',smcfcs(1,3),slianl,snoanl,len)
        call monitr('smcanl4',smcfcs(1,4),slianl,snoanl,len)
        call monitr('stcanl3',stcfcs(1,3),slianl,snoanl,len)
        call monitr('stcanl4',stcfcs(1,4),slianl,snoanl,len)
        endif
        call monitr('tg3dif',tg3fcs,slianl,snoanl,len)
        call monitr('zordif',zorfcs,slianl,snoanl,len)
!       if (gaus) then
          call monitr('cvadif',cvfcs ,slianl,snoanl,len)
          call monitr('cvbdif',cvbfcs,slianl,snoanl,len)
          call monitr('cvtdif',cvtfcs,slianl,snoanl,len)
!       endif
        call monitr('slidif',slifcs,slianl,snoanl,len)
!       call monitr('plrdif',plrfcs,slianl,snoanl,len)
        call monitr('cnpdif',cnpfcs,slianl,snoanl,len)
        call monitr('vegdif',vegfcs,slianl,snoanl,len)
        call monitr('vetdif',vetfcs,slianl,snoanl,len)
        call monitr('sotdif',sotfcs,slianl,snoanl,len)
!cwu [+2l] add sih, sic
        call monitr('sihdif',sihfcs,slianl,snoanl,len)
        call monitr('sicdif',sicfcs,slianl,snoanl,len)
!clu [+4l] add vmn, vmx, slp, abs
        call monitr('vmndif',vmnfcs,slianl,snoanl,len)
        call monitr('vmxdif',vmxfcs,slianl,snoanl,len)
        call monitr('slpdif',slpfcs,slianl,snoanl,len)
        call monitr('absdif',absfcs,slianl,snoanl,len)
       endif
      endif
!
!
      do i=1,len
        tsffcs(i) = tsfanl(i)
        snofcs(i) = snoanl(i)
        tg3fcs(i) = tg3anl(i)
        zorfcs(i) = zoranl(i)
!       plrfcs(i) = plranl(i)
!       albfcs(i) = albanl(i)
        slifcs(i) = slianl(i)
        aisfcs(i) = aisanl(i)
        cvfcs(i)  = cvanl(i)
        cvbfcs(i) = cvbanl(i)
        cvtfcs(i) = cvtanl(i)
        cnpfcs(i) = cnpanl(i)
        vegfcs(i) = veganl(i)
        vetfcs(i) = vetanl(i)
        sotfcs(i) = sotanl(i)
!clu [+4l] add vmn, vmx, slp, abs
        vmnfcs(i) = vmnanl(i)
        vmxfcs(i) = vmxanl(i)
        slpfcs(i) = slpanl(i)
        absfcs(i) = absanl(i)
      enddo
      do j = 1,lsoil
        do i = 1,len
          smcfcs(i,j) = smcanl(i,j)
          if (slifcs(i) .gt. 0.0) then
             stcfcs(i,j) = stcanl(i,j)
          else
             stcfcs(i,j) = tsffcs(i)
          endif
        enddo
      enddo
      do j = 1,4
        do i = 1,len
          albfcs(i,j) = albanl(i,j)
        enddo
      enddo
      do j = 1,2
        do i = 1,len
          alffcs(i,j) = alfanl(i,j)
        enddo
      enddo

!cwu [+20l] update sihfcs, sicfcs. remove sea ice over non-ice points
      crit=aislim
      do i=1,len
        sihfcs(i) = sihanl(i)
        sitfcs(i) = tsffcs(i)
        if (slifcs(i).ge.2.) then
          if (sicfcs(i).gt.crit) then
            tsffcs(i) = (sicanl(i)*tsffcs(i)
     &                + (sicfcs(i)-sicanl(i))*tgice)/sicfcs(i)
            sitfcs(i) = (tsffcs(i)-tgice*(1.0-sicfcs(i))) / sicfcs(i)
          else
            tsffcs(i) = tsfanl(i)
!           tsffcs(i) = tgice
            sihfcs(i) = sihnew
          endif
        endif
        sicfcs(i) = sicanl(i)
      enddo
      do i=1,len
        if (slifcs(i).lt.1.5) then
          sihfcs(i) = 0.
          sicfcs(i) = 0.
          sitfcs(i) = tsffcs(i)
        else if ((slifcs(i).ge.1.5).and.(sicfcs(i).lt.crit)) then
          print *,'warning: check, slifcs and sicfcs',
     &            slifcs(i),sicfcs(i)
        endif
      enddo

!
! ensure the consistency between slc and smc
!
       do k=1, lsoil
        fixratio(k) = .false.
        if (fsmcl(k).lt.99999.) fixratio(k) = .true.
       enddo

       if(me .eq. 0) then
       print *,'dbgx --fixratio:',(fixratio(k),k=1,lsoil)
       endif

       do k=1, lsoil
        if(fixratio(k)) then
         do i = 1, len
           if(swratio(i,k) .eq. -999.) then
            slcfcs(i,k) = smcfcs(i,k)
           else
            slcfcs(i,k) = swratio(i,k) * smcfcs(i,k)
           endif
           if (slifcs(i) .ne. 1.0) slcfcs(i,k) = 1.0  ! flag value for non-land points.
         enddo
        endif
       enddo
! set liquid soil moisture to a flag value of 1.0
       if (landice) then
         do i = 1, len
           if (slifcs(i) .eq. 1.0 .and. 
     &         nint(vetfcs(i)) == veg_type_landice) then
             do k=1, lsoil
               slcfcs(i,k) = 1.0
             enddo
           endif
         enddo
       end if
!
! ensure the consistency between snwdph and sheleg
!
      if(fsnol .lt. 99999.) then  
       if(me .eq. 0) then
       print *,'dbgx -- scale snwdph from sheleg'
       endif
       do i = 1, len
        if(slifcs(i).eq.1.) swdfcs(i) = 10.* snofcs(i)
       enddo
      endif

! sea ice model only uses the liquid equivalent depth.
! so update the physical depth only for display purposes.
! use the same 3:1 ratio used by ice model.

      do i = 1, len
        if (slifcs(i).ne.1) swdfcs(i) = 3.*snofcs(i)
      enddo

      do i = 1, len
        if(slifcs(i).eq.1.) then
        if(snofcs(i).ne.0. .and. swdfcs(i).eq.0.) then
          print *,'dbgx --scale snwdph from sheleg',
     +        i, swdfcs(i), snofcs(i)
          swdfcs(i) = 10.* snofcs(i)
        endif
        endif
      enddo
! landice mods  - impose same minimum snow depth at
!                 landice as noah lsm.  also ensure
!                 lower thermal boundary condition
!                 and skin t is no warmer than freezing
!                 after adjustment to terrain.
       if (landice) then
         do i = 1, len
           if (slifcs(i) .eq. 1.0 .and. 
     &         nint(vetfcs(i)) == veg_type_landice) then
             snofcs(i) = max(snofcs(i),100.0)  ! in mm
             swdfcs(i) = max(swdfcs(i),1000.0) ! in mm
             tg3fcs(i) = min(tg3fcs(i),273.15)
             tsffcs(i) = min(tsffcs(i),273.15)
           endif
         enddo
       end if
!
!     if(lprnt) print *,' tsffcsf=',tsffcs(iprnt)
      return
      end subroutine sfccycle_sub
      subroutine count(slimsk,sno,ijmax)
      use machine , only : kind_io8,kind_io4
      implicit none
      real (kind=kind_io8) rl3,rl1,rl0,rl2,rl6,rl7,rl4,rl5
      integer l8,l7,l1,l2,ijmax,l0,l3,l5,l6,l4,ij
!
      real (kind=kind_io8) slimsk(1),sno(1)
!
!  count number of points for the four surface conditions
!
      l0 = 0
      l1 = 0
      l2 = 0
      l3 = 0
      l4 = 0
      do ij=1,ijmax
        if(slimsk(ij).eq.0.) l1 = l1 + 1
        if(slimsk(ij).eq.1. .and. sno(ij).le.0.) l0 = l0 + 1
        if(slimsk(ij).eq.2. .and. sno(ij).le.0.) l2 = l2 + 1
        if(slimsk(ij).eq.1. .and. sno(ij).gt.0.) l3 = l3 + 1
        if(slimsk(ij).eq.2. .and. sno(ij).gt.0.) l4 = l4 + 1
      enddo
      l5  = l0 + l3
      l6  = l2 + l4
      l7  = l1 + l6
      l8  = l1 + l5 + l6
      rl0 = float(l0) / float(l8)*100.
      rl3 = float(l3) / float(l8)*100.
      rl1 = float(l1) / float(l8)*100.
      rl2 = float(l2) / float(l8)*100.
      rl4 = float(l4) / float(l8)*100.
      rl5 = float(l5) / float(l8)*100.
      rl6 = float(l6) / float(l8)*100.
      rl7 = float(l7) / float(l8)*100.
      print *,'1) no. of not snow-covered land points   ',l0,' ',rl0,' '
      print *,'2) no. of snow covered land points       ',l3,' ',rl3,' '
      print *,'3) no. of open sea points                ',l1,' ',rl1,' '
      print *,'4) no. of not snow-covered seaice points ',l2,' ',rl2,' '
      print *,'5) no. of snow covered sea ice points    ',l4,' ',rl4,' '
      print *,' '
      print *,'6) no. of land points                    ',l5,' ',rl5,' '
      print *,'7) no. sea points (including sea ice)    ',l7,' ',rl7,' '
      print *,'   (no. of sea ice points)          (',l6,')',' ',rl6,' '
      print *,' '
      print *,'9) no. of total grid points               ',l8
!     print *,' '
!     print *,' '

!
!     if(lprnt) print *,' tsffcsf=',tsffcs(iprnt)
      return
      end
      subroutine monitr(lfld,fld,slimsk,sno,ijmax)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer ij,n,ijmax
!
      real (kind=kind_io8) fld(ijmax), slimsk(ijmax),sno(ijmax)
!
      real (kind=kind_io8) rmax(5),rmin(5)
      character(len=*) lfld
!
!  find max/min
!
      do n=1,5
        rmax(n) = -9.e20
        rmin(n) =  9.e20
      enddo
!
      do ij=1,ijmax
         if(slimsk(ij).eq.0.) then
            rmax(1) = max(rmax(1), fld(ij))
            rmin(1) = min(rmin(1), fld(ij))
         elseif(slimsk(ij).eq.1.) then
            if(sno(ij).le.0.) then
               rmax(2) = max(rmax(2), fld(ij))
               rmin(2) = min(rmin(2), fld(ij))
            else
               rmax(4) = max(rmax(4), fld(ij))
               rmin(4) = min(rmin(4), fld(ij))
            endif
         else
            if(sno(ij).le.0.) then
               rmax(3) = max(rmax(3), fld(ij))
               rmin(3) = min(rmin(3), fld(ij))
            else
               rmax(5) = max(rmax(5), fld(ij))
               rmin(5) = min(rmin(5), fld(ij))
            endif
         endif
      enddo
!
      print 100,lfld
      print 101,rmax(1),rmin(1)
      print 102,rmax(2),rmin(2), rmax(4), rmin(4)
      print 103,rmax(3),rmin(3), rmax(5), rmin(5)
!
!     print 102,rmax(2),rmin(2)
!     print 103,rmax(3),rmin(3)
!     print 104,rmax(4),rmin(4)
!     print 105,rmax(5),rmin(5)
  100 format('0  *** ',a8,' ***')
  101 format(' open sea  ......... max=',e12.4,' min=',e12.4)
  102 format(' land nosnow/snow .. max=',e12.4,' min=',e12.4
     &,                          ' max=',e12.4,' min=',e12.4)
  103 format(' seaice nosnow/snow  max=',e12.4,' min=',e12.4
     &,                          ' max=',e12.4,' min=',e12.4)
!
! 100 format('0',2x,'*** ',a8,' ***')
! 102 format(2x,' land without snow ..... max=',e12.4,' min=',e12.4)
! 103 format(2x,' seaice without snow ... max=',e12.4,' min=',e12.4)
! 104 format(2x,' land with snow ........ max=',e12.4,' min=',e12.4)
! 105 format(2x,' sea ice with snow ..... max=',e12.4,' min=',e12.4)
!
      return
      end
      subroutine dayoyr(iyr,imo,idy,ldy)
      implicit none
      integer ldy,i,idy,iyr,imo
!
!  this routine figures out the day of the year given imo and idy
!
      integer month(13)
      data month/0,31,28,31,30,31,30,31,31,30,31,30,31/
      if(mod(iyr,4).eq.0) month(3) = 29
      ldy = idy
      do i = 1, imo
        ldy = ldy + month(i)
      enddo
      return
      end
      subroutine hmskrd(lugb,imsk,jmsk,fnmskh,
     &                  kpds5,slmskh,gausm,blnmsk,bltmsk,me)
      use machine , only : kind_io8,kind_io4
      use sfccyc_module, only : mdata, xdata, ydata
      implicit none
      integer kpds5,me,i,imsk,jmsk,lugb
!
      character*500 fnmskh
!
      real (kind=kind_io8) slmskh(mdata)
      logical gausm
      real (kind=kind_io8) blnmsk,bltmsk
!
      imsk = xdata
      jmsk = ydata

      if (me .eq. 0) then
        write(6,*)' imsk=',imsk,' jmsk=',jmsk,' xdata=',xdata,' ydata='
     &,             ydata
      endif

      call fixrdg(lugb,imsk,jmsk,fnmskh,
     &            kpds5,slmskh,gausm,blnmsk,bltmsk,me)

!     print *,'in sfc_sub, aft fixrdg,slmskh=',maxval(slmskh),
!    &  minval(slmskh),'mdata=',mdata,'imsk*jmsk=',imsk*jmsk

      do i=1,imsk*jmsk
         slmskh(i) = nint(slmskh(i))
      enddo
!
      return
      end
      subroutine fixrdg(lugb,idim,jdim,fngrib,
     &                  kpds5,gdata,gaus,blno,blto,me)
      use machine      , only : kind_io8,kind_io4
      use sfccyc_module, only : mdata
      implicit none
      integer lgrib,n,lskip,jret,j,ndata,lugi,jdim,idim,lugb,
     &        iret, me,kpds5,kdata,i,w3kindreal,w3kindint
!
      character*(*) fngrib
!
      real (kind=kind_io8) gdata(idim*jdim)
      logical gaus
      real (kind=kind_io8) blno,blto
      real (kind=kind_io8), allocatable :: data8(:)
      real (kind=kind_io4), allocatable :: data4(:)
!
      logical, allocatable :: lbms(:)
!
      integer kpds(200),kgds(200)
      integer jpds(200),jgds(200), kpds0(200)
!
      allocate(data8(1:idim*jdim))
      allocate(lbms(1:mdata))
      kpds = 0
      kgds = 0
      jpds = 0
      jgds = 0
      kpds0 = 0
!
!     if(me .eq. 0) then
!     write(6,*) ' '
!     write(6,*) '************************************************'
!     endif
!
      close(lugb)
      call baopenr(lugb,fngrib,iret)
      if (iret .ne. 0) then
        write(6,*) ' error in opening file ',trim(fngrib)
        print *,'error in opening file ',trim(fngrib)
        call abort
      endif
      if (me .eq. 0) write(6,*) ' file ',trim(fngrib),
     &              ' opened. unit=',lugb
      lugi    = 0
      lskip   = -1
      n       = 0
      jpds    = -1
      jgds    = -1
      jpds(5) = kpds5
      kpds    = jpds
!
      call getgbh(lugb,lugi,lskip,jpds,jgds,lgrib,ndata,
     &            lskip,kpds,kgds,iret)
!
      if(me .eq. 0) then
        write(6,*) ' first grib record.'
        write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
        write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
        write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
      endif
!
      kpds0=jpds
      kpds0(4)=-1
      kpds0(18)=-1
      if(iret.ne.0) then
        write(6,*) ' error in getgbh. iret: ', iret
        if (iret == 99) write(6,*) ' field not found.'
        call abort
      endif
!
      jpds = kpds0
      lskip = -1
      kdata=idim*jdim
      call w3kind(w3kindreal,w3kindint)
      if (w3kindreal == 8) then
        call getgb(lugb,lugi,kdata,lskip,jpds,jgds,ndata,lskip,
     &             kpds,kgds,lbms,data8,jret)
      else if (w3kindreal == 4) then
        allocate(data4(1:idim*jdim))
        call getgb(lugb,lugi,kdata,lskip,jpds,jgds,ndata,lskip,
     &             kpds,kgds,lbms,data4,jret)
        data8 = real(data4, kind=kind_io8)
        deallocate(data4)
      else
        write(0,*)' Invalid w3kindreal --- aborting'
        call abort
      endif
!
      if(jret == 0) then
        if(ndata.eq.0) then
          write(6,*) ' error in getgb'
          write(6,*) ' kpds=',kpds
          write(6,*) ' kgds=',kgds
          call abort
        endif
        idim = kgds(2)
        jdim = kgds(3)
        gaus = kgds(1).eq.4
        blno = kgds(5)*1.d-3
        blto = kgds(4)*1.d-3
        gdata(1:idim*jdim) = data8(1:idim*jdim)
        if (me == 0) write(6,*) 'idim,jdim=',idim,jdim
     &,                ' gaus=',gaus,' blno=',blno,' blto=',blto
      else
        if (me ==. 0) write(6,*) 'idim,jdim=',idim,jdim
     &,                ' gaus=',gaus,' blno=',blno,' blto=',blto
        write(6,*) ' error in getgb : jret=',jret
        write(6,*) ' kpds(13)=',kpds(13),' kpds(15)=',kpds(15)
        call abort
      endif
!
      deallocate(data8)
      deallocate(lbms)
      return
      end
      subroutine getarea(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr
     &,                  me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer j,me,kgds11
      real (kind=kind_io8) f0lon,f0lat,elon,dlon,dlat,rslat,wlon,rnlat
!
!  get area of the grib record
!
      integer kgds(22)
      logical ijordr
!
      if (me .eq. 0) then
       write(6,*) ' kgds( 1-12)=',(kgds(j),j= 1,12)
       write(6,*) ' kgds(13-22)=',(kgds(j),j=13,22)
      endif
!
      if(kgds(1).eq.0) then                      !  lat/lon grid
!
        if (me .eq. 0) write(6,*) 'lat/lon grid'
        dlat   = float(kgds(10)) * 0.001
        dlon   = float(kgds( 9)) * 0.001
        f0lon  = float(kgds(5))  * 0.001
        f0lat  = float(kgds(4))  * 0.001
        kgds11 = kgds(11)
        if(kgds11.ge.128) then
          wlon = f0lon - dlon*(kgds(2)-1)
          elon = f0lon
          if(dlon*kgds(2).gt.359.99) then
            wlon =f0lon - dlon*kgds(2)
          endif
          dlon   = -dlon
          kgds11 = kgds11 - 128
        else
          wlon = f0lon
          elon = f0lon + dlon*(kgds(2)-1)
          if(dlon*kgds(2).gt.359.99) then
            elon = f0lon + dlon*kgds(2)
          endif
        endif
        if(kgds11.ge.64) then
          rnlat  = f0lat + dlat*(kgds(3)-1)
          rslat  = f0lat
          kgds11 = kgds11 - 64
        else
          rnlat = f0lat
          rslat = f0lat - dlat*(kgds(3)-1)
          dlat  = -dlat
        endif
        if(kgds11.ge.32) then
          ijordr = .false.
        else
          ijordr = .true.
        endif

        if(wlon.gt.180.) wlon = wlon - 360.
        if(elon.gt.180.) elon = elon - 360.
        wlon  = nint(wlon*1000.)  * 0.001
        elon  = nint(elon*1000.)  * 0.001
        rslat = nint(rslat*1000.) * 0.001
        rnlat = nint(rnlat*1000.) * 0.001
        return
!
      elseif(kgds(1).eq.1) then                  !  mercator projection
        write(6,*) 'mercator grid'
        write(6,*) 'cannot process'
        call abort
!
      elseif(kgds(1).eq.2) then                  !  gnomonic projection
        write(6,*) 'gnomonic grid'
        write(6,*) 'error!! gnomonic projection not coded'
        call abort
!
      elseif(kgds(1).eq.3) then                  !  lambert conformal
        write(6,*) 'lambert conformal'
        write(6,*) 'cannot process'
        call abort
      elseif(kgds(1).eq.4) then                  !  gaussian grid
!
        if (me .eq. 0) write(6,*) 'gaussian grid'
        dlat   = 99.
        dlon   = float(kgds( 9)) / 1000.0
        f0lon  = float(kgds(5))  / 1000.0
        f0lat  = 99.
        kgds11 = kgds(11)
        if(kgds11.ge.128) then
          wlon = f0lon
          elon = f0lon
          if(dlon*kgds(2).gt.359.99) then
            wlon = f0lon - dlon*kgds(2)
          endif
          dlon   = -dlon
          kgds11 = kgds11-128
        else
          wlon = f0lon
          elon = f0lon + dlon*(kgds(2)-1)
          if(dlon*kgds(2).gt.359.99) then
            elon = f0lon + dlon*kgds(2)
          endif
        endif
        if(kgds11.ge.64) then
          rnlat  = 99.
          rslat  = 99.
          kgds11 = kgds11 - 64
        else
          rnlat = 99.
          rslat = 99.
          dlat  = -99.
        endif
        if(kgds11.ge.32) then
          ijordr = .false.
        else
          ijordr = .true.
        endif
        return
!
      elseif(kgds(1).eq.5) then                  !  polar strereographic
        write(6,*) 'polar stereographic grid'
        write(6,*) 'cannot process'
        call abort
        return
!
      elseif(kgds(1).eq.13) then                 !  oblique lambert conformal
        write(6,*) 'oblique lambert conformal grid'
        write(6,*) 'cannot process'
        call abort
!
      elseif(kgds(1).eq.50) then                 !  spherical coefficient
        write(6,*) 'spherical coefficient'
        write(6,*) 'cannot process'
        call abort
        return
!
      elseif(kgds(1).eq.90) then                 !  space view perspective
!                                                  (orthographic grid)
        write(6,*) 'space view perspective grid'
        write(6,*) 'cannot process'
        call abort
        return
!
      else                                       !  unknown projection.  abort.
        write(6,*) 'error!! unknown map projection'
        write(6,*) 'kgds(1)=',kgds(1)
        print *,'error!! unknown map projection'
        print *,'kgds(1)=',kgds(1)
        call abort
      endif
!
      return
      end
      subroutine subst(data,imax,jmax,dlon,dlat,ijordr)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,j,ii,jj,jmax,imax,iret
      real (kind=kind_io8) dlat,dlon
!
      logical ijordr
!
      real (kind=kind_io8) data(imax,jmax)
      real (kind=kind_io8), allocatable ::  work(:,:)
!
      if(.not.ijordr.or.
     &  (ijordr.and.(dlat.gt.0..or.dlon.lt.0.))) then
        allocate (work(imax,jmax))

        if(.not.ijordr) then
          do j=1,jmax
            do i=1,imax
              work(i,j) = data(j,i)
            enddo
          enddo
        else
          do j=1,jmax
            do i=1,imax
              work(i,j) = data(i,j)
            enddo
          enddo
        endif
        if (dlat > 0.0) then
          if (dlon > 0.0) then
            do j=1,jmax
              jj = jmax - j + 1
              do i=1,imax
                data(i,jj) = work(i,j)
              enddo
            enddo
          else
            do i=1,imax
              data(imax-i+1,jj) = work(i,j)
            enddo
          endif
        else
          if (dlon > 0.0) then
            do j=1,jmax
              do i=1,imax
                data(i,j) = work(i,j)
              enddo
            enddo
          else
            do j=1,jmax
              do i=1,imax
                data(imax-i+1,j) = work(i,j)
              enddo
            enddo
          endif
        endif
        deallocate (work, stat=iret)
      endif
      return
      end
      subroutine la2ga(regin,imxin,jmxin,rinlon,rinlat,rlon,rlat,inttyp,
     &                 gauout,len,lmask,rslmsk,slmask
     &,                outlat,outlon,me)
      use machine , only : kind_io8,kind_io4
      implicit none
      real (kind=kind_io8) wei4,wei3,wei2,sum2,sum1,sum3,wei1,sum4,
     &                     wsum,tem,wsumiv,sums,sumn,wi2j2,x,y,wi1j1,
     &                     wi1j2,wi2j1,rlat,rlon,aphi,
     &                     rnume,alamd,denom
      integer jy,ifills,ix,len,inttyp,me,i,j,jmxin,imxin,jq,jx,j1,j2,
     &        ii,i1,i2,kmami,it
      integer nx,kxs,kxt
      integer, allocatable, save :: imxnx(:)
      integer, allocatable       :: ifill(:)
!
!  interpolation from lat/lon or gaussian grid to other lat/lon grid
!
      real (kind=kind_io8) outlon(len),outlat(len),gauout(len),
     &                     slmask(len)
      real (kind=kind_io8) regin (imxin,jmxin),rslmsk(imxin,jmxin)
!
      real (kind=kind_io8)    rinlat(jmxin),  rinlon(imxin)
      integer iindx1(len),    iindx2(len)
      integer jindx1(len),    jindx2(len)
      real (kind=kind_io8)    ddx(len),       ddy(len),   wrk(len)
!
      logical lmask
!
      logical first
      integer   num_threads
      data first /.true./
      save num_threads, first
!
      integer len_thread_m, len_thread, i1_t, i2_t
      integer num_parthds
!
      if (first) then
         num_threads = num_parthds()
         first = .false.
         if (.not. allocated(imxnx)) allocate (imxnx(num_threads))
      endif
!
      if (me == 0) print *,' num_threads =',num_threads,' me=',me
!
!     if(me .eq. 0) then
!     print *,'rlon=',rlon,' me=',me
!     print *,'rlat=',rlat,' me=',me,' imxin=',imxin,' jmxin=',jmxin
!     endif
!
!     do j=1,jmxin
!       if(rlat.gt.0.) then
!         rinlat(j) = rlat - float(j-1)*dlain
!       else
!         rinlat(j) = rlat + float(j-1)*dlain
!       endif
!     enddo
!
!     if (me .eq. 0) then
!       print *,'rinlat='
!       print *,(rinlat(j),j=1,jmxin)
!       print *,'rinlon='
!       print *,(rinlon(i),i=1,imxin)
!
!       print *,'outlat='
!       print *,(outlat(j),j=1,len)
!       print *,(outlon(j),j=1,len)
!     endif
!
!     do i=1,imxin
!       rinlon(i) = rlon + float(i-1)*dloin
!     enddo
!
!     print *,'rinlon='
!     print *,(rinlon(i),i=1,imxin)
!
      len_thread_m  = (len+num_threads-1) / num_threads

      if (inttyp /=1) allocate (ifill(num_threads))
!
!$omp parallel do default(none) 
!$omp+private(i1_t,i2_t,len_thread,it,i,ii,i1,i2)
!$omp+private(j,j1,j2,jq,ix,jy,nx,kxs,kxt,kmami)
!$omp+private(alamd,denom,rnume,aphi,x,y,wsum,wsumiv,sum1,sum2)
!$omp+private(sum3,sum4,wi1j1,wi2j1,wi1j2,wi2j2,wei1,wei2,wei3,wei4)
!$omp+private(sumn,sums)
!$omp+shared(imxin,jmxin,ifill)
!$omp+shared(outlon,outlat,wrk,iindx1,rinlon,jindx1,rinlat,ddx,ddy)
!$omp+shared(rlon,rlat,regin,gauout,imxnx)
!$omp+private(tem)
!$omp+shared(num_threads,len_thread_m,len,lmask,iindx2,jindx2,rslmsk)
!$omp+shared(inttyp,me,slmask)
!
      do it=1,num_threads   ! start of threaded loop ...................
        i1_t       = (it-1)*len_thread_m+1
        i2_t       = min(i1_t+len_thread_m-1,len)
        len_thread = i2_t-i1_t+1
!
!       find i-index for interpolation
!
        do i=i1_t, i2_t
          alamd = outlon(i)
          if (alamd .lt. rlon)   alamd = alamd + 360.0
          if (alamd .gt. 360.0+rlon) alamd = alamd - 360.0
          wrk(i)    = alamd
          iindx1(i) = imxin
        enddo
        do i=i1_t,i2_t
          do ii=1,imxin
            if(wrk(i) .ge. rinlon(ii)) iindx1(i) = ii
          enddo
        enddo
        do i=i1_t,i2_t
          i1 = iindx1(i)
          if (i1 .lt. 1) i1 = imxin
          i2 = i1 + 1
          if (i2 .gt. imxin) i2 = 1
          iindx1(i) = i1
          iindx2(i) = i2
          denom     = rinlon(i2) - rinlon(i1)
          if(denom.lt.0.) denom = denom + 360.
          rnume = wrk(i) - rinlon(i1)
          if(rnume.lt.0.) rnume = rnume + 360.
          ddx(i) = rnume / denom
        enddo
!
!  find j-index for interplation
!
        if(rlat.gt.0.) then
          do j=i1_t,i2_t
            jindx1(j)=0
          enddo
          do jx=1,jmxin
            do j=i1_t,i2_t
              if(outlat(j).le.rinlat(jx)) jindx1(j) = jx
            enddo
          enddo
          do j=i1_t,i2_t
            jq = jindx1(j)
            aphi=outlat(j)
            if(jq.ge.1 .and. jq .lt. jmxin) then
              j2=jq+1
              j1=jq
             ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
            elseif (jq .eq. 0) then
              j2=1
              j1=1
              if(abs(90.-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            else
              j2=jmxin
              j1=jmxin
              if(abs(-90.-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            endif
            jindx1(j)=j1
            jindx2(j)=j2
          enddo
        else
          do j=i1_t,i2_t
            jindx1(j) = jmxin+1
          enddo
          do jx=jmxin,1,-1
            do j=i1_t,i2_t
              if(outlat(j).le.rinlat(jx)) jindx1(j) = jx
            enddo
          enddo
          do j=i1_t,i2_t
            jq = jindx1(j)
            aphi=outlat(j)
            if(jq.gt.1 .and. jq .le. jmxin) then
              j2=jq
              j1=jq-1
              ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
            elseif (jq .eq. 1) then
              j2=1
              j1=1
              if(abs(-90.-rinlat(j1)).gt.0.001) then
                 ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
              else
                 ddy(j)=0.0
              endif
            else
              j2=jmxin
              j1=jmxin
              if(abs(90.-rinlat(j1)).gt.0.001) then
                 ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
              else
                 ddy(j)=0.0
              endif
            endif
            jindx1(j)=j1
            jindx2(j)=j2
          enddo
        endif
!
!     if (me .eq. 0 .and. inttyp .eq. 1) then
!       print *,'la2ga'
!       print *,'iindx1'
!       print *,(iindx1(n),n=1,len)
!       print *,'iindx2'
!       print *,(iindx2(n),n=1,len)
!       print *,'jindx1'
!       print *,(jindx1(n),n=1,len)
!       print *,'jindx2'
!       print *,(jindx2(n),n=1,len)
!       print *,'ddy'
!       print *,(ddy(n),n=1,len)
!       print *,'ddx'
!       print *,(ddx(n),n=1,len)
!     endif
!
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        if (lmask) then
          wei1 = 0.
          wei2 = 0.
          wei3 = 0.
          wei4 = 0.
          do i=1,imxin
            sum1 = sum1 + regin(i,1) * rslmsk(i,1)
            sum2 = sum2 + regin(i,jmxin) * rslmsk(i,jmxin)
            wei1 = wei1 + rslmsk(i,1)
            wei2 = wei2 + rslmsk(i,jmxin)
!
            sum3 = sum3 + regin(i,1) * (1.0-rslmsk(i,1))
            sum4 = sum4 + regin(i,jmxin) * (1.0-rslmsk(i,jmxin))
            wei3 = wei3 + (1.0-rslmsk(i,1))
            wei4 = wei4 + (1.0-rslmsk(i,jmxin))
          enddo
!
          if(wei1.gt.0.) then
            sum1 = sum1 / wei1
          else
            sum1 = 0.
          endif
          if(wei2.gt.0.) then
            sum2 = sum2 / wei2
          else
            sum2 = 0.
          endif
          if(wei3.gt.0.) then
            sum3 = sum3 / wei3
          else
            sum3 = 0.
          endif
          if(wei4.gt.0.) then
            sum4 = sum4 / wei4
          else
            sum4 = 0.
          endif
        else
          do i=1,imxin
            sum1 = sum1 + regin(i,1)
            sum2 = sum2 + regin(i,jmxin)
          enddo
          sum1 = sum1 / imxin
          sum2 = sum2 / imxin
          sum3 = sum1
          sum4 = sum2
        endif
!
!     print *,' sum1=',sum1,' sum2=',sum2
!    *,' sum3=',sum3,' sum4=',sum4
!     print *,' rslmsk=',(rslmsk(i,1),i=1,imxin)
!     print *,' slmask=',(slmask(i),i=1,imxout)
!    *,' j1=',jindx1(1),' j2=',jindx2(1)
!
!
!  inttyp=1  take the closest point value
!
        if(inttyp.eq.1) then

          do i=i1_t,i2_t
            jy = jindx1(i)
            if(ddy(i) .ge. 0.5) jy = jindx2(i)
            ix = iindx1(i)
            if(ddx(i) .ge. 0.5) ix = iindx2(i)
!
!cggg start
!
            if (.not. lmask) then

              gauout(i) = regin(ix,jy)

            else

              if(slmask(i).eq.rslmsk(ix,jy)) then

                gauout(i) = regin(ix,jy)

              else

                i1 = ix
                j1 = jy

! spiral around until matching mask is found.
                do nx=1,jmxin*imxin/2
                  kxs=sqrt(4*nx-2.5)
                  kxt=nx-int(kxs**2/4+1)
                  select case(mod(kxs,4))
                  case(1)
                    ix=i1-kxs/4+kxt
                    jx=j1-kxs/4
                  case(2)
                    ix=i1+1+kxs/4
                    jx=j1-kxs/4+kxt
                  case(3)
                    ix=i1+1+kxs/4-kxt
                    jx=j1+1+kxs/4
                  case default
                    ix=i1-kxs/4
                    jx=j1+kxs/4-kxt
                  end select
                  if(jx.lt.1) then
                    ix=ix+imxin/2
                    jx=2-jx
                  elseif(jx.gt.jmxin) then
                    ix=ix+imxin/2
                    jx=2*jmxin-jx
                  endif
                  ix=modulo(ix-1,imxin)+1
                  if(slmask(i).eq.rslmsk(ix,jx)) then
                    gauout(i) = regin(ix,jx)
                    go to 81
                  endif
                enddo

!cggg here, set the gauout value to be 0, and let's sarah's land
!cggg routine assign a default.

              if (num_threads == 1) then
                print*,'no matching mask found ',i,i1,j1,ix,jx
                print*,'set to default value.'
              endif
              gauout(i) = 0.0


   81  continue

              end if

            end if

!cggg end

          enddo
          kmami=1
          if (me == 0 .and. num_threads == 1)
     &                  call maxmin(gauout(i1_t),len_thread,kmami)
        else  ! nearest neighbor interpolation

!
!  quasi-bilinear interpolation
!
          ifill(it) = 0
          imxnx(it) = 0
          do i=i1_t,i2_t
            y  = ddy(i)
            j1 = jindx1(i)
            j2 = jindx2(i)
            x  = ddx(i)
            i1 = iindx1(i)
            i2 = iindx2(i)
!
            wi1j1 = (1.-x) * (1.-y)
            wi2j1 =     x  *( 1.-y)
            wi1j2 = (1.-x) *      y
            wi2j2 =     x  *      y
!
            tem = 4.*slmask(i) - rslmsk(i1,j1) - rslmsk(i2,j1)
     &                         - rslmsk(i1,j2) - rslmsk(i2,j2)
            if(lmask .and. abs(tem) .gt. 0.01) then
              if(slmask(i).eq.1.) then
                  wi1j1 = wi1j1 * rslmsk(i1,j1)
                  wi2j1 = wi2j1 * rslmsk(i2,j1)
                  wi1j2 = wi1j2 * rslmsk(i1,j2)
                  wi2j2 = wi2j2 * rslmsk(i2,j2)
              else
                  wi1j1 = wi1j1 * (1.0-rslmsk(i1,j1))
                  wi2j1 = wi2j1 * (1.0-rslmsk(i2,j1))
                  wi1j2 = wi1j2 * (1.0-rslmsk(i1,j2))
                  wi2j2 = wi2j2 * (1.0-rslmsk(i2,j2))
              endif
            endif
!
            wsum   = wi1j1 + wi2j1 + wi1j2 + wi2j2
            wrk(i) = wsum
            if(wsum.ne.0.) then
              wsumiv = 1./wsum
!
              if(j1.ne.j2) then
                gauout(i) = (wi1j1*regin(i1,j1) + wi2j1*regin(i2,j1) +
     &                       wi1j2*regin(i1,j2) + wi2j2*regin(i2,j2))
     &                    *wsumiv
              else
!
                if (rlat .gt. 0.0) then
                  if (slmask(i) .eq. 1.0) then
                    sumn = sum1
                    sums = sum2
                  else
                    sumn = sum3
                    sums = sum4
                  endif
                  if( j1 .eq. 1) then
                    gauout(i) = (wi1j1*sumn        +wi2j1*sumn        +
     &                           wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))
     &                        * wsumiv
                  elseif (j1 .eq. jmxin) then
                    gauout(i) = (wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+
     &                           wi1j2*sums        +wi2j2*sums        )
     &                        * wsumiv
                  endif
!       print *,' slmask=',slmask(i),' sums=',sums,' sumn=',sumn
!    &  ,' regin=',regin(i1,j2),regin(i2,j2),' j1=',j1,' j2=',j2
!    &  ,' wij=',wi1j1, wi2j1, wi1j2, wi2j2,wsumiv
                else
                  if (slmask(i) .eq. 1.0) then
                    sums = sum1
                    sumn = sum2
                  else
                    sums = sum3
                    sumn = sum4
                  endif
                  if( j1 .eq. 1) then
                    gauout(i) = (wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+
     &                           wi1j2*sums        +wi2j2*sums        )
     &                        * wsumiv
                  elseif (j1 .eq. jmxin) then
                    gauout(i) = (wi1j1*sumn        +wi2j1*sumn        +
     &                           wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))
     &                        * wsumiv
                  endif
                endif
              endif            ! if j1 .ne. j2
            endif
          enddo
          do i=i1_t,i2_t
            j1 = jindx1(i)
            j2 = jindx2(i)
            i1 = iindx1(i)
            i2 = iindx2(i)
            if(wrk(i) .eq. 0.0) then
              if(.not.lmask) then
                if (num_threads == 1)
     &            write(6,*) ' la2ga called with lmask=.true. but bad',
     &                     ' rslmsk or slmask given'
                call abort
              endif
              ifill(it) = ifill(it) + 1
              if(ifill(it) <= 2 ) then
                if (me == 0 .and. num_threads == 1) then
                  write(6,*) 'i1,i2,j1,j2=',i1,i2,j1,j2
                  write(6,*) 'rslmsk=',rslmsk(i1,j1),rslmsk(i1,j2),
     &                                 rslmsk(i2,j1),rslmsk(i2,j2)
!                 write(6,*) 'i,j=',i,j,' slmask(i)=',slmask(i)
                  write(6,*) 'i=',i,' slmask(i)=',slmask(i)
     &,           ' outlon=',outlon(i),' outlat=',outlat(i)
                endif
              endif
! spiral around until matching mask is found.
              do nx=1,jmxin*imxin/2
                kxs=sqrt(4*nx-2.5)
                kxt=nx-int(kxs**2/4+1)
                select case(mod(kxs,4))
                case(1)
                  ix=i1-kxs/4+kxt
                  jx=j1-kxs/4
                case(2)
                  ix=i1+1+kxs/4
                  jx=j1-kxs/4+kxt
                case(3)
                  ix=i1+1+kxs/4-kxt
                  jx=j1+1+kxs/4
                case default
                  ix=i1-kxs/4
                  jx=j1+kxs/4-kxt
                end select
                if(jx.lt.1) then
                  ix=ix+imxin/2
                  jx=2-jx
                elseif(jx.gt.jmxin) then
                  ix=ix+imxin/2
                  jx=2*jmxin-jx
                endif
                ix=modulo(ix-1,imxin)+1
                if(slmask(i).eq.rslmsk(ix,jx)) then
                  gauout(i) = regin(ix,jx)
                  imxnx(it) = max(imxnx(it),nx)
                  go to 71
                endif
              enddo
!
              if (num_threads == 1) then
                write(6,*) ' error!!! no filling value found in la2ga'
!               write(6,*) ' i ix jx slmask(i) rslmsk ',
!    &                       i,ix,jx,slmask(i),rslmsk(ix,jx)
              endif
              call abort
!
   71         continue
            endif
!
          enddo
        endif
      enddo            ! end of threaded loop ...................
!$omp end parallel do
!
      if(inttyp /= 1)then
        ifills = 0
        do it=1,num_threads
          ifills = ifills + ifill(it)
        enddo

        if(ifills.gt.1) then
          if (me .eq. 0) then
          write(6,*) ' unable to interpolate.  filled with nearest',
     &               ' point value at ',ifills,' points'
!    &               ' point value at ',ifills,' points  imxnx=',imxnx(:)
          endif
        endif
        deallocate (ifill)
      endif
!
      kmami=1
      if (me .eq. 0) call maxmin(gauout,len,kmami)
!
      return
      end subroutine la2ga
      subroutine maxmin(f,imax,kmax)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,iimin,iimax,kmax,imax,k
      real (kind=kind_io8) fmin,fmax
!
      real (kind=kind_io8) f(imax,kmax)
!
      do k=1,kmax
!
        fmax = f(1,k)
        fmin = f(1,k)
!
        do i=1,imax
          if(fmax.le.f(i,k)) then
            fmax  = f(i,k)
            iimax = i
          endif
          if(fmin.ge.f(i,k)) then
            fmin  = f(i,k)
            iimin = i
          endif
        enddo
!
      write(6,100) k,fmax,iimax,fmin,iimin
  100 format(2x,'level=',i2,' max=',e11.4,' at i=',i7,
     &                      ' min=',e11.4,' at i=',i7)
!
      enddo
!
      return
      end
      subroutine filanl(tsfanl,tsfan2,wetanl,snoanl,zoranl,albanl,
     &                  aisanl,
     &                  tg3anl,cvanl ,cvbanl,cvtanl,
     &                  cnpanl,smcanl,stcanl,slianl,scvanl,veganl,
     &                  vetanl,sotanl,alfanl,
!cwu [+1l] add ()anl for sih, sic
     &                  sihanl,sicanl,
!clu [+1l] add ()anl for vmn, vmx, slp, abs
     &                  vmnanl,vmxanl,slpanl,absanl,
     &                  tsfclm,tsfcl2,wetclm,snoclm,zorclm,albclm,
     &                  aisclm,
     &                  tg3clm,cvclm ,cvbclm,cvtclm,
     &                  cnpclm,smcclm,stcclm,sliclm,scvclm,vegclm,
     &                  vetclm,sotclm,alfclm,
!cwu [+1l] add ()clm for sih, sic
     &                  sihclm,sicclm,
!clu [+1l] add ()clm for vmn, vmx, slp, abs
     &                  vmnclm,vmxclm,slpclm,absclm,
     &                  len,lsoil)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,j,len,lsoil
!
      real (kind=kind_io8) tsfanl(len),tsfan2(len),wetanl(len),
     &     snoanl(len),
     &     zoranl(len),albanl(len,4),aisanl(len),
     &     tg3anl(len),
     &     cvanl (len),cvbanl(len),cvtanl(len),
     &     cnpanl(len),
     &     smcanl(len,lsoil),stcanl(len,lsoil),
     &     slianl(len),scvanl(len),veganl(len),
     &     vetanl(len),sotanl(len),alfanl(len,2)
!cwu [+1l] add ()anl for sih, sic
     &,    sihanl(len),sicanl(len)
!clu [+1l] add ()anl for vmn, vmx, slp, abs
     &,    vmnanl(len),vmxanl(len),slpanl(len),absanl(len)
      real (kind=kind_io8) tsfclm(len),tsfcl2(len),wetclm(len),
     &     snoclm(len),
     &     zorclm(len),albclm(len,4),aisclm(len),
     &     tg3clm(len),
     &     cvclm (len),cvbclm(len),cvtclm(len),
     &     cnpclm(len),
     &     smcclm(len,lsoil),stcclm(len,lsoil),
     &     sliclm(len),scvclm(len),vegclm(len),
     &     vetclm(len),sotclm(len),alfclm(len,2)
!cwu [+1l] add ()clm for sih, sic
     &,    sihclm(len),sicclm(len)
!clu [+1l] add ()clm for vmn, vmx, slp, abs
     &,    vmnclm(len),vmxclm(len),slpclm(len),absclm(len)
!
      do i=1,len
        tsfanl(i)   = tsfclm(i)      !  tsf at t
        tsfan2(i)   = tsfcl2(i)      !  tsf at t-deltsfc
        wetanl(i)   = wetclm(i)      !  soil wetness
        snoanl(i)   = snoclm(i)      !  snow
        scvanl(i)   = scvclm(i)      !  snow cover
        aisanl(i)   = aisclm(i)      !  seaice
        slianl(i)   = sliclm(i)      !  land/sea/snow mask
        zoranl(i)   = zorclm(i)      !  surface roughness
!       plranl(i)   = plrclm(i)      !  maximum stomatal resistance
        tg3anl(i)   = tg3clm(i)      !  deep soil temperature
        cnpanl(i)   = cnpclm(i)      !  canopy water content
        veganl(i)   = vegclm(i)      !  vegetation cover
        vetanl(i)   = vetclm(i)      !  vegetation type
        sotanl(i)   = sotclm(i)      !  soil type
        cvanl(i)    = cvclm(i)       !  cv
        cvbanl(i)   = cvbclm(i)      !  cvb
        cvtanl(i)   = cvtclm(i)      !  cvt
!cwu [+4l] add sih, sic
        sihanl(i)   = sihclm(i)      !  sea ice thickness
        sicanl(i)   = sicclm(i)      !  sea ice concentration
!clu [+4l] add vmn, vmx, slp, abs
        vmnanl(i)   = vmnclm(i)      !  min vegetation cover
        vmxanl(i)   = vmxclm(i)      !  max vegetation cover 
        slpanl(i)   = slpclm(i)      !  slope type
        absanl(i)   = absclm(i)      !  max snow albedo
      enddo
!
      do j=1,lsoil
        do i=1,len
          smcanl(i,j) = smcclm(i,j)  !   layer soil wetness
          stcanl(i,j) = stcclm(i,j)  !   soil temperature
        enddo
      enddo
      do j=1,4
        do i=1,len
          albanl(i,j) = albclm(i,j)  !  albedo
        enddo
      enddo
      do j=1,2
        do i=1,len
          alfanl(i,j) = alfclm(i,j)  !  vegetation fraction for albedo
        enddo
      enddo
!
      return
      end
      subroutine analy(lugb,iy,im,id,ih,fh,len,lsoil,
     &                 slmask,fntsfa,fnweta,fnsnoa,fnzora,fnalba,fnaisa,
     &                 fntg3a,fnscva,fnsmca,fnstca,fnacna,fnvega,
     &                 fnveta,fnsota,
!clu [+1l] add fn()a for vmn, vmx, slp, abs
     &                 fnvmna,fnvmxa,fnslpa,fnabsa,
     &                 tsfanl,wetanl,snoanl,zoranl,albanl,aisanl,
     &                 tg3anl,cvanl ,cvbanl,cvtanl,
     &                 smcanl,stcanl,slianl,scvanl,acnanl,veganl,
     &                 vetanl,sotanl,alfanl,tsfan0,
!clu [+1l] add ()anl for vmn, vmx, slp, abs
     &                 vmnanl,vmxanl,slpanl,absanl,
!cggg snow mods start    &        kpdtsf,kpdwet,kpdsno,kpdzor,kpdalb,kpdais,
     &                 kpdtsf,kpdwet,kpdsno,kpdsnd,kpdzor,kpdalb,kpdais,
!cggg snow mods end
     &                 kpdtg3,kpdscv,kpdacn,kpdsmc,kpdstc,kpdveg,
     &                 kprvet,kpdsot,kpdalf,
!clu [+1l] add kpd() for vmn, vmx, slp, abs
     &                 kpdvmn,kpdvmx,kpdslp,kpdabs,
     &                 irttsf,irtwet,irtsno,irtzor,irtalb,irtais,
     &                 irttg3,irtscv,irtacn,irtsmc,irtstc,irtveg,
     &                 irtvet,irtsot,irtalf
!clu [+1l] add irt() for vmn, vmx, slp, abs
     &,                irtvmn,irtvmx,irtslp,irtabs
     &,                imsk, jmsk, slmskh, outlat, outlon
     &,                gaus, blno, blto, me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer irtsmc,irtacn,irtstc,irtvet,irtveg,irtscv,irtzor,irtsno,
     &        irtalb,irttg3,irtais,iret,me,kk,kpdvet,i,irtalf,irtsot,
!cggg snow mods start     & imsk,jmsk,irtwet,lsoil,len, kpdtsf,kpdsno,kpdwet,iy,
     &        imsk,jmsk,irtwet,lsoil,len,kpdtsf,kpdsno,kpdsnd,kpdwet,iy,
!cggg snow mods end
     &        lugb,im,ih,id,kpdveg,kpdstc,kprvet,irttsf,kpdsot,kpdsmc,
     &        kpdais,kpdzor,kpdtg3,kpdacn,kpdscv,j
!clu [+1l] add kpd() and irt() for vmn, vmx, slp, abs
     &,       kpdvmn,kpdvmx,kpdslp,kpdabs,irtvmn,irtvmx,irtslp,irtabs
      real (kind=kind_io8) blto,blno,fh
!
      real (kind=kind_io8)    slmask(len)
      real (kind=kind_io8)    slmskh(imsk,jmsk)
      real (kind=kind_io8)    outlat(len), outlon(len)
      integer kpdalb(4),   kpdalf(2)
!cggg snow mods start
      integer kpds(1000),kgds(1000),jpds(1000),jgds(1000)
      integer lugi, lskip, lgrib, ndata
!cggg snow mods end
!
      character*500 fntsfa,fnweta,fnsnoa,fnzora,fnalba,fnaisa,
     &             fntg3a,fnscva,fnsmca,fnstca,fnacna,fnvega,
     &             fnveta,fnsota
!clu [+1l] add fn()a for vmn, vmx, slp, abs
     &,            fnvmna,fnvmxa,fnslpa,fnabsa

      real (kind=kind_io8) tsfanl(len), wetanl(len),   snoanl(len),
     &     zoranl(len), albanl(len,4), aisanl(len),
     &     tg3anl(len), acnanl(len),
     &     cvanl (len), cvbanl(len),   cvtanl(len),
     &     slianl(len), scvanl(len),   veganl(len),
     &     vetanl(len), sotanl(len),   alfanl(len,2),
     &     smcanl(len,lsoil), stcanl(len,lsoil),
     &     tsfan0(len)
!clu [+1l] add ()anl for vmn, vmx, slp, abs
     &,    vmnanl(len),vmxanl(len),slpanl(len),absanl(len)
!
      logical gaus
!
! tsf
!
      irttsf = 1
      if(fntsfa(1:8).ne.'        ') then
        call fixrda(lugb,fntsfa,kpdtsf,slmask,
     &             iy,im,id,ih,fh,tsfanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irttsf = iret
        if(iret.eq.1) then
          write(6,*) 't surface analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old t surface analysis provided, indicating proper',
     &            ' file name is given.  no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'t surface analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no tsf analysis available.  climatology used'
        endif
      endif
!
! tsf0
!
!     if(fntsfa(1:8).ne.'        ') then
!       call fixrda(lugb,fntsfa,kpdtsf,slmask,
!    &             iy,im,id,ih,0.,tsfan0,len,iret
!    &,            imsk, jmsk, slmskh, gaus,blno, blto
!    &,            outlat, outlon, me)
!       if(iret.eq.1) then
!         write(6,*) 't surface at ft=0 analysis read error'
!         call abort
!       elseif(iret.eq.-1) then
!         write(6,*) 'could not find t surface analysis at ft=0'
!         call abort
!       else
!         print *,'t surface analysis at ft=0 found.'
!       endif
!     else
!       do i=1,len
!         tsfan0(i)=-999.9
!       enddo
!     endif
!
!  albedo
!
      irtalb=0
      if(fnalba(1:8).ne.'        ') then
        do kk = 1, 4
          call fixrda(lugb,fnalba,kpdalb(kk),slmask,
     &               iy,im,id,ih,fh,albanl(1,kk),len,iret
     &,              imsk, jmsk, slmskh, gaus,blno, blto
     &,              outlat, outlon, me)
          irtalb=iret
          if(iret.eq.1) then
            write(6,*) 'albedo analysis read error'
            call abort
          elseif(iret.eq.-1) then
            if (me .eq. 0) then
            print *,'old albedo analysis provided, indicating proper',
     &              ' file name is given.  no error suspected.'
            write(6,*) 'forecast guess will be used'
            endif
          else
            if (me .eq. 0 .and. kk .eq. 4)
     &                  print *,'albedo analysis provided.'
          endif
        enddo
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no albedo analysis available.  climatology used'
        endif
      endif
!
!  vegetation fraction for albedo
!
      irtalf=0
      if(fnalba(1:8).ne.'        ') then
        do kk = 1, 2
          call fixrda(lugb,fnalba,kpdalf(kk),slmask,
     &               iy,im,id,ih,fh,alfanl(1,kk),len,iret
     &,              imsk, jmsk, slmskh, gaus,blno, blto
     &,              outlat, outlon, me)
          irtalf=iret
          if(iret.eq.1) then
            write(6,*) 'albedo analysis read error'
            call abort
          elseif(iret.eq.-1) then
            if (me .eq. 0) then
            print *,'old albedo analysis provided, indicating proper',
     &              ' file name is given.  no error suspected.'
            write(6,*) 'forecast guess will be used'
            endif
          else
            if (me .eq. 0 .and. kk .eq. 4)
     &                  print *,'albedo analysis provided.'
          endif
        enddo
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no vegfalbedo analysis available.  climatology used'
        endif
      endif
!
!  soil wetness
!
      irtwet=0
      irtsmc=0
      if(fnweta(1:8).ne.'        ') then
        call fixrda(lugb,fnweta,kpdwet,slmask,
     &             iy,im,id,ih,fh,wetanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtwet=iret
        if(iret.eq.1) then
          write(6,*) 'bucket wetness analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old wetness analysis provided, indicating proper',
     &            ' file name is given.  no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'bucket wetness analysis provided.'
        endif
      elseif(fnsmca(1:8).ne.'        ') then
        call fixrda(lugb,fnsmca,kpdsmc,slmask,
     &             iy,im,id,ih,fh,smcanl(1,1),len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        call fixrda(lugb,fnsmca,kpdsmc,slmask,
     &             iy,im,id,ih,fh,smcanl(1,2),len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtsmc=iret
        if(iret.eq.1) then
          write(6,*) 'layer soil wetness analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old layer soil wetness analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'layer soil wetness analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no soil wetness analysis available.  climatology used'
        endif
      endif
!
!  read in snow depth/snow cover
!
      irtscv=0
      if(fnsnoa(1:8).ne.'        ') then
        do i=1,len
          scvanl(i)=0.
        enddo
!cggg snow mods start
!cggg need to determine if the snow data is on the gaussian grid
!cggg or not.  if gaussian, then data is a depth, not liq equiv
!cggg depth. if not gaussian, then data is from hua-lu's
!cggg program and is a liquid equiv.  need to communicate
!cggg this to routine fixrda via the 3rd argument which is
!cggg the grib parameter id number.
        call baopenr(lugb,fnsnoa,iret)
        if (iret .ne. 0) then
          write(6,*) ' error in opening file ',trim(fnsnoa)
          print *,'error in opening file ',trim(fnsnoa)
          call abort
        endif
        lugi=0
        lskip=-1
        jpds=-1
        jgds=-1
        kpds=jpds
        call getgbh(lugb,lugi,lskip,jpds,jgds,lgrib,ndata,
     &              lskip,kpds,kgds,iret)
        close(lugb)
        if (iret .ne. 0) then
          write(6,*) ' error reading header of file: ',trim(fnsnoa)
          print *,'error reading header of file: ',trim(fnsnoa)
          call abort
        endif
        if (kgds(1) == 4) then  ! gaussian data is depth
          call fixrda(lugb,fnsnoa,kpdsnd,slmask,
     &                iy,im,id,ih,fh,snoanl,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
          snoanl=snoanl*100.  ! convert from meters to liq. eq.
                              ! depth in mm using 10:1 ratio
        else                    ! lat/lon data is liq equv. depth
          call fixrda(lugb,fnsnoa,kpdsno,slmask,
     &                iy,im,id,ih,fh,snoanl,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
        endif
!cggg snow mods end
        irtscv=iret
        if(iret.eq.1) then
          write(6,*) 'snow depth analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old snow depth analysis provided, indicating proper',
     &            ' file name is given.  no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'snow depth analysis provided.'
        endif
        irtsno=0
      elseif(fnscva(1:8).ne.'        ') then
        do i=1,len
          snoanl(i)=0.
        enddo
        call fixrda(lugb,fnscva,kpdscv,slmask,
     &             iy,im,id,ih,fh,scvanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtsno=iret
        if(iret.eq.1) then
          write(6,*) 'snow cover analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old snow cover analysis provided, indicating proper',
     &            ' file name is given.  no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'snow cover analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no snow/snocov analysis available.  climatology used'
        endif
      endif
!
!  sea ice mask
!
      irtacn=0
      irtais=0
      if(fnacna(1:8).ne.'        ') then
        call fixrda(lugb,fnacna,kpdacn,slmask,
     &             iy,im,id,ih,fh,acnanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtacn=iret
        if(iret.eq.1) then
          write(6,*) 'ice concentration analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old ice concentration analysis provided',
     &            ' indicating proper file name is given'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'ice concentration analysis provided.'
        endif
      elseif(fnaisa(1:8).ne.'        ') then
        call fixrda(lugb,fnaisa,kpdais,slmask,
     &             iy,im,id,ih,fh,aisanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtais=iret
        if(iret.eq.1) then
          write(6,*) 'ice mask analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old ice-mask analysis provided, indicating proper',
     &            ' file name is given.  no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'ice mask analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no sea-ice analysis available.  climatology used'
        endif
      endif
!
!  surface roughness
!
      irtzor=0
      if(fnzora(1:8).ne.'        ') then
        call fixrda(lugb,fnzora,kpdzor,slmask,
     &             iy,im,id,ih,fh,zoranl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtzor=iret
        if(iret.eq.1) then
          write(6,*) 'roughness analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old roughness analysis provided, indicating proper',
     &            ' file name is given.  no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'roughness analysis provided.'
        endif
      else
          if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no srfc roughness analysis available. climatology used'
        endif
      endif
!
!  deep soil temperature
!
      irttg3=0
      irtstc=0
      if(fntg3a(1:8).ne.'        ') then
        call fixrda(lugb,fntg3a,kpdtg3,slmask,
     &             iy,im,id,ih,fh,tg3anl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irttg3=iret
        if(iret.eq.1) then
          write(6,*) 'deep soil tmp analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old deep soil temp analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'deep soil tmp analysis provided.'
        endif
      elseif(fnstca(1:8).ne.'        ') then
        call fixrda(lugb,fnstca,kpdstc,slmask,
     &             iy,im,id,ih,fh,stcanl(1,1),len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        call fixrda(lugb,fnstca,kpdstc,slmask,
     &             iy,im,id,ih,fh,stcanl(1,2),len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtstc=iret
        if(iret.eq.1) then
          write(6,*) 'layer soil tmp analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old deep soil temp analysis provided',
     &            'iindicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'layer soil tmp analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no deep soil temp analy available.  climatology used'
        endif
      endif
!
!  vegetation cover
!
      irtveg=0
      if(fnvega(1:8).ne.'        ') then
        call fixrda(lugb,fnvega,kpdveg,slmask,
     &             iy,im,id,ih,fh,veganl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtveg=iret
        if(iret.eq.1) then
          write(6,*) 'vegetation cover analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old vegetation cover analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'gegetation cover analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no vegetation cover anly available. climatology used'
        endif
      endif
!
!  vegetation type
!
      irtvet=0
      if(fnveta(1:8).ne.'        ') then
        call fixrda(lugb,fnveta,kpdvet,slmask,
     &             iy,im,id,ih,fh,vetanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtvet=iret
        if(iret.eq.1) then
          write(6,*) 'vegetation type analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old vegetation type analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'vegetation type analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no vegetation type anly available. climatology used'
        endif
      endif
!
!  soil type
!
      irtsot=0
      if(fnsota(1:8).ne.'        ') then
        call fixrda(lugb,fnsota,kpdsot,slmask,
     &             iy,im,id,ih,fh,sotanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtsot=iret
        if(iret.eq.1) then
          write(6,*) 'soil type analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old soil type analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'soil type analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no soil type anly available. climatology used'
        endif
      endif

!clu [+120l]--------------------------------------------------------------
!
!  min vegetation cover
!
      irtvmn=0
      if(fnvmna(1:8).ne.'        ') then
        call fixrda(lugb,fnvmna,kpdvmn,slmask,
     &             iy,im,id,ih,fh,vmnanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtvmn=iret
        if(iret.eq.1) then
          write(6,*) 'shdmin analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old shdmin analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'shdmin analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no shdmin anly available. climatology used'
        endif
      endif

!
!  max vegetation cover
!
      irtvmx=0
      if(fnvmxa(1:8).ne.'        ') then
        call fixrda(lugb,fnvmxa,kpdvmx,slmask,
     &             iy,im,id,ih,fh,vmxanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtvmx=iret
        if(iret.eq.1) then
          write(6,*) 'shdmax analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old shdmax analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'shdmax analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no shdmax anly available. climatology used'
        endif
      endif

!
!  slope type
!
      irtslp=0
      if(fnslpa(1:8).ne.'        ') then
        call fixrda(lugb,fnslpa,kpdslp,slmask,
     &             iy,im,id,ih,fh,slpanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtslp=iret
        if(iret.eq.1) then
          write(6,*) 'slope type analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old slope type analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'slope type analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no slope type anly available. climatology used'
        endif
      endif

!
!  max snow albedo
!
      irtabs=0
      if(fnabsa(1:8).ne.'        ') then
        call fixrda(lugb,fnabsa,kpdabs,slmask,
     &             iy,im,id,ih,fh,absanl,len,iret
     &,            imsk, jmsk, slmskh, gaus,blno, blto
     &,            outlat, outlon, me)
        irtabs=iret
        if(iret.eq.1) then
          write(6,*) 'snoalb analysis read error'
          call abort
        elseif(iret.eq.-1) then
          if (me .eq. 0) then
          print *,'old snoalb analysis provided',
     &            ' indicating proper file name is given.'
          print *,' no error suspected.'
          write(6,*) 'forecast guess will be used'
          endif
        else
          if (me .eq. 0) print *,'snoalb analysis provided.'
        endif
      else
        if (me .eq. 0) then
!       print *,'************************************************'
        print *,'no snoalb anly available. climatology used'
        endif
      endif

!clu ----------------------------------------------------------------------
!
      return
      end
      subroutine filfcs(tsffcs,wetfcs,snofcs,zorfcs,albfcs,
     &                  tg3fcs,cvfcs ,cvbfcs,cvtfcs,
     &                  cnpfcs,smcfcs,stcfcs,slifcs,aisfcs,
     &                  vegfcs, vetfcs, sotfcs, alffcs,
!cwu [+1l] add ()fcs for sih, sic
     &                  sihfcs,sicfcs,
!clu [+1l] add ()fcs for vmn, vmx, slp, abs
     &                  vmnfcs,vmxfcs,slpfcs,absfcs,
     &                  tsfanl,wetanl,snoanl,zoranl,albanl,
     &                  tg3anl,cvanl ,cvbanl,cvtanl,
     &                  cnpanl,smcanl,stcanl,slianl,aisanl,
     &                  veganl, vetanl, sotanl, alfanl,
!cwu [+1l] add ()anl for sih, sic
     &                  sihanl,sicanl,
!clu [+1l] add ()anl for vmn, vmx, slp, abs
     &                  vmnanl,vmxanl,slpanl,absanl,
     &                  len,lsoil)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,j,len,lsoil
      real (kind=kind_io8) tsffcs(len),wetfcs(len),snofcs(len),
     &     zorfcs(len),albfcs(len,4),aisfcs(len),
     &     tg3fcs(len),
     &     cvfcs (len),cvbfcs(len),cvtfcs(len),
     &     cnpfcs(len),
     &     smcfcs(len,lsoil),stcfcs(len,lsoil),
     &     slifcs(len),vegfcs(len),
     &     vetfcs(len),sotfcs(len),alffcs(len,2)
!cwu [+1l] add ()fcs for sih, sic
     &,    sihfcs(len),sicfcs(len)
!clu [+1l] add ()fcs for vmn, vmx, slp, abs
     &,    vmnfcs(len),vmxfcs(len),slpfcs(len),absfcs(len)
      real (kind=kind_io8) tsfanl(len),wetanl(len),snoanl(len),
     &     zoranl(len),albanl(len,4),aisanl(len),
     &     tg3anl(len),
     &     cvanl (len),cvbanl(len),cvtanl(len),
     &     cnpanl(len),
     &     smcanl(len,lsoil),stcanl(len,lsoil),
     &     slianl(len),veganl(len),
     &     vetanl(len),sotanl(len),alfanl(len,2)
!cwu [+1l] add ()anl for sih, sic
     &,    sihanl(len),sicanl(len)
!clu [+1l] add ()anl for vmn, vmx, slp, abs
     &,    vmnanl(len),vmxanl(len),slpanl(len),absanl(len)
!
      write(6,*) '  this is a dead start run, tsfc over land is',
     &           ' set as lowest sigma level temperture if given.'
      write(6,*) '  if not, set to climatological tsf over land is used'
!
!
      do i=1,len
        tsffcs(i)   = tsfanl(i)      !  tsf
        albfcs(i,1) = albanl(i,1)    !  albedo
        albfcs(i,2) = albanl(i,2)    !  albedo
        albfcs(i,3) = albanl(i,3)    !  albedo
        albfcs(i,4) = albanl(i,4)    !  albedo
        wetfcs(i)   = wetanl(i)      !  soil wetness
        snofcs(i)   = snoanl(i)      !  snow
        aisfcs(i)   = aisanl(i)      !  seaice
        slifcs(i)   = slianl(i)      !  land/sea/snow mask
        zorfcs(i)   = zoranl(i)      !  surface roughness
!       plrfcs(i)   = plranl(i)      !  maximum stomatal resistance
        tg3fcs(i)   = tg3anl(i)      !  deep soil temperature
        cnpfcs(i)   = cnpanl(i)      !  canopy water content
        cvfcs(i)    = cvanl(i)       !  cv
        cvbfcs(i)   = cvbanl(i)      !  cvb
        cvtfcs(i)   = cvtanl(i)      !  cvt
        vegfcs(i)   = veganl(i)      !  vegetation cover
        vetfcs(i)   = vetanl(i)      !  vegetation type
        sotfcs(i)   = sotanl(i)      !  soil type
        alffcs(i,1) = alfanl(i,1)    !  vegetation fraction for albedo
        alffcs(i,2) = alfanl(i,2)    !  vegetation fraction for albedo
!cwu [+2l] add sih, sic
        sihfcs(i)   = sihanl(i)      !  sea ice thickness
        sicfcs(i)   = sicanl(i)      !  sea ice concentration
!clu [+4l] add vmn, vmx, slp, abs
        vmnfcs(i)   = vmnanl(i)      !  min vegetation cover
        vmxfcs(i)   = vmxanl(i)      !  max vegetation cover
        slpfcs(i)   = slpanl(i)      !  slope type
        absfcs(i)   = absanl(i)      !  max snow albedo
      enddo
!
      do j=1,lsoil
        do i=1,len
          smcfcs(i,j) = smcanl(i,j)  !   layer soil wetness
          stcfcs(i,j) = stcanl(i,j)  !   soil temperature
        enddo
      enddo
!
      return
      end
      subroutine bktges(smcfcs,slianl,stcfcs,len,lsoil)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,j,len,lsoil,k
      real (kind=kind_io8) smcfcs(len,lsoil), stcfcs(len,lsoil),
     &                     slianl(len)
!
!  note that smfcs comes in with the original unit (cm?) (not grib file)
!
      do i = 1, len
        smcfcs(i,1) = (smcfcs(i,1)/150.) * .37 + .1
      enddo
      do k = 2, lsoil
        do i = 1, len
          smcfcs(i,k) = smcfcs(i,1)
        enddo
      enddo
      if(lsoil.gt.2) then
        do k = 3, lsoil
          do i = 1, len
            stcfcs(i,k) = stcfcs(i,2)
          enddo
        enddo
      endif
!
      return
      end
      subroutine rof01(aisfld,len,op,crit)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) aisfld(len),crit
      character*2 op
!
      if(op.eq.'ge') then
        do i=1,len
          if(aisfld(i).ge.crit) then
            aisfld(i)=1.
          else
            aisfld(i)=0.
          endif
        enddo
      elseif(op.eq.'gt') then
        do i=1,len
          if(aisfld(i).gt.crit) then
            aisfld(i)=1.
          else
            aisfld(i)=0.
          endif
        enddo
      elseif(op.eq.'le') then
        do i=1,len
          if(aisfld(i).le.crit) then
            aisfld(i)=1.
          else
            aisfld(i)=0.
          endif
        enddo
      elseif(op.eq.'lt') then
        do i=1,len
          if(aisfld(i).lt.crit) then
            aisfld(i)=1.
          else
            aisfld(i)=0.
          endif
        enddo
      else
        write(6,*) ' illegal operator in rof01.  op=',op
        call abort
      endif
!
      return
      end
      subroutine tsfcor(tsfc,orog,slmask,umask,len,rlapse)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) rlapse,umask
      real (kind=kind_io8) tsfc(len), orog(len), slmask(len)
!
      do i=1,len
        if(slmask(i).eq.umask) then
          tsfc(i) = tsfc(i) - orog(i)*rlapse
        endif
      enddo
      return
      end
      subroutine snodpth(scvanl,slianl,tsfanl,snoclm,
     &                   glacir,snwmax,snwmin,landice,len,snoanl, me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,me,len
      logical, intent(in) :: landice
      real (kind=kind_io8) sno,snwmax,snwmin
!
      real (kind=kind_io8) scvanl(len), slianl(len), tsfanl(len),
     &     snoclm(len), snoanl(len), glacir(len)
!
      if (me .eq. 0) write(6,*) 'snodpth'
!
!  use surface temperature to get snow depth estimate
!
      do i=1,len
        sno = 0.0
!
!  over land
!
        if(slianl(i).eq.1.) then
          if(scvanl(i).eq.1.0) then
            if(tsfanl(i).lt.243.0) then
              sno = snwmax
            elseif(tsfanl(i).lt.273.0) then
              sno = snwmin+(snwmax-snwmin)*(273.0-tsfanl(i))/30.0
            else
              sno = snwmin
            endif
          endif
!
!  if glacial points has snow in climatology, set sno to snomax
!
          if (.not.landice) then
            if(glacir(i).eq.1.0) then
              sno = snoclm(i)
              if(sno.eq.0.) sno=snwmax
            endif
          endif
        endif
!
!  over sea ice
!
!  snow over sea ice is cycled as of 01/01/94.....hua-lu pan
!
        if(slianl(i).eq.2.0) then
          sno=snoclm(i)
          if(sno.eq.0.) sno=snwmax
        endif
!
        snoanl(i) = sno
      enddo
      return
      end subroutine snodpth
      subroutine merge(len,lsoil,iy,im,id,ih,fh,deltsfc,
     &                 sihfcs,sicfcs,
     &                 vmnfcs,vmxfcs,slpfcs,absfcs,
     &                 tsffcs,wetfcs,snofcs,zorfcs,albfcs,aisfcs,
     &                 cvfcs ,cvbfcs,cvtfcs,
     &                 cnpfcs,smcfcs,stcfcs,slifcs,vegfcs,
     &                 vetfcs,sotfcs,alffcs,
     &                 sihanl,sicanl,                 
     &                 vmnanl,vmxanl,slpanl,absanl,
     &                 tsfanl,tsfan2,wetanl,snoanl,zoranl,albanl,aisanl,
     &                 cvanl ,cvbanl,cvtanl,
     &                 cnpanl,smcanl,stcanl,slianl,veganl,
     &                 vetanl,sotanl,alfanl,
     &                 ctsfl,calbl,caisl,csnol,csmcl,czorl,cstcl,cvegl,
     &                 ctsfs,calbs,caiss,csnos,csmcs,czors,cstcs,cvegs,
     &                 ccv,ccvb,ccvt,ccnp,cvetl,cvets,csotl,csots,
     &                 calfl,calfs,
     &                 csihl,csihs,csicl,csics,
     &                 cvmnl,cvmns,cvmxl,cvmxs,cslpl,cslps,cabsl,cabss,
     &                 irttsf,irtwet,irtsno,irtzor,irtalb,irtais,
     &                 irttg3,irtscv,irtacn,irtsmc,irtstc,irtveg,
     &                 irtvmn,irtvmx,irtslp,irtabs,
     &                 irtvet,irtsot,irtalf, landice, me)
      use machine , only : kind_io8,kind_io4
      use sfccyc_module, only : veg_type_landice, soil_type_landice
      implicit none
      integer k,i,im,id,iy,len,lsoil,ih,irtacn,irtsmc,irtscv,irtais,
     &        irttg3,irtstc,irtalf,me,irtsot,irtveg,irtvet, irtzor,
     &        irtalb,irtsno,irttsf,irtwet,j
     &,       irtvmn,irtvmx,irtslp,irtabs
      logical, intent(in)  :: landice
      real (kind=kind_io8) rvegs,rvets,rzors,raiss,rsnos,rsots,rcnp,
     &                     rcvt,rcv,rcvb,rsnol,rzorl,raisl,ralbl,
     &                     ralfl,rvegl,ralbs,ralfs,rtsfs,rvetl,rsotl,
     &                     qzors,qvegs,qsnos,qalfs,qaiss,qvets,qcvt,
     &                     qcnp,qcvb,qsots,qcv,qaisl,qsnol,qalfl,
     &                     qtsfl,qalbl,qzorl,qtsfs,qalbs,qsotl,qvegl,
     &                     qvetl,rtsfl,calbs,caiss,ctsfs,czorl,cvegl,
     &                     csnos,ccvb,ccvt,ccv,czors,cvegs,caisl,csnol,
     &                     calbl,fh,ctsfl,ccnp,csots,calfl,csotl,cvetl,
     &                     cvets,calfs,deltsfc,
     &                     csihl,csihs,csicl,csics,
     &                     rsihl,rsihs,rsicl,rsics,
     &                     qsihl,qsihs,qsicl,qsics
     &,                    cvmnl,cvmns,cvmxl,cvmxs,cslpl,cslps
     &,                    cabsl,cabss,rvmnl,rvmns,rvmxl,rvmxs
     &,                    rslpl,rslps,rabsl,rabss,qvmnl,qvmns
     &,                    qvmxl,qvmxs,qslpl,qslps,qabsl,qabss
!
      real (kind=kind_io8) tsffcs(len), wetfcs(len),   snofcs(len),
     &     zorfcs(len), albfcs(len,4), aisfcs(len),
     &     cvfcs (len), cvbfcs(len),   cvtfcs(len),
     &     cnpfcs(len),
     &     smcfcs(len,lsoil),stcfcs(len,lsoil),
     &     slifcs(len), vegfcs(len),
     &     vetfcs(len), sotfcs(len),   alffcs(len,2)
     &,    sihfcs(len), sicfcs(len)
     &,    vmnfcs(len),vmxfcs(len),slpfcs(len),absfcs(len)
      real (kind=kind_io8) tsfanl(len),tsfan2(len),
     &     wetanl(len),snoanl(len),
     &     zoranl(len), albanl(len,4), aisanl(len),
     &     cvanl (len), cvbanl(len),   cvtanl(len),
     &     cnpanl(len),
     &     smcanl(len,lsoil),stcanl(len,lsoil),
     &     slianl(len), veganl(len),
     &     vetanl(len), sotanl(len),   alfanl(len,2)
     &,    sihanl(len),sicanl(len)           
     &,    vmnanl(len),vmxanl(len),slpanl(len),absanl(len)
!
      real (kind=kind_io8) csmcl(lsoil), csmcs(lsoil),
     &                     cstcl(lsoil), cstcs(lsoil)
      real (kind=kind_io8) rsmcl(lsoil), rsmcs(lsoil),
     &                     rstcl(lsoil), rstcs(lsoil)
      real (kind=kind_io8) qsmcl(lsoil), qsmcs(lsoil),
     &                     qstcl(lsoil), qstcs(lsoil)
      logical first
      integer   num_threads
      data first /.true./
      save num_threads, first
!
      integer len_thread_m, i1_t, i2_t, it
      integer num_parthds
!
      if (first) then
         num_threads = num_parthds()
         first = .false.
      endif
!
!  coeeficients of blending forecast and interpolated clim
!  (or analyzed) fields over sea or land(l) (not for clouds)
!  1.0 = use of forecast
!  0.0 = replace with interpolated analysis
!
!  merging coefficients are defined by parameter statement in calling program
!  and therefore they should not be modified in this program.
!
      rtsfl = ctsfl
      ralbl = calbl
      ralfl = calfl
      raisl = caisl
      rsnol = csnol
!clu  rsmcl = csmcl
      rzorl = czorl
      rvegl = cvegl
      rvetl = cvetl
      rsotl = csotl
      rsihl = csihl
      rsicl = csicl
      rvmnl = cvmnl
      rvmxl = cvmxl
      rslpl = cslpl
      rabsl = cabsl
!
      rtsfs = ctsfs
      ralbs = calbs
      ralfs = calfs
      raiss = caiss
      rsnos = csnos
!     rsmcs = csmcs
      rzors = czors
      rvegs = cvegs
      rvets = cvets
      rsots = csots
      rsihs = csihs
      rsics = csics
      rvmns = cvmns
      rvmxs = cvmxs
      rslps = cslps
      rabss = cabss
!
      rcv  = ccv
      rcvb = ccvb
      rcvt = ccvt
      rcnp = ccnp
!
      do k=1,lsoil
        rsmcl(k) = csmcl(k)
        rsmcs(k) = csmcs(k)
        rstcl(k) = cstcl(k)
        rstcs(k) = cstcs(k)
      enddo
      if (fh-deltsfc < -0.001 .and. irttsf == 1) then
        rtsfs = 1.0
        rtsfl = 1.0
!       do k=1,lsoil
!         rsmcl(k) = 1.0
!         rsmcs(k) = 1.0
!         rstcl(k) = 1.0
!         rstcs(k) = 1.0
!       enddo
      endif
!
!  if analysis file name is given but no matching analysis date found,
!  use guess (these are flagged by irt???=1).
!
      if(irttsf == -1) then
        rtsfl = 1.
        rtsfs = 1.
      endif
      if(irtalb == -1) then
        ralbl = 1.
        ralbs = 1.
        ralfl = 1.
        ralfs = 1.
      endif
      if(irtais == -1) then
        raisl = 1.
        raiss = 1.
      endif
      if(irtsno == -1 .or. irtscv == -1) then
        rsnol = 1.
        rsnos = 1.
      endif
      if(irtsmc == -1 .or. irtwet == -1) then
!       rsmcl = 1.
!       rsmcs = 1.
        do k=1,lsoil
          rsmcl(k) = 1.
          rsmcs(k) = 1.
        enddo
      endif
      if(irtstc.eq.-1) then
        do k=1,lsoil
          rstcl(k) = 1.
          rstcs(k) = 1.
        enddo
      endif
      if(irtzor == -1) then
        rzorl = 1.
        rzors = 1.
      endif
      if(irtveg == -1) then
        rvegl = 1.
        rvegs = 1.
      endif
      if(irtvet.eq.-1) then
        rvetl = 1.
        rvets = 1.
      endif
      if(irtsot == -1) then
        rsotl = 1.
        rsots = 1.
      endif

      if(irtacn == -1) then
        rsicl = 1.
        rsics = 1.
      endif
      if(irtvmn == -1) then
        rvmnl = 1.
        rvmns = 1.
      endif
      if(irtvmx == -1) then
        rvmxl = 1.
        rvmxs = 1.
      endif
      if(irtslp == -1) then
        rslpl = 1.
        rslps = 1.
      endif
      if(irtabs == -1) then
        rabsl = 1.
        rabss = 1.
      endif
!
      if(raiss == 1. .or. irtacn == -1) then
        if (me == 0) print *,'use forecast land-sea-ice mask'
        do i = 1, len
          aisanl(i) = aisfcs(i)
          slianl(i) = slifcs(i)
        enddo
      endif
!
      if (me == 0) then
      write(6,100) rtsfl,ralbl,raisl,rsnol,rsmcl,rzorl,rvegl
  100 format('rtsfl,ralbl,raisl,rsnol,rsmcl,rzorl,rvegl=',10f7.3)
      write(6,101) rtsfs,ralbs,raiss,rsnos,rsmcs,rzors,rvegs
  101 format('rtsfs,ralbs,raiss,rsnos,rsmcs,rzors,rvegs=',10f7.3)
!     print *,' ralfl=',ralfl,' ralfs=',ralfs,' rsotl=',rsotl
!    *,' rsots=',rsots,' rvetl=',rvetl,' rvets=',rvets
      endif
!
      qtsfl = 1. - rtsfl
      qalbl = 1. - ralbl
      qalfl = 1. - ralfl
      qaisl = 1. - raisl
      qsnol = 1. - rsnol
!     qsmcl = 1. - rsmcl
      qzorl = 1. - rzorl
      qvegl = 1. - rvegl
      qvetl = 1. - rvetl
      qsotl = 1. - rsotl
      qsihl = 1. - rsihl
      qsicl = 1. - rsicl
      qvmnl = 1. - rvmnl
      qvmxl = 1. - rvmxl
      qslpl = 1. - rslpl
      qabsl = 1. - rabsl
!
      qtsfs = 1. - rtsfs
      qalbs = 1. - ralbs
      qalfs = 1. - ralfs
      qaiss = 1. - raiss
      qsnos = 1. - rsnos
!     qsmcs = 1. - rsmcs
      qzors = 1. - rzors
      qvegs = 1. - rvegs
      qvets = 1. - rvets
      qsots = 1. - rsots
      qsihs = 1. - rsihs
      qsics = 1. - rsics
      qvmns = 1. - rvmns
      qvmxs = 1. - rvmxs
      qslps = 1. - rslps
      qabss = 1. - rabss
!
      qcv   = 1. - rcv
      qcvb  = 1. - rcvb
      qcvt  = 1. - rcvt
      qcnp  = 1. - rcnp
!
      do k=1,lsoil
        qsmcl(k) = 1. - rsmcl(k)
        qsmcs(k) = 1. - rsmcs(k)
        qstcl(k) = 1. - rstcl(k)
        qstcs(k) = 1. - rstcs(k)
      enddo
!
!  merging
!
      if(me .eq. 0) then
        print *, 'dbgx-- csmcl:', (csmcl(k),k=1,lsoil)
        print *, 'dbgx-- rsmcl:', (rsmcl(k),k=1,lsoil)
        print *, 'dbgx-- csnol, csnos:',csnol,csnos
        print *, 'dbgx-- rsnol, rsnos:',rsnol,rsnos
      endif

!     print *, rtsfs, qtsfs, raiss , qaiss
!    *,        rsnos , qsnos, rzors , qzors, rvegs , qvegs
!    *,        rvets , qvets, rsots , qsots
!    *,        rcv, rcvb, rcvt, qcv, qcvb, qcvt
!    *,        ralbs, qalbs, ralfs, qalfs
!     print *, rtsfl, qtsfl, raisl , qaisl
!    *,        rsnol , qsnol, rzorl , qzorl, rvegl , qvegl
!    *,        rvetl , qvetl, rsotl , qsotl
!    *,        ralbl, qalbl, ralfl, qalfl
!
!
      len_thread_m  = (len+num_threads-1) / num_threads

!$omp parallel do private(i1_t,i2_t,it,i)
      do it=1,num_threads   ! start of threaded loop ...................
        i1_t       = (it-1)*len_thread_m+1
        i2_t       = min(i1_t+len_thread_m-1,len)
        do i=i1_t,i2_t
          if(slianl(i).eq.0.) then
            vetanl(i) = vetfcs(i)*rvets + vetanl(i)*qvets
            sotanl(i) = sotfcs(i)*rsots + sotanl(i)*qsots
          else
            vetanl(i) = vetfcs(i)*rvetl + vetanl(i)*qvetl
            sotanl(i) = sotfcs(i)*rsotl + sotanl(i)*qsotl
          endif
        enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(i1_t,i2_t,it,i,k)
!
      do it=1,num_threads   ! start of threaded loop ...................
        i1_t       = (it-1)*len_thread_m+1
        i2_t       = min(i1_t+len_thread_m-1,len)
!
      do i=i1_t,i2_t
        if(slianl(i).eq.0.) then
!.... tsffc2 is the previous anomaly + today's climatology
!         tsffc2 = (tsffcs(i)-tsfan2(i))+tsfanl(i)
!         tsfanl(i) = tsffc2    *rtsfs+tsfanl(i)*qtsfs
!
          tsfanl(i) = tsffcs(i)*rtsfs + tsfanl(i)*qtsfs
!         albanl(i) = albfcs(i)*ralbs + albanl(i)*qalbs
          aisanl(i) = aisfcs(i)*raiss + aisanl(i)*qaiss
          snoanl(i) = snofcs(i)*rsnos + snoanl(i)*qsnos
          
          zoranl(i) = zorfcs(i)*rzors + zoranl(i)*qzors
          veganl(i) = vegfcs(i)*rvegs + veganl(i)*qvegs
          sihanl(i) = sihfcs(i)*rsihs + sihanl(i)*qsihs
          sicanl(i) = sicfcs(i)*rsics + sicanl(i)*qsics
          vmnanl(i) = vmnfcs(i)*rvmns + vmnanl(i)*qvmns
          vmxanl(i) = vmxfcs(i)*rvmxs + vmxanl(i)*qvmxs
          slpanl(i) = slpfcs(i)*rslps + slpanl(i)*qslps
          absanl(i) = absfcs(i)*rabss + absanl(i)*qabss
        else
          tsfanl(i) = tsffcs(i)*rtsfl + tsfanl(i)*qtsfl
!         albanl(i) = albfcs(i)*ralbl + albanl(i)*qalbl
          aisanl(i) = aisfcs(i)*raisl + aisanl(i)*qaisl
          if(rsnol.ge.0)then
            snoanl(i) = snofcs(i)*rsnol + snoanl(i)*qsnol
          else  ! envelope method
            if(snoanl(i).ne.0)then
             snoanl(i) = max(-snoanl(i)/rsnol,
     &                   min(-snoanl(i)*rsnol, snofcs(i)))
            endif
          endif
          zoranl(i) = zorfcs(i)*rzorl + zoranl(i)*qzorl
          veganl(i) = vegfcs(i)*rvegl + veganl(i)*qvegl
          vmnanl(i) = vmnfcs(i)*rvmnl + vmnanl(i)*qvmnl
          vmxanl(i) = vmxfcs(i)*rvmxl + vmxanl(i)*qvmxl
          slpanl(i) = slpfcs(i)*rslpl + slpanl(i)*qslpl
          absanl(i) = absfcs(i)*rabsl + absanl(i)*qabsl
          sihanl(i) = sihfcs(i)*rsihl + sihanl(i)*qsihl
          sicanl(i) = sicfcs(i)*rsicl + sicanl(i)*qsicl
        endif

        cnpanl(i) = cnpfcs(i)*rcnp + cnpanl(i)*qcnp
!
!  snow over sea ice is cycled
!
        if(slianl(i).eq.2.) then
          snoanl(i) = snofcs(i)
        endif
!
      enddo

! at landice points, set the soil type, slope type and
! greenness fields to flag values.

      if (landice) then
        do i=i1_t,i2_t
          if (nint(slianl(i)) == 1) then
            if (nint(vetanl(i)) == veg_type_landice) then
              sotanl(i) = soil_type_landice
              veganl(i) = 0.0
              slpanl(i) = 9.0
              vmnanl(i) = 0.0
              vmxanl(i) = 0.0
            endif
          end if  ! if land
        enddo
      endif

      do i=i1_t,i2_t
        cvanl(i)  = cvfcs(i)*rcv   + cvanl(i)*qcv
        cvbanl(i) = cvbfcs(i)*rcvb + cvbanl(i)*qcvb
        cvtanl(i) = cvtfcs(i)*rcvt + cvtanl(i)*qcvt
      enddo
!
      do k = 1, 4
        do i=i1_t,i2_t
          if(slianl(i).eq.0.) then
            albanl(i,k) = albfcs(i,k)*ralbs + albanl(i,k)*qalbs
          else
            albanl(i,k) = albfcs(i,k)*ralbl + albanl(i,k)*qalbl
          endif
        enddo
      enddo
!
      do k = 1, 2
        do i=i1_t,i2_t
          if(slianl(i).eq.0.) then
            alfanl(i,k) = alffcs(i,k)*ralfs + alfanl(i,k)*qalfs
          else
            alfanl(i,k) = alffcs(i,k)*ralfl + alfanl(i,k)*qalfl
          endif
        enddo
      enddo
!
      do k = 1, lsoil
        do i=i1_t,i2_t
          if(slianl(i).eq.0.) then
            smcanl(i,k) = smcfcs(i,k)*rsmcs(k) + smcanl(i,k)*qsmcs(k)
            stcanl(i,k) = stcfcs(i,k)*rstcs(k) + stcanl(i,k)*qstcs(k)
          else
! soil moisture not used at landice points, so
! don't bother merging it.  also, for now don't allow nudging
! to raise subsurface temperature above freezing.
            stcanl(i,k) = stcfcs(i,k)*rstcl(k) + stcanl(i,k)*qstcl(k)
            if (landice .and. slianl(i) == 1.0 .and.
     &          nint(vetanl(i)) == veg_type_landice) then
              smcanl(i,k) = 1.0  ! use value as flag
              stcanl(i,k) = min(stcanl(i,k), 273.15)
            else
              smcanl(i,k) = smcfcs(i,k)*rsmcl(k) + smcanl(i,k)*qsmcl(k)
            end if
          endif
        enddo
      enddo
!
      enddo            ! end of threaded loop ...................
!$omp end parallel do
      return
      end subroutine merge
      subroutine newice(slianl,slifcs,tsfanl,tsffcs,len,lsoil,
!cwu [+1l] add sihnew,sicnew,sihanl,sicanl
     &                   sihnew,sicnew,sihanl,sicanl,      
     &                   albanl,snoanl,zoranl,smcanl,stcanl,
     &                   albsea,snosea,zorsea,smcsea,smcice,
     &                   tsfmin,tsfice,albice,zorice,tgice,
     &                   rla,rlo,me)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      real (kind=kind_io8), parameter :: one=1.0
      real (kind=kind_io8) tgice,albice,zorice,tsfice,albsea,snosea,
     &                     smcice,tsfmin,zorsea,smcsea
!cwu [+1l] add sicnew,sihnew
     &,                    sicnew,sihnew   
      integer i,me,kount1,kount2,k,len,lsoil
      real (kind=kind_io8) slianl(len),   slifcs(len),
     &                     tsffcs(len),tsfanl(len)
      real (kind=kind_io8) albanl(len,4), snoanl(len), zoranl(len)
      real (kind=kind_io8) smcanl(len,lsoil), stcanl(len,lsoil)
!cwu [+1l] add sihanl & sicanl
      real (kind=kind_io8) sihanl(len), sicanl(len)
!
      real (kind=kind_io8) rla(len), rlo(len)
!
      if (me .eq. 0) write(6,*) 'newice'
!
      kount1 = 0
      kount2 = 0
      do i=1,len
        if(slifcs(i).ne.slianl(i)) then
          if(slifcs(i).eq.1..or.slianl(i).eq.1.) then
            print *,'inconsistency in slifcs or slianl'
            print 910,rla(i),rlo(i),slifcs(i),slianl(i),
     &                tsffcs(i),tsfanl(i)
  910       format(2x,'at lat=',f5.1,' lon=',f5.1,' slifcs=',f4.1,
     &          ' slimsk=',f4.1,' tsffcs=',f5.1,' set to tsfanl=',f5.1)
            call abort
          endif
!
!  interpolated climatology indicates melted sea ice
!
          if(slianl(i).eq.0..and.slifcs(i).eq.2.) then
            tsfanl(i)   = tsfmin
            albanl(i,1) = albsea
            albanl(i,2) = albsea
            albanl(i,3) = albsea
            albanl(i,4) = albsea
            snoanl(i)   = snosea
            zoranl(i)   = zorsea
            do k = 1, lsoil
              smcanl(i,k) = smcsea
!cwu [+1l] set stcanl to tgice (over sea-ice)
              stcanl(i,k) = tgice 
            enddo
!cwu [+2l] set siganl and sicanl
            sihanl(i) = 0.
            sicanl(i) = 0.
            kount1 = kount1 + 1
          endif
!
!  interplated climatoloyg/analysis indicates new sea ice
!
          if(slianl(i).eq.2..and.slifcs(i).eq.0.) then
            tsfanl(i)   = tsfice
            albanl(i,1) = albice
            albanl(i,2) = albice
            albanl(i,3) = albice
            albanl(i,4) = albice
            snoanl(i)   = 0.
            zoranl(i)   = zorice
            do k = 1, lsoil
              smcanl(i,k) = smcice
              stcanl(i,k) = tgice
            enddo
!cwu [+2l] add sihanl & sicanl
            sihanl(i) = sihnew
            sicanl(i) = min(one, max(sicnew,sicanl(i)))
            kount2 = kount2 + 1
          endif
        endif
      enddo
!
      if (me .eq. 0) then
      if(kount1.gt.0) then
        write(6,*) 'sea ice melted.  tsf,alb,zor are filled',
     &             ' at ',kount1,' points'
      endif
      if(kount2.gt.0) then
        write(6,*) 'sea ice formed.  tsf,alb,zor are filled',
     &             ' at ',kount2,' points'
      endif
      endif
!
      return
      end
      subroutine qcsnow(snoanl,slmask,aisanl,glacir,len,snoval,
     &                  landice,me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer kount,i,len,me
      logical, intent(in)  :: landice
      real (kind=kind_io8) per,snoval
      real (kind=kind_io8) snoanl(len),slmask(len),
     &                     aisanl(len),glacir(len)
      if (me .eq. 0) then
        write(6,*) ' '
        write(6,*) 'qc of snow'
      endif
      if (.not.landice) then
        kount=0
        do i=1,len
          if(glacir(i).ne.0..and.snoanl(i).eq.0.) then
!         if(glacir(i).ne.0..and.snoanl(i).lt.snoval*0.5) then
            snoanl(i) = snoval
            kount     = kount + 1
          endif
        enddo
        per = float(kount) / float(len)*100.
        if(kount.gt.0) then
          if (me .eq. 0) then
          print *,'snow filled over glacier points at ',kount,
     &            ' points (',per,'percent)'
          endif
        endif
      endif ! landice check
      kount = 0
      do i=1,len
        if(slmask(i).eq.0.and.aisanl(i).eq.0) then
          snoanl(i) = 0.
          kount     = kount + 1
        endif
      enddo
      per = float(kount) / float(len)*100.
      if(kount.gt.0) then
        if (me .eq. 0) then
        print *,'snow set to zero over open sea at ',kount,
     &          ' points (',per,'percent)'
        endif
      endif
      return
      end subroutine qcsnow
      subroutine qcsice(ais,glacir,amxice,aicice,aicsea,sllnd,slmask,
     &                  rla,rlo,len,me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer kount1,kount,i,me,len
      real (kind=kind_io8) per,aicsea,aicice,sllnd
!
      real (kind=kind_io8) ais(len), glacir(len),
     &                     amxice(len), slmask(len)
      real (kind=kind_io8) rla(len), rlo(len)
!
!  check sea-ice cover mask against land-sea mask
!
      if (me .eq. 0) write(6,*) 'qc of sea ice'
      kount  = 0
      kount1 = 0
      do i=1,len
        if(ais(i).ne.aicice.and.ais(i).ne.aicsea) then
          print *,'sea ice mask not ',aicice,' or ',aicsea
          print *,'ais(i),aicice,aicsea,rla(i),rlo(i,=',
     &             ais(i),aicice,aicsea,rla(i),rlo(i)
          call abort
        endif
        if(slmask(i).eq.0..and.glacir(i).eq.1..and.
!       if(slmask(i).eq.0..and.glacir(i).eq.2..and.
     &     ais(i).ne.1.) then
          kount1 = kount1 + 1
          ais(i) = 1.
        endif
        if(slmask(i).eq.sllnd.and.ais(i).eq.aicice) then
          kount  = kount + 1
          ais(i) = aicsea
        endif
      enddo
!     enddo
      per = float(kount) / float(len)*100.
      if(kount.gt.0) then
        if(me .eq. 0) then
        print *,' sea ice over land mask at ',kount,' points (',per,
     &          'percent)'
        endif
      endif
      per = float(kount1) / float(len)*100.
      if(kount1.gt.0) then
        if(me .eq. 0) then
        print *,' sea ice set over glacier points over ocean at ',
     &          kount1,' points (',per,'percent)'
        endif
      endif
!     kount=0
!     do j=1,jdim
!     do i=1,idim
!       if(amxice(i,j).ne.0..and.ais(i,j).eq.0.) then
!         ais(i,j)=0.
!         kount=kount+1
!       endif
!     enddo
!     enddo
!     per=float(kount)/float(idim*jdim)*100.
!     if(kount.gt.0) then
!       print *,' sea ice exceeds maxice at ',kount,' points (',per,
!    &          'percent)'
!     endif
!
!  remove isolated open ocean surrounded by sea ice and/or land
!
!  remove isolated open ocean surrounded by sea ice and/or land
!
!     ij = 0
!     do j=1,jdim
!       do i=1,idim
!         ij = ij + 1
!         ip = i  + 1
!         im = i  - 1
!         jp = j  + 1
!         jm = j  - 1
!         if(jp.gt.jdim) jp = jdim - 1
!         if(jm.lt.1)    jm = 2
!         if(ip.gt.idim) ip = 1
!         if(im.lt.1)    im = idim
!         if(slmask(i,j).eq.0..and.ais(i,j).eq.0.) then
!           if((slmask(ip,jp).eq.1..or.ais(ip,jp).eq.1.).and.
!    &         (slmask(i ,jp).eq.1..or.ais(i ,jp).eq.1.).and.
!    &         (slmask(im,jp).eq.1..or.ais(im,jp).eq.1.).and.
!    &         (slmask(ip,j ).eq.1..or.ais(ip,j ).eq.1.).and.
!    &         (slmask(im,j ).eq.1..or.ais(im,j ).eq.1.).and.
!    &         (slmask(ip,jm).eq.1..or.ais(ip,jm).eq.1.).and.
!    &         (slmask(i ,jm).eq.1..or.ais(i ,jm).eq.1.).and.
!    &         (slmask(im,jm).eq.1..or.ais(im,jm).eq.1.)) then
!               ais(i,j) = 1.
!             write(6,*) ' isolated open sea point surrounded by',
!    &                   ' sea ice or land modified to sea ice',
!    &                   ' at lat=',rla(i,j),' lon=',rlo(i,j)
!           endif
!         endif
!       enddo
!     enddo
      return
      end
      subroutine setlsi(slmask,aisfld,len,aicice,slifld)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) aicice
      real (kind=kind_io8) slmask(len), slifld(len), aisfld(len)
!
!  set surface condition indicator slimsk
!
      do i=1,len
        slifld(i) = slmask(i)
!       if(aisfld(i).eq.aicice) slifld(i) = 2.0
        if(aisfld(i).eq.aicice .and. slmask(i) .eq. 0.0)
     &                                slifld(i) = 2.0
      enddo
      return
      end
      subroutine scale(fld,len,scl)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) fld(len),scl
      do i=1,len
        fld(i) = fld(i) * scl
      enddo
      return
      end
      subroutine qcmxmn(ttl,fld,slimsk,sno,iceflg,
     &                  fldlmx,fldlmn,fldomx,fldomn,fldimx,fldimn,
     &                  fldjmx,fldjmn,fldsmx,fldsmn,epsfld,
     &                  rla,rlo,len,mode,percrit,lgchek,me)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      real (kind=kind_io8) permax,per,fldimx,fldimn,fldjmx,fldomn,
     &                     fldlmx,fldlmn,fldomx,fldjmn,percrit,
     &                     fldsmx,fldsmn,epsfld
      integer kmaxi,kmini,kmaxj,kmino,kmaxl,kminl,kmaxo,mmprt,kminj,
     &        ij,nprt,kmaxs,kmins,i,me,len,mode
      parameter(mmprt=2)
!
      character*8 ttl
      logical iceflg(len)
      real (kind=kind_io8) fld(len),slimsk(len),sno(len),
     &                     rla(len), rlo(len)
      integer iwk(len)
      logical lgchek
!
      logical first
      integer   num_threads
      data first /.true./
      save num_threads, first
!
      integer len_thread_m, i1_t, i2_t, it
      integer num_parthds
!
      if (first) then
         num_threads = num_parthds()
         first = .false.
      endif
!
!  check against land-sea mask and ice cover mask
!
      if(me .eq. 0) then
!     print *,' '
      print *,'performing qc of ',ttl,' mode=',mode,
     &        '(0=count only, 1=replace)'
      endif
!
      len_thread_m  = (len+num_threads-1) / num_threads
      kmaxl = 0
      kminl = 0
      kmaxo = 0
      kmino = 0
      kmaxi = 0
      kmini = 0
      kmaxj = 0
      kminj = 0
      kmaxs = 0
      kmins = 0
!$omp parallel do private(i1_t,i2_t,it,i)
!$omp+private(nprt,ij,iwk)
!$omp+reduction(+:kmaxs,kmins,kmaxl,kminl,kmaxo)
!$omp+reduction(+:kmino,kmaxi,kmini,kmaxj,kminj)
!$omp+shared(mode,epsfld)
!$omp+shared(fldlmx,fldlmn,fldomx,fldjmn,fldsmx,fldsmn)
!$omp+shared(fld,slimsk,sno,rla,rlo)
      do it=1,num_threads   ! start of threaded loop
        i1_t = (it-1)*len_thread_m+1
        i2_t = min(i1_t+len_thread_m-1,len)
!
!  lower bound check over bare land
!
        if (fldlmn .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.1..and.sno(i).le.0..and.
     &         fld(i).lt.fldlmn-epsfld) then
               kminl=kminl+1
               iwk(kminl) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kminl)
            do i=1,nprt
              ij = iwk(i)
              print 8001,rla(ij),rlo(ij),fld(ij),fldlmn
 8001         format(' bare land min. check. lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e13.6, ' to ',e13.6)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kminl
              fld(iwk(i)) = fldlmn
            enddo
          endif
        endif
!
!  upper bound check over bare land
!
        if (fldlmx .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.1..and.sno(i).le.0..and.
     &         fld(i).gt.fldlmx+epsfld) then
               kmaxl=kmaxl+1
               iwk(kmaxl) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kmaxl)
            do i=1,nprt
              ij = iwk(i)
              print 8002,rla(ij),rlo(ij),fld(ij),fldlmx
 8002         format(' bare land max. check. lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e13.6, ' to ',e13.6)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kmaxl
              fld(iwk(i)) = fldlmx
            enddo
          endif
        endif
!
!  lower bound check over snow covered land
!
        if (fldsmn .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.1..and.sno(i).gt.0..and.
     &         fld(i).lt.fldsmn-epsfld) then
               kmins=kmins+1
               iwk(kmins) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kmins)
            do i=1,nprt
              ij = iwk(i)
              print 8003,rla(ij),rlo(ij),fld(ij),fldsmn
 8003         format(' sno covrd land min. check. lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e11.4, ' to ',e11.4)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kmins
              fld(iwk(i)) = fldsmn
            enddo
          endif
        endif
!
!  upper bound check over snow covered land
!
        if (fldsmx .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.1..and.sno(i).gt.0..and.
     &         fld(i).gt.fldsmx+epsfld) then
               kmaxs=kmaxs+1
               iwk(kmaxs) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kmaxs)
            do i=1,nprt
              ij = iwk(i)
              print 8004,rla(ij),rlo(ij),fld(ij),fldsmx
 8004         format(' snow land max. check. lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e11.4, ' to ',e11.4)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kmaxs
              fld(iwk(i)) = fldsmx
            enddo
          endif
        endif
!
!  lower bound check over open ocean
!
        if (fldomn .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.0..and.
     &         fld(i).lt.fldomn-epsfld) then
               kmino=kmino+1
               iwk(kmino) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kmino)
            do i=1,nprt
              ij = iwk(i)
              print 8005,rla(ij),rlo(ij),fld(ij),fldomn
 8005         format(' open ocean min. check. lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e11.4,' to ',e11.4)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kmino
              fld(iwk(i)) = fldomn
            enddo
          endif
      endif
!
!  upper bound check over open ocean
!
        if (fldomx .ne. 999.0) then
          do i=i1_t,i2_t
            if(fldomx.ne.999..and.slimsk(i).eq.0..and.
     &         fld(i).gt.fldomx+epsfld) then
               kmaxo=kmaxo+1
               iwk(kmaxo) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kmaxo)
            do i=1,nprt
              ij = iwk(i)
              print 8006,rla(ij),rlo(ij),fld(ij),fldomx
 8006         format(' open ocean max. check. lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e11.4, ' to ',e11.4)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kmaxo
              fld(iwk(i)) = fldomx
            enddo
          endif
        endif
!
!  lower bound check over sea ice without snow
!
        if (fldimn .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.2..and.sno(i).le.0..and.
     &         fld(i).lt.fldimn-epsfld) then
               kmini=kmini+1
               iwk(kmini) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kmini)
            do i=1,nprt
              ij = iwk(i)
              print 8007,rla(ij),rlo(ij),fld(ij),fldimn
 8007         format(' seaice no snow min. check lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e11.4, ' to ',e11.4)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kmini
              fld(iwk(i)) = fldimn
            enddo
          endif
        endif
!
!  upper bound check over sea ice without snow
!
        if (fldimx .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.2..and.sno(i).le.0..and.
     &         fld(i).gt.fldimx+epsfld  .and. iceflg(i)) then
!    &         fld(i).gt.fldimx+epsfld) then
               kmaxi=kmaxi+1
               iwk(kmaxi) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kmaxi)
            do i=1,nprt
              ij = iwk(i)
              print 8008,rla(ij),rlo(ij),fld(ij),fldimx
 8008         format(' seaice no snow max. check lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e11.4, ' to ',e11.4)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kmaxi
              fld(iwk(i)) = fldimx
            enddo
          endif
        endif
!
!  lower bound check over sea ice with snow
!
        if (fldjmn .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.2..and.sno(i).gt.0..and.
     &         fld(i).lt.fldjmn-epsfld) then
               kminj=kminj+1
               iwk(kminj) = i
            endif
          enddo
          if(me == 0 . and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kminj)
            do i=1,nprt
              ij = iwk(i)
              print 8009,rla(ij),rlo(ij),fld(ij),fldjmn
 8009         format(' sea ice snow min. check lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e11.4, ' to ',e11.4)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kminj
              fld(iwk(i)) = fldjmn
            enddo
          endif
        endif
!
!  upper bound check over sea ice with snow
!
        if (fldjmx .ne. 999.0) then
          do i=i1_t,i2_t
            if(slimsk(i).eq.2..and.sno(i).gt.0..and.
     &         fld(i).gt.fldjmx+epsfld  .and. iceflg(i)) then
!    &         fld(i).gt.fldjmx+epsfld) then
               kmaxj=kmaxj+1
               iwk(kmaxj) = i
            endif
          enddo
          if(me == 0 .and. it == 1 .and. num_threads == 1) then
            nprt = min(mmprt,kmaxj)
            do i=1,nprt
              ij = iwk(i)
              print 8010,rla(ij),rlo(ij),fld(ij),fldjmx
 8010         format(' seaice snow max check lat=',f5.1,
     &             ' lon=',f6.1,' fld=',e11.4, ' to ',e11.4)
            enddo
          endif
          if (mode .eq. 1) then
            do i=1,kmaxj
              fld(iwk(i)) = fldjmx
            enddo
          endif
        endif
      enddo            ! end of threaded loop
!$omp end parallel do
!
!  print results
!
      if(me .eq. 0) then
!     write(6,*) 'summary of qc'
      permax=0.
      if(kminl.gt.0) then
        per=float(kminl)/float(len)*100.
        print 9001,fldlmn,kminl,per
 9001   format(' bare land min check.  modified to ',f8.1,
     &         ' at ',i5,' points ',f8.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kmaxl.gt.0) then
        per=float(kmaxl)/float(len)*100.
        print 9002,fldlmx,kmaxl,per
 9002   format(' bare land max check. modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kmino.gt.0) then
        per=float(kmino)/float(len)*100.
        print 9003,fldomn,kmino,per
 9003   format(' open ocean min check.  modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kmaxo.gt.0) then
        per=float(kmaxo)/float(len)*100.
        print 9004,fldomx,kmaxo,per
 9004   format(' open sea max check. modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kmins.gt.0) then
        per=float(kmins)/float(len)*100.
        print 9009,fldsmn,kmins,per
 9009   format(' snow covered land min check. modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kmaxs.gt.0) then
        per=float(kmaxs)/float(len)*100.
        print 9010,fldsmx,kmaxs,per
 9010   format(' snow covered land max check. modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kmini.gt.0) then
        per=float(kmini)/float(len)*100.
        print 9005,fldimn,kmini,per
 9005   format(' bare ice min check.  modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kmaxi.gt.0) then
        per=float(kmaxi)/float(len)*100.
        print 9006,fldimx,kmaxi,per
 9006   format(' bare ice max check. modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kminj.gt.0) then
        per=float(kminj)/float(len)*100.
        print 9007,fldjmn,kminj,per
 9007   format(' snow covered ice min check.  modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
      if(kmaxj.gt.0) then
        per=float(kmaxj)/float(len)*100.
        print 9008,fldjmx,kmaxj,per
 9008   format(' snow covered ice max check. modified to ',f8.1,
     &         ' at ',i5,' points ',f4.1,'percent')
        if(per.gt.permax) permax=per
      endif
!     commented on 06/30/99  -- moorthi
!     if(lgchek) then
!       if(permax.gt.percrit) then
!         write(6,*) ' too many bad points.  aborting ....'
!         call abort
!       endif
!     endif
!
      endif
!
      return
      end
      subroutine setzro(fld,eps,len)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) fld(len),eps
      do i=1,len
        if(abs(fld(i)).lt.eps) fld(i) = 0.
      enddo
      return
      end
      subroutine getscv(snofld,scvfld,len)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) snofld(len),scvfld(len)
!
      do i=1,len
        scvfld(i) = 0.
        if(snofld(i).gt.0.) scvfld(i) = 1.
      enddo
      return
      end
      subroutine getstc(tsffld,tg3fld,slifld,len,lsoil,stcfld,tsfimx)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer k,i,len,lsoil
      real (kind=kind_io8) factor,tsfimx
      real (kind=kind_io8) tsffld(len), tg3fld(len), slifld(len)
      real (kind=kind_io8) stcfld(len,lsoil)
!
!  layer soil temperature
!
      do k = 1, lsoil
        do i = 1, len
          if(slifld(i).eq.1.0) then
            factor = ((k-1) * 2 + 1) / (2. * lsoil)
            stcfld(i,k) = factor*tg3fld(i)+(1.-factor)*tsffld(i)
          elseif(slifld(i).eq.2.0) then
            factor = ((k-1) * 2 + 1) / (2. * lsoil)
            stcfld(i,k) = factor*tsfimx+(1.-factor)*tsffld(i)
          else
            stcfld(i,k) = tg3fld(i)
          endif
        enddo
      enddo
      if(lsoil.gt.2) then
        do k = 3, lsoil
          do i = 1, len
            stcfld(i,k) = stcfld(i,2)
          enddo
        enddo
      endif
      return
      end
      subroutine getsmc(wetfld,len,lsoil,smcfld,me)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer k,i,len,lsoil,me
      real (kind=kind_io8) wetfld(len), smcfld(len,lsoil)
!
      if (me .eq. 0) write(6,*) 'getsmc'
!
!  layer soil wetness
!
      do k = 1, lsoil
        do i = 1, len
          smcfld(i,k) = (wetfld(i)*1000./150.)*.37 + .1
        enddo
      enddo
      return
      end
      subroutine usesgt(sig1t,slianl,tg3anl,len,lsoil,tsfanl,stcanl,
     &                  tsfimx)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len,lsoil
      real (kind=kind_io8) tsfimx
      real (kind=kind_io8) sig1t(len), slianl(len), tg3anl(len)
      real (kind=kind_io8) tsfanl(len), stcanl(len,lsoil)
!
!  soil temperature
!
      if(sig1t(1).gt.0.) then
        do i=1,len
          if(slianl(i).ne.0.) then
            tsfanl(i) = sig1t(i)
          endif
        enddo
      endif
      call getstc(tsfanl,tg3anl,slianl,len,lsoil,stcanl,tsfimx)
!
      return
      end
      subroutine snosfc(snoanl,tsfanl,tsfsmx,len,me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer kount,i,len,me
      real (kind=kind_io8) per,tsfsmx
      real (kind=kind_io8) snoanl(len), tsfanl(len)
!
      if (me .eq. 0) write(6,*) 'set snow temp to tsfsmx if greater'
      kount=0
      do i=1,len
        if(snoanl(i).gt.0.) then
          if(tsfanl(i).gt.tsfsmx) tsfanl(i)=tsfsmx
          kount = kount + 1
        endif
      enddo
      if(kount.gt.0) then
        if(me .eq. 0) then
        per=float(kount)/float(len)*100.
        write(6,*) 'snow sfc.  tsf set to ',tsfsmx,' at ',
     &              kount, ' points ',per,'percent'
        endif
      endif
      return
      end
      subroutine albocn(albclm,slmask,albomx,len)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) albomx
      real (kind=kind_io8) albclm(len,4), slmask(len)
      do i=1,len
        if(slmask(i).eq.0) then
          albclm(i,1) = albomx
          albclm(i,2) = albomx
          albclm(i,3) = albomx
          albclm(i,4) = albomx
        endif
      enddo
      return
      end
      subroutine qcmxice(glacir,amxice,len,me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,kount,len,me
      real (kind=kind_io8) glacir(len),amxice(len),per
      if (me .eq. 0) write(6,*) 'qc of maximum ice extent'
      kount=0
      do i=1,len
        if(glacir(i).eq.1..and.amxice(i).eq.0.) then
          amxice(i) = 0.
          kount     = kount + 1
        endif
      enddo
      if(kount.gt.0) then
        per = float(kount) / float(len)*100.
        if(me .eq. 0) write(6,*) ' max ice limit less than glacier'
     &,            ' coverage at ', kount, ' points ',per,'percent'
      endif
      return
      end
      subroutine qcsli(slianl,slifcs,len,me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,kount,len,me
      real (kind=kind_io8) slianl(len), slifcs(len),per
      if (me .eq. 0) then
      write(6,*) ' '
      write(6,*) 'qcsli'
      endif
      kount=0
      do i=1,len
        if(slianl(i).eq.1..and.slifcs(i).eq.0.) then
          kount      = kount + 1
          slifcs(i) = 1.
        endif
        if(slianl(i).eq.0..and.slifcs(i).eq.1.) then
          kount      = kount + 1
          slifcs(i) = 0.
        endif
        if(slianl(i).eq.2..and.slifcs(i).eq.1.) then
          kount      = kount + 1
          slifcs(i) = 0.
        endif
        if(slianl(i).eq.1..and.slifcs(i).eq.2.) then
          kount      = kount + 1
          slifcs(i) = 1.
        endif
      enddo
      if(kount.gt.0) then
        per=float(kount)/float(len)*100.
        if(me .eq. 0) then
        write(6,*) ' inconsistency of slmask between forecast and',
     &             ' analysis corrected at ',kount, ' points ',per,
     &             'percent'
        endif
      endif
      return
      end
!     subroutine nntprt(data,imax,fact)
!     real (kind=kind_io8) data(imax)
!     ilast=0
!     i1=1
!     i2=80
!1112 continue
!     if(i2.ge.imax) then
!       ilast=1
!       i2=imax
!     endif
!     write(6,*) ' '
!     do j=1,jmax
!       write(6,1111) (nint(data(imax*(j-1)+i)*fact),i=i1,i2)
!     enddo
!     if(ilast.eq.1) return
!     i1=i1+80
!     i2=i1+79
!     if(i2.ge.imax) then
!       ilast=1
!       i2=imax
!     endif
!     go to 1112
!1111 format(80i1)
!     return
!     end
      subroutine qcbyfc(tsffcs,snofcs,qctsfs,qcsnos,qctsfi,
     &                  len,lsoil,snoanl,aisanl,slianl,tsfanl,albanl,
     &                  zoranl,smcanl,
     &                  smcclm,tsfsmx,albomx,zoromx, me)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer kount,me,k,i,lsoil,len
      real (kind=kind_io8) zoromx,per,albomx,qctsfi,qcsnos,qctsfs,tsfsmx
      real (kind=kind_io8) tsffcs(len), snofcs(len)
      real (kind=kind_io8) snoanl(len), aisanl(len),
     &     slianl(len), zoranl(len),
     &     tsfanl(len), albanl(len,4),
     &     smcanl(len,lsoil)
      real (kind=kind_io8) smcclm(len,lsoil)
!
      if (me .eq. 0) write(6,*) 'qc of snow and sea-ice analysis'
!
! qc of snow analysis
!
!  questionable snow cover
!
      kount = 0
      do i=1,len
        if(slianl(i).gt.0..and.
     &     tsffcs(i).gt.qctsfs.and.snoanl(i).gt.0.) then
          kount      = kount + 1
          snoanl(i) = 0.
          tsfanl(i) = tsffcs(i)
        endif
      enddo
      if(kount.gt.0) then
        per=float(kount)/float(len)*100.
        if (me .eq. 0) then
        write(6,*) ' guess surface temp .gt. ',qctsfs,
     &             ' but snow analysis indicates snow cover'
        write(6,*) ' snow analysis set to zero',
     &             ' at ',kount, ' points ',per,'percent'
        endif
      endif
!
!  questionable no snow cover
!
      kount = 0
      do i=1,len
        if(slianl(i).gt.0..and.
     &     snofcs(i).gt.qcsnos.and.snoanl(i).lt.0.) then
          kount      = kount + 1
          snoanl(i) = snofcs(i)
          tsfanl(i) = tsffcs(i)
        endif
      enddo
      if(kount.gt.0) then
        per=float(kount)/float(len)*100.
        if (me .eq. 0) then
        write(6,*) ' guess snow depth .gt. ',qcsnos,
     &             ' but snow analysis indicates no snow cover'
        write(6,*) ' snow analysis set to guess value',
     &             ' at ',kount, ' points ',per,'percent'
        endif
      endif
!
!  questionable sea ice cover ! this qc is disable to correct error in
!  surface temparature over observed sea ice points
!
!     kount = 0
!     do i=1,len
!       if(slianl(i).eq.2..and.
!    &     tsffcs(i).gt.qctsfi.and.aisanl(i).eq.1.) then
!         kount        = kount + 1
!         aisanl(i)   = 0.
!         slianl(i)   = 0.
!         tsfanl(i)   = tsffcs(i)
!         snoanl(i)   = 0.
!         zoranl(i)   = zoromx
!         albanl(i,1) = albomx
!         albanl(i,2) = albomx
!         albanl(i,3) = albomx
!         albanl(i,4) = albomx
!         do k=1,lsoil
!           smcanl(i,k) = smcclm(i,k)
!         enddo
!       endif
!     enddo
!     if(kount.gt.0) then
!       per=float(kount)/float(len)*100.
!       if (me .eq. 0) then
!       write(6,*) ' guess surface temp .gt. ',qctsfi,
!    &             ' but sea-ice analysis indicates sea-ice'
!       write(6,*) ' sea-ice analysis set to zero',
!    &             ' at ',kount, ' points ',per,'percent'
!       endif
!     endif
!
      return
      end
      subroutine setrmsk(kpds5,slmask,igaul,jgaul,wlon,rnlat,
     &                   data,imax,jmax,rlnout,rltout,lmask,rslmsk
     &,                  gaus,blno, blto, kgds1, kpds4, lbms)
      use machine , only : kind_io8,kind_io4
      use sfccyc_module
      implicit none
      real (kind=kind_io8) blno,blto,wlon,rnlat,crit,data_max
      integer i,j,ijmax,jgaul,igaul,kpds5,jmax,imax, kgds1, kspla
      integer, intent(in)   :: kpds4
      logical*1, intent(in) :: lbms(imax,jmax)
      real*4                :: dummy(imax,jmax)

      real (kind=kind_io8)    slmask(igaul,jgaul)
      real (kind=kind_io8)    data(imax,jmax),rslmsk(imax,jmax)
     &,                       rlnout(imax), rltout(jmax)
      real (kind=kind_io8)    a(jmax), w(jmax), radi, dlat, dlon
      logical lmask, gaus
!
!     set the longitude and latitudes for the grib file
!
      if (kgds1 .eq. 4) then         ! grib file on gaussian grid
        kspla=4
        call splat(kspla, jmax, a, w)
!
        radi = 180.0 / (4.*atan(1.))
        do  j=1,jmax
          rltout(j) = acos(a(j)) * radi
        enddo
!
        if (rnlat .gt. 0.0) then
          do j=1,jmax
            rltout(j) = 90. - rltout(j)
          enddo
        else
          do j=1,jmax
            rltout(j) = -90. + rltout(j)
          enddo
        endif
      elseif (kgds1 .eq. 0) then     ! grib file on lat/lon grid
        dlat = -(rnlat+rnlat) / float(jmax-1)
        do j=1,jmax
         rltout(j) = rnlat + (j-1) * dlat
        enddo
      else                           ! grib file on some other grid
        call abort
      endif
      dlon = 360.0 / imax
      do i=1,imax
        rlnout(i) = wlon + (i-1)*dlon
      enddo
!
!
      ijmax  = imax*jmax
      rslmsk = 0.
!
!  surface temperature
!
      if(kpds5.eq.kpdtsf) then
!       lmask=.false.
        call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
     &,            rlnout, rltout, gaus, blno, blto)
!    &,            dlon, dlat, gaus, blno, blto)
        crit=0.5
        call rof01(rslmsk,ijmax,'ge',crit)
        lmask=.true.
!
!  bucket soil wetness
!
      elseif(kpds5.eq.kpdwet) then
        call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
     &,            rlnout, rltout, gaus, blno, blto)
!    &,            dlon, dlat, gaus, blno, blto)
        crit=0.5
        call rof01(rslmsk,ijmax,'ge',crit)
        lmask=.true.
!       write(6,*) 'wet rslmsk'
!       znnt=1.
!       call nntprt(rslmsk,ijmax,znnt)
!
!  snow depth
!
      elseif(kpds5.eq.kpdsnd) then
        if(kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.
              end if
            enddo
          enddo
          lmask=.true.
        else
          lmask=.false.
        end if
!
! snow liq equivalent depth
!
      elseif(kpds5.eq.kpdsno) then
        call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
     &,            rlnout, rltout, gaus, blno, blto)
!    &,            dlon, dlat, gaus, blno, blto)
        crit=0.5
        call rof01(rslmsk,ijmax,'ge',crit)
        lmask=.true.
!       write(6,*) 'sno rslmsk'
!       znnt=1.
!       call nntprt(rslmsk,ijmax,znnt)
!
!  soil moisture
!
      elseif(kpds5.eq.kpdsmc) then
        if(kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.
              end if
            enddo
          enddo
          lmask=.true.
        else
          call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
     &,            rlnout, rltout, gaus, blno, blto)
          crit=0.5
          call rof01(rslmsk,ijmax,'ge',crit)
          lmask=.true.
        endif
!
!  surface roughness
!
      elseif(kpds5.eq.kpdzor) then
        do j=1,jmax
          do i=1,imax
            rslmsk(i,j)=data(i,j)
          enddo
        enddo
        crit=9.9
        call rof01(rslmsk,ijmax,'lt',crit)
        lmask=.true.
!       write(6,*) 'zor rslmsk'
!       znnt=1.
!       call nntprt(rslmsk,ijmax,znnt)
!
!  albedo
!
!     elseif(kpds5.eq.kpdalb) then
!       do j=1,jmax
!         do i=1,imax
!           rslmsk(i,j)=data(i,j)
!         enddo
!       enddo
!       crit=99.
!       call rof01(rslmsk,ijmax,'lt',crit)
!       lmask=.true.
!       write(6,*) 'alb rslmsk'
!       znnt=1.
!       call nntprt(rslmsk,ijmax,znnt)
!
!  albedo
!
!cbosu  new snowfree albedo database has bitmap, use it.
      elseif(kpds5.eq.kpdalb(1)) then
        if (kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.  
              end if
            enddo
          enddo
          lmask = .true.
        else  ! no bitmap. old database has no water flag.
          lmask=.false.
        end if
      elseif(kpds5.eq.kpdalb(2)) then
!cbosu
        if (kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.  
              end if
            enddo
          enddo
          lmask = .true.
        else  ! no bitmap. old database has no water flag.
          lmask=.false.
        end if
      elseif(kpds5.eq.kpdalb(3)) then
!cbosu
        if (kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.  
              end if
            enddo
          enddo
          lmask = .true.
        else  ! no bitmap. old database has no water flag.
          lmask=.false.
        end if
      elseif(kpds5.eq.kpdalb(4)) then
!cbosu
        if (kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.  
              end if
            enddo
          enddo
          lmask = .true.
        else  ! no bitmap. old database has no water flag.
          lmask=.false.
        end if
!
!  vegetation fraction for albedo
!
      elseif(kpds5.eq.kpdalf(1)) then
!       rslmsk=data
!       crit=0.
!       call rof01(rslmsk,ijmax,'gt',crit)
!       lmask=.true.
        lmask=.false.
      elseif(kpds5.eq.kpdalf(2)) then
!       rslmsk=data
!       crit=0.
!       call rof01(rslmsk,ijmax,'gt',crit)
!       lmask=.true.
        lmask=.false.
!
!  sea ice
!
      elseif(kpds5.eq.kpdais) then
        lmask=.false.
!       call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
!    &,            dlon, dlat, gaus, blno, blto)
!       crit=0.5
!       call rof01(rslmsk,ijmax,'ge',crit)
!
        data_max = 0.0
        do j=1,jmax
          do i=1,imax
              rslmsk(i,j) = data(i,j)
              data_max= max(data_max,data(i,j))
          enddo
        enddo
        crit=1.0
        if (data_max .gt. crit) then
          call rof01(rslmsk,ijmax,'gt',crit)
          lmask=.true.
        else
          lmask=.false.
        endif
!       write(6,*) 'acn rslmsk'
!       znnt=1.
!       call nntprt(rslmsk,ijmax,znnt)
!
!  deep soil temperature
!
      elseif(kpds5.eq.kpdtg3) then
        lmask=.false.
!       call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
!    &,            rlnout, rltout, gaus, blno, blto)
!    &,            dlon, dlat, gaus, blno, blto)
!       crit=0.5
!       call rof01(rslmsk,ijmax,'ge',crit)
!       lmask=.true.
!
!  plant resistance
!
!     elseif(kpds5.eq.kpdplr) then
!       call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
!    &,            rlnout, rltout, gaus, blno, blto)
!    &,            dlon, dlat, gaus, blno, blto)
!       crit=0.5
!       call rof01(rslmsk,ijmax,'ge',crit)
!       lmask=.true.
!
!       write(6,*) 'plr rslmsk'
!       znnt=1.
!       call nntprt(rslmsk,ijmax,znnt)
!
!  glacier points
!
      elseif(kpds5.eq.kpdgla) then
        lmask=.false.
!
!  max ice extent
!
      elseif(kpds5.eq.kpdmxi) then
        lmask=.false.
!
!  snow cover
!
      elseif(kpds5.eq.kpdscv) then
        call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
     &,            rlnout, rltout, gaus, blno, blto)
!    &,            dlon, dlat, gaus, blno, blto)
        crit=0.5
        call rof01(rslmsk,ijmax,'ge',crit)
        lmask=.true.
!       write(6,*) 'scv rslmsk'
!       znnt=1.
!       call nntprt(rslmsk,ijmax,znnt)
!
!  sea ice concentration
!
      elseif(kpds5.eq.kpdacn) then
        lmask=.false.
        call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
     &,            rlnout, rltout, gaus, blno, blto)
!    &,            dlon, dlat, gaus, blno, blto)
        crit=0.5
        call rof01(rslmsk,ijmax,'ge',crit)
        lmask=.true.
!       write(6,*) 'acn rslmsk'
!       znnt=1.
!       call nntprt(rslmsk,ijmax,znnt)
!
!  vegetation cover
!
      elseif(kpds5.eq.kpdveg) then
!cggg
        if (kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,jmax-j+1) = 1.  ! need to flip grid in n/s direction
              end if
            enddo
          enddo
          lmask = .true.
        else  ! no bitmap, set mask the old way.

          call ga2la(slmask,igaul,jgaul,rslmsk,imax,jmax,wlon,rnlat
     &,              rlnout, rltout, gaus, blno, blto)
          crit=0.5
          call rof01(rslmsk,ijmax,'ge',crit)
          lmask=.true.

       end if
!
!  soil type
!
      elseif(kpds5.eq.kpdsot) then

        if (kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.
              end if
            enddo
          enddo
!  soil type is zero over water, use this to get a bitmap.
        else
          do j = 1, jmax
          do i = 1, imax
            rslmsk(i,j) = data(i,j)
          enddo
          enddo
          crit=0.1
          call rof01(rslmsk,ijmax,'gt',crit)
        endif
        lmask=.true.
!
!  vegetation type
!
      elseif(kpds5.eq.kpdvet) then

        if (kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.
              end if
            enddo
          enddo
!  veg type is zero over water, use this to get a bitmap.
        else
          do j = 1, jmax
          do i = 1, imax
            rslmsk(i,j) = data(i,j)
          enddo
          enddo
          crit=0.1
          call rof01(rslmsk,ijmax,'gt',crit)
        endif
        lmask=.true.
!
!      these are for four new data type added by clu -- not sure its correct!
!
      elseif(kpds5.eq.kpdvmn) then
!
!cggg  greenness is zero over water, use this to get a bitmap.
!
        do j = 1, jmax
          do i = 1, imax
            rslmsk(i,j) = data(i,j)
          enddo
        enddo
!
        crit=0.1
        call rof01(rslmsk,ijmax,'gt',crit)
        lmask=.true.
!cggg        lmask=.false.
!
      elseif(kpds5.eq.kpdvmx) then
!
!cggg  greenness is zero over water, use this to get a bitmap.
!
        do j = 1, jmax
          do i = 1, imax
            rslmsk(i,j) = data(i,j)
          enddo
        enddo
!
        crit=0.1
        call rof01(rslmsk,ijmax,'gt',crit)
        lmask=.true.
!cggg        lmask=.false.
!
      elseif(kpds5.eq.kpdslp) then
!
!cggg slope type is zero over water, use this to get a bitmap.
!
        do j = 1, jmax
          do i = 1, imax
            rslmsk(i,j) = data(i,j)
          enddo
        enddo
!
        crit=0.1
        call rof01(rslmsk,ijmax,'gt',crit)
        lmask=.true.
!cggg        lmask=.false.
!
!cbosu new maximum snow albedo database has bitmap
      elseif(kpds5.eq.kpdabs) then
        if (kpds4 == 192) then  ! use the bitmap
          rslmsk = 0.
          do j = 1, jmax
            do i = 1, imax
              if (lbms(i,j)) then
                rslmsk(i,j) = 1.  
              end if
            enddo
          enddo
          lmask = .true.
        else  ! no bitmap. old database has zero over water
          do j = 1, jmax
            do i = 1, imax
              rslmsk(i,j) = data(i,j)
            enddo
          enddo
          crit=0.1
          call rof01(rslmsk,ijmax,'gt',crit)
          lmask=.true.
        end if
      endif
!
      return
      end
      subroutine ga2la(gauin,imxin,jmxin,regout,imxout,jmxout,
     &                 wlon,rnlat,rlnout,rltout,gaus,blno, blto)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i1,i2,j2,ishft,i,jj,j1,jtem,jmxout,imxin,jmxin,imxout,
     &        j,iret
      real (kind=kind_io8) alamd,dxin,aphi,x,sum1,sum2,y,dlati,wlon,
     &                     rnlat,dxout,dphi,dlat,facns,tem,blno,
     &                     blto
!
!  interpolation from lat/lon grid to other lat/lon grid
!
      real (kind=kind_io8) gauin (imxin,jmxin), regout(imxout,jmxout)
     &,                    rlnout(imxout), rltout(jmxout)
      logical gaus
!
      real, allocatable :: gaul(:)
      real (kind=kind_io8) ddx(imxout),ddy(jmxout)
      integer iindx1(imxout), iindx2(imxout),
     &        jindx1(jmxout), jindx2(jmxout)
      integer jmxsav,n,kspla
      data    jmxsav/0/
      save    jmxsav, gaul, dlati
      real (kind=kind_io8) radi
      real (kind=kind_io8) a(jmxin), w(jmxin)
!
!
      logical first
      integer   num_threads
      data first /.true./
      save num_threads, first
!
      integer len_thread_m, j1_t, j2_t, it
      integer num_parthds
!
      if (first) then
         num_threads = num_parthds()
         first = .false.
      endif
!
      if (jmxin .ne. jmxsav) then
        if (jmxsav .gt. 0) deallocate (gaul, stat=iret)
        allocate (gaul(jmxin))
        jmxsav = jmxin
        if (gaus) then
cjfe      call gaulat(gaul,jmxin)
cjfe
!
          kspla=4
          call splat(kspla, jmxin, a, w)
!
          radi = 180.0 / (4.*atan(1.))
          do  n=1,jmxin
            gaul(n) = acos(a(n)) * radi
          enddo
cjfe
          do j=1,jmxin
            gaul(j) = 90. - gaul(j)
          enddo
        else
          dlat = -2*blto / float(jmxin-1)
          dlati = 1 / dlat
          do j=1,jmxin
           gaul(j) = blto + (j-1) * dlat
          enddo
        endif
      endif
!
!
      dxin  = 360. / float(imxin )
!
      do i=1,imxout
        alamd = rlnout(i)
        i1     = floor((alamd-blno)/dxin) + 1
        ddx(i) = (alamd-blno)/dxin-(i1-1)
        iindx1(i) = modulo(i1-1,imxin) + 1
        iindx2(i) = modulo(i1  ,imxin) + 1
      enddo
!
!
      len_thread_m  = (jmxout+num_threads-1) / num_threads
!
      if (gaus) then
!
!$omp parallel do private(j1_t,j2_t,it,j1,j2,jj)
!$omp+private(aphi)
!$omp+shared(num_threads,len_thread_m)
!$omp+shared(jmxin,jmxout,gaul,rltout,jindx1,ddy)
!
        do it=1,num_threads   ! start of threaded loop ...................
          j1_t       = (it-1)*len_thread_m+1
          j2_t       = min(j1_t+len_thread_m-1,jmxout)
!
          j2=1
          do 40 j=j1_t,j2_t
            aphi=rltout(j)
            do 50 jj=1,jmxin
              if(aphi.lt.gaul(jj)) go to 50
              j2=jj
              go to 42
   50       continue
   42       continue
            if(j2.gt.2) go to 43
            j1=1
            j2=2
            go to 44
   43       continue
            if(j2.le.jmxin) go to 45
            j1=jmxin-1
            j2=jmxin
            go to 44
   45       continue
            j1=j2-1
   44       continue
            jindx1(j)=j1
            jindx2(j)=j2
            ddy(j)=(aphi-gaul(j1))/(gaul(j2)-gaul(j1))
   40     continue
        enddo             ! end of threaded loop ...................
!$omp   end parallel do
!
      else
!$omp parallel do private(j1_t,j2_t,it,j1,j2,jtem)
!$omp+private(aphi)
!$omp+shared(num_threads,len_thread_m)
!$omp+shared(jmxin,jmxout,gaul,rltout,jindx1,ddy,dlati,blto)
!
        do it=1,num_threads   ! start of threaded loop ...................
          j1_t       = (it-1)*len_thread_m+1
          j2_t       = min(j1_t+len_thread_m-1,jmxout)
!
          j2=1
          do 400 j=j1_t,j2_t
            aphi=rltout(j)
            jtem = (aphi - blto) * dlati + 1
            if (jtem .ge. 1 .and. jtem .lt. jmxin) then
              j1 = jtem
              j2 = j1 + 1
              ddy(j)=(aphi-gaul(j1))/(gaul(j2)-gaul(j1))
            elseif (jtem .eq. jmxin) then
              j1 = jmxin
              j2 = jmxin
              ddy(j)=1.0
            else
              j1 = 1
              j2 = 1
              ddy(j)=1.0
            endif
!
            jindx1(j) = j1
            jindx2(j) = j2
  400     continue
        enddo             ! end of threaded loop ...................
!$omp   end parallel do
      endif
!
!     write(6,*) 'ga2la'
!     write(6,*) 'iindx1'
!     write(6,*) (iindx1(n),n=1,imxout)
!     write(6,*) 'iindx2'
!     write(6,*) (iindx2(n),n=1,imxout)
!     write(6,*) 'jindx1'
!     write(6,*) (jindx1(n),n=1,jmxout)
!     write(6,*) 'jindx2'
!     write(6,*) (jindx2(n),n=1,jmxout)
!     write(6,*) 'ddy'
!     write(6,*) (ddy(n),n=1,jmxout)
!     write(6,*) 'ddx'
!     write(6,*) (ddx(n),n=1,jmxout)
!
!
!$omp parallel do private(j1_t,j2_t,it,i,i1,i2)
!$omp+private(j,j1,j2,x,y)
!$omp+shared(num_threads,len_thread_m)
!$omp+shared(imxout,iindx1,jindx1,ddx,ddy,gauin,regout)
!
      do it=1,num_threads   ! start of threaded loop ...................
        j1_t       = (it-1)*len_thread_m+1
        j2_t       = min(j1_t+len_thread_m-1,jmxout)
!
        do  j=j1_t,j2_t
          y  = ddy(j)
          j1 = jindx1(j)
          j2 = jindx2(j)
          do i=1,imxout
            x  = ddx(i)
            i1 = iindx1(i)
            i2 = iindx2(i)
            regout(i,j) = (1.-x)*((1.-y)*gauin(i1,j1) + y*gauin(i1,j2))
     &                  +     x *((1.-y)*gauin(i2,j1) + y*gauin(i2,j2))
          enddo
        enddo
      enddo             ! end of threaded loop ...................
!$omp end parallel do
!
      sum1 = 0.
      sum2 = 0.
      do i=1,imxin
        sum1 = sum1 + gauin(i,1)
        sum2 = sum2 + gauin(i,jmxin)
      enddo
      sum1 = sum1 / float(imxin)
      sum2 = sum2 / float(imxin)
!
      if (gaus) then
        if (rnlat .gt. 0.0) then
          do i=1,imxout
            regout(i,     1) = sum1
            regout(i,jmxout) = sum2
          enddo
        else
          do i=1,imxout
            regout(i,     1) = sum2
            regout(i,jmxout) = sum1
          enddo
        endif
      else
        if (blto .lt. 0.0) then
          if (rnlat .gt. 0.0) then
            do i=1,imxout
              regout(i,     1) = sum2
              regout(i,jmxout) = sum1
            enddo
          else
            do i=1,imxout
              regout(i,     1) = sum1
              regout(i,jmxout) = sum2
            enddo
          endif
        else
          if (rnlat .lt. 0.0) then
            do i=1,imxout
              regout(i,     1) = sum2
              regout(i,jmxout) = sum1
            enddo
          else
            do i=1,imxout
              regout(i,     1) = sum1
              regout(i,jmxout) = sum2
            enddo
          endif
        endif
      endif
!
      return
      end
      subroutine landtyp(vegtype,soiltype,slptype,slmask,len)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) vegtype(len),soiltype(len),slmask(len)
     +,                    slptype(len)  
!
!  make sure that the soil type and veg type are non-zero over land
!
      do i = 1, len
        if (slmask(i) .eq. 1) then
          if (vegtype(i)  .eq. 0.)  vegtype(i)  = 7
          if (soiltype(i) .eq. 0.)  soiltype(i) = 2
          if (slptype(i)  .eq. 0.)  slptype(i)  = 1
        endif
      enddo
      return
      end subroutine landtyp
      subroutine gaulat(gaul,k)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer n,k
      real (kind=kind_io8) radi
      real (kind=kind_io8) a(k), w(k), gaul(k)
!
      call splat(4, k, a, w)
!
      radi = 180.0 / (4.*atan(1.))
      do  n=1,k
        gaul(n) = acos(a(n)) * radi
      enddo
!
!     print *,'gaussian lat (deg) for jmax=',k
!     print *,(gaul(n),n=1,k)
!
      return
   70 write(6,6000)
 6000 format(//5x,'error in gauaw'//)
      stop
      end
!-----------------------------------------------------------------------
      subroutine anomint(tsfan0,tsfclm,tsfcl0,tsfanl,len)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,len
      real (kind=kind_io8) tsfanl(len), tsfan0(len),
     &                     tsfclm(len), tsfcl0(len)
!
!  time interpolation of anomalies
!  add initial anomaly to date interpolated climatology
!
      write(6,*) 'anomint'
      do i=1,len
        tsfanl(i) = tsfan0(i) - tsfcl0(i) + tsfclm(i)
      enddo
      return
      end
      subroutine clima(lugb,iy,im,id,ih,fh,len,lsoil,
     &                 slmask,fntsfc,fnwetc,fnsnoc,fnzorc,fnalbc,fnaisc,
     &                 fntg3c,fnscvc,fnsmcc,fnstcc,fnacnc,fnvegc,
     &                 fnvetc,fnsotc,
     &                 fnvmnc,fnvmxc,fnslpc,fnabsc,
     &                 tsfclm,tsfcl2,wetclm,snoclm,zorclm,albclm,aisclm,
     &                 tg3clm,cvclm ,cvbclm,cvtclm,
     &                 cnpclm,smcclm,stcclm,sliclm,scvclm,acnclm,vegclm,
     &                 vetclm,sotclm,alfclm,
     &                 vmnclm,vmxclm,slpclm,absclm,
     &                 kpdtsf,kpdwet,kpdsno,kpdzor,kpdalb,kpdais,
     &                 kpdtg3,kpdscv,kpdacn,kpdsmc,kpdstc,kpdveg,
     &                 kpdvet,kpdsot,kpdalf,tsfcl0,
     &                 kpdvmn,kpdvmx,kpdslp,kpdabs,
     &                 deltsfc, lanom
     &,                imsk, jmsk, slmskh, outlat, outlon
     &,                gaus, blno, blto, me,lprnt,iprnt, fnalbc2, ialb)
!
      use machine , only : kind_io8,kind_io4
      implicit none
      real (kind=kind_io8) rjday,wei1x,wei2x,rjdayh,wei2m,wei1m,wei1s,
     &                     wei2s,fh,stcmon1s,blto,blno,deltsfc,rjdayh2
      real (kind=kind_io8) wei1y,wei2y
      integer jdoy,jday,jh,jdow,mmm,mmp,mm,iret,monend,i,k,jm,jd,iy4,
     &        jy,mon1,is2,isx,kpd9,is1,l,nn,mon2,mon,is,kpdsno,
     &        kpdzor,kpdtsf,kpdwet,kpdscv,kpdacn,kpdais,kpdtg3,im,id,
     &        lugb,iy,len,lsoil,ih,kpdsmc,iprnt,me,m1,m2,k1,k2,
     &        kpdvet,kpdsot,kpdstc,kpdveg,jmsk,imsk,j,ialb
     &,       kpdvmn,kpdvmx,kpdslp,kpdabs,landice_cat
      integer kpdalb(4), kpdalf(2)
!
      character*500 fntsfc,fnwetc,fnsnoc,fnzorc,fnalbc,fnaisc,
     &             fntg3c,fnscvc,fnsmcc,fnstcc,fnacnc,fnvegc,
     &             fnvetc,fnsotc,fnalbc2
     &,            fnvmnc,fnvmxc,fnslpc,fnabsc
      real (kind=kind_io8) tsfclm(len),tsfcl2(len),
     &     wetclm(len),snoclm(len),
     &     zorclm(len),albclm(len,4),aisclm(len),
     &     tg3clm(len),acnclm(len),
     &     cvclm (len),cvbclm(len),cvtclm(len),
     &     cnpclm(len),
     &     smcclm(len,lsoil),stcclm(len,lsoil),
     &     sliclm(len),scvclm(len),vegclm(len),
     &     vetclm(len),sotclm(len),alfclm(len,2)
     &,    vmnclm(len),vmxclm(len),slpclm(len),absclm(len)
      real (kind=kind_io8) slmskh(imsk,jmsk)
      real (kind=kind_io8) outlat(len), outlon(len)
!
      real (kind=kind_io8) slmask(len), tsfcl0(len)
      real (kind=kind_io8), allocatable :: slmask_noice(:)
!
      logical lanom, gaus, first
!
! set z0 based on sib vegetation type
      real (kind=kind_io8) z0_sib(13)
      data z0_sib /2.653, 0.826, 0.563, 1.089, 0.854, 0.856,
     &             0.035, 0.238, 0.065, 0.076, 0.011, 0.125,
     &             0.011 /
! set z0 based on igbp vegetation type
      real (kind=kind_io8) z0_igbp_min(20), z0_igbp_max(20)
      real (kind=kind_io8) z0_season(12)
      data z0_igbp_min /1.089, 2.653, 0.854, 0.826, 0.800, 0.050,
     &                  0.030, 0.856, 0.856, 0.150, 0.040, 0.130,
     &                  1.000, 0.250, 0.011, 0.011, 0.001, 0.076,
     &                  0.050, 0.030/
      data z0_igbp_max /1.089, 2.653, 0.854, 0.826, 0.800, 0.050,
     &                  0.030, 0.856, 0.856, 0.150, 0.040, 0.130,
     &                  1.000, 0.250, 0.011, 0.011, 0.001, 0.076,
     &                  0.050, 0.030/
!
! dayhf : julian day of the middle of each month
!
      real (kind=kind_io8), parameter :: dayhf(13) = 
     &                (/ 15.5, 45.0, 74.5,105.0,135.5,166.0,
     &          196.5,227.5,258.0,288.5,319.0,349.5,380.5 /)
!
      real (kind=kind_io8) fha(5)
      real(4) fha4(5)
      integer w3kindreal,w3kindint
      integer ida(8),jda(8),ivtyp, kpd7
!
      real (kind=kind_io8), allocatable :: tsf(:,:),sno(:,:),
     &                     zor(:,:),wet(:,:),
     &                     ais(:,:), acn(:,:),   scv(:,:), smc(:,:,:),
     &                     tg3(:),   alb(:,:,:), alf(:,:),
     &                     vet(:),   sot(:),     tsf2(:),
     &                     veg(:,:), stc(:,:,:)
     &,                    vmn(:), vmx(:),  slp(:), abs(:)
!
      integer mon1s, mon2s, sea1s, sea2s, sea1, sea2, hyr1, hyr2
      data first/.true./
      data mon1s/0/, mon2s/0/, sea1s/0/, sea2s/0/
!
      save first, tsf, sno, zor, wet,  ais, acn, scv, smc, tg3,
     &     alb,   alf, vet, sot, tsf2, veg, stc
     &,    vmn,   vmx, slp, abs,
     &     mon1s, mon2s, sea1s, sea2s, k1, k2, m1, m2
!
      logical lprnt
!
      do i=1,len
        tsfclm(i) = 0.0
        tsfcl2(i) = 0.0
        snoclm(i) = 0.0
        wetclm(i) = 0.0
        zorclm(i) = 0.0
        aisclm(i) = 0.0
        tg3clm(i) = 0.0
        acnclm(i) = 0.0
        cvclm(i)  = 0.0
        cvbclm(i) = 0.0
        cvtclm(i) = 0.0
        cnpclm(i) = 0.0
        sliclm(i) = 0.0
        scvclm(i) = 0.0
        vmnclm(i) = 0.0
        vmxclm(i) = 0.0
        slpclm(i) = 0.0
        absclm(i) = 0.0
      enddo
      do k=1,lsoil
        do i=1,len
          smcclm(i,k) = 0.0
          stcclm(i,k) = 0.0
        enddo
      enddo
      do k=1,4
        do i=1,len
          albclm(i,k) = 0.0
        enddo
      enddo
      do k=1,2
        do i=1,len
          alfclm(i,k) = 0.0
        enddo
      enddo
!
      iret   = 0
      monend = 9999
!
      if (first) then
!
!    allocate variables to be saved
!
       allocate (tsf(len,2), sno(len,2),      zor(len,2),
     &           wet(len,2), ais(len,2),      acn(len,2),
     &           scv(len,2), smc(len,lsoil,2),
     &           tg3(len),   alb(len,4,2),    alf(len,2),
     &           vet(len),   sot(len), tsf2(len),
!clu [+1l] add vmn, vmx, slp, abs
     &           vmn(len),   vmx(len), slp(len), abs(len),
     &           veg(len,2), stc(len,lsoil,2))
!
!     get tsf climatology for the begining of the forecast
!
        if (fh .gt. 0.0) then
!cbosu
          if (me == 0) print*,'bosu fh gt 0'

          iy4=iy
          if(iy.lt.101) iy4=1900+iy4
          fha=0
          ida=0
          jda=0
!         fha(2)=nint(fh)
          ida(1)=iy
          ida(2)=im
          ida(3)=id
          ida(5)=ih
          call w3kind(w3kindreal,w3kindint)
          if(w3kindreal == 4) then
            fha4=fha
            call w3movdat(fha4,ida,jda)
          else
            call w3movdat(fha,ida,jda)
          endif
          jy=jda(1)
          jm=jda(2)
          jd=jda(3)
          jh=jda(5)
          if (me .eq. 0) write(6,*) ' forecast jy,jm,jd,jh',
     &                   jy,jm,jd,jh
          jdow = 0
          jdoy = 0
          jday = 0
          call w3doxdat(jda,jdow,jdoy,jday)
          rjday=jdoy+jda(5)/24.
          if(rjday.lt.dayhf(1)) rjday=rjday+365.
!
          if (me .eq. 0) write(6,*) 'forecast jy,jm,jd,jh=',jy,jm,jd,jh
!
!         for monthly mean climatology
!
          monend = 12
          do mm=1,monend
            mmm=mm
            mmp=mm+1
            if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
               mon1=mmm
               mon2=mmp
               go to 10
            endif
          enddo
          print *,'wrong rjday',rjday
          call abort
   10     continue
          wei1m = (dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
          wei2m = (rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
          if(mon2.eq.13) mon2=1
          if (me .eq. 0) print *,'rjday,mon1,mon2,wei1m,wei2m=',
     &                   rjday,mon1,mon2,wei1m,wei2m
!
!       read monthly mean climatology of tsf
!
          kpd7 = -1
          do nn=1,2
            mon = mon1
            if (nn .eq. 2) mon = mon2
            call fixrdc(lugb,fntsfc,kpdtsf,kpd7,mon,slmask,
     &                 tsf(1,nn),len,iret
     &,                imsk, jmsk, slmskh, gaus,blno, blto
     &,                outlat, outlon, me)
          enddo
!
!  tsf at the begining of forecast i.e. fh=0
!
          do i=1,len
            tsfcl0(i) = wei1m * tsf(i,1) + wei2m * tsf(i,2)
          enddo
        endif
      endif
!
!  compute current jy,jm,jd,jh of forecast and the day of the year
!
      iy4=iy
      if(iy.lt.101) iy4=1900+iy4
      fha    = 0
      ida    = 0
      jda    = 0
      fha(2) = nint(fh)
      ida(1) = iy
      ida(2) = im
      ida(3) = id
      ida(5) = ih
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        fha4=fha
        call w3movdat(fha4,ida,jda)
      else
        call w3movdat(fha,ida,jda)
      endif
      jy     = jda(1)
      jm     = jda(2)
      jd     = jda(3)
      jh     = jda(5)
!     if (me .eq. 0) write(6,*) ' forecast jy,jm,jd,jh,rjday=',
!    &               jy,jm,jd,jh,rjday
      jdow   = 0
      jdoy   = 0
      jday   = 0
      call w3doxdat(jda,jdow,jdoy,jday)
      rjday  = jdoy+jda(5)/24.
      if(rjday.lt.dayhf(1)) rjday=rjday+365.

      if (me .eq. 0) write(6,*) ' forecast jy,jm,jd,jh,rjday=',
     &               jy,jm,jd,jh,rjday
!
      if (me .eq. 0) write(6,*) 'forecast jy,jm,jd,jh=',jy,jm,jd,jh
!
!     for monthly mean climatology
!
      monend = 12
      do mm=1,monend
         mmm=mm
         mmp=mm+1
         if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
            mon1=mmm
            mon2=mmp
            go to 20
         endif
      enddo
      print *,'wrong rjday',rjday
      call abort
   20 continue
      wei1m=(dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
      wei2m=(rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
      if(mon2.eq.13) mon2=1
      if (me .eq. 0) print *,'rjday,mon1,mon2,wei1m,wei2m=',
     &               rjday,mon1,mon2,wei1m,wei2m
!
!     for seasonal mean climatology
!
      monend = 4
      is     = im/3 + 1
      if (is.eq.5) is = 1
      do mm=1,monend
        mmm = mm*3 - 2
        mmp = (mm+1)*3 - 2
        if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
          sea1 = mmm
          sea2 = mmp
          go to 30
        endif
      enddo
      print *,'wrong rjday',rjday
      call abort
   30 continue
      wei1s = (dayhf(sea2)-rjday)/(dayhf(sea2)-dayhf(sea1))
      wei2s = (rjday-dayhf(sea1))/(dayhf(sea2)-dayhf(sea1))
      if(sea2.eq.13) sea2=1
      if (me .eq. 0) print *,'rjday,sea1,sea2,wei1s,wei2s=',
     &               rjday,sea1,sea2,wei1s,wei2s
!
!     for summer and winter values (maximum and minimum).
!
      monend = 2
      is     = im/6 + 1
      if (is.eq.3) is = 1
      do mm=1,monend
        mmm = mm*6 - 5
        mmp = (mm+1)*6 - 5
        if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
          hyr1 = mmm
          hyr2 = mmp
          go to 31
        endif
      enddo
      print *,'wrong rjday',rjday
      call abort
   31 continue
      wei1y = (dayhf(hyr2)-rjday)/(dayhf(hyr2)-dayhf(hyr1))
      wei2y = (rjday-dayhf(hyr1))/(dayhf(hyr2)-dayhf(hyr1))
      if(hyr2.eq.13) hyr2=1
      if (me .eq. 0) print *,'rjday,hyr1,hyr2,wei1y,wei2y=',
     &               rjday,hyr1,hyr2,wei1y,wei2y
!
!  start reading in climatology and interpolate to the date
!
      first_time : if (first) then
!cbosu
        if (me == 0) print*,'bosu first time thru'
!
!     annual mean climatology
!
!  fraction of vegetation field for albedo --  there are two
!  fraction fields in this version: strong zenith angle dependent
!  and weak zenith angle dependent
!
        kpd9 = -1
cjfe
        alf=0.
cjfe

        kpd7=-1
        if (ialb == 1) then
!cbosu    still need facsf and facwf.  read them from the production
!cbosu    file
!cbosu
          call fixrdc(lugb,fnalbc2,kpdalf(1),kpd7,kpd9,slmask
     &,               alf,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
        else
          call fixrdc(lugb,fnalbc,kpdalf(1),kpd7,kpd9,slmask
     &,               alf,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
        endif
        do i = 1, len
          if(slmask(i).eq.1.) then
            alf(i,2) = 100. - alf(i,1)
          endif
        enddo
!
!  deep soil temperature
!
        if(fntg3c(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fntg3c,kpdtg3,kpd7,kpd9,slmask,
     &                tg3,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
        endif
!
!  vegetation type
!
!  when using the new gldas soil moisture climatology, a veg type
!  dataset must be selected.
!
        if(fnvetc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnvetc,kpdvet,kpd7,kpd9,slmask,
     &                vet,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological vegetation',
     &                              ' type read in.'
          landice_cat=13
          if (maxval(vet)> 13.0) landice_cat=15
        elseif(index(fnsmcc,'soilmgldas') /= 0) then ! new soil moisture climo
          if (me .eq. 0) write(6,*) 'fatal error: must choose'
          if (me .eq. 0) write(6,*) 'climatological veg type when'
          if (me .eq. 0) write(6,*) 'using new gldas soil moisture.'
          call abort
        endif
!
!  soil type
!
        if(fnsotc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnsotc,kpdsot,kpd7,kpd9,slmask,
     &                sot,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological soil type read in.'
        endif

!
!  min vegetation cover
!
        if(fnvmnc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnvmnc,kpdvmn,kpd7,kpd9,slmask,
     &                vmn,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological shdmin read in.'
        endif
!
!  max vegetation cover
!
        if(fnvmxc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnvmxc,kpdvmx,kpd7,kpd9,slmask,
     &                vmx,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological shdmax read in.'
        endif
!
!  slope type
!
        if(fnslpc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnslpc,kpdslp,kpd7,kpd9,slmask,
     &                slp,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological slope read in.'
        endif
!
!  max snow albeod
!
        if(fnabsc(1:8).ne.'        ') then
          kpd7=-1
          call fixrdc(lugb,fnabsc,kpdabs,kpd7,kpd9,slmask,
     &                abs,len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
          if (me .eq. 0) write(6,*) 'climatological snoalb read in.'
        endif
!clu ----------------------------------------------------------------------
!
        is1 = sea1/3 + 1
        is2 = sea2/3 + 1
        if (is1 .eq. 5) is1 = 1
        if (is2 .eq. 5) is2 = 1
        do nn=1,2
!
!     seasonal mean climatology
          if(nn.eq.1) then
             isx=is1
          else
             isx=is2
          endif
          if(isx.eq.1) kpd9 = 12
          if(isx.eq.2) kpd9 = 3
          if(isx.eq.3) kpd9 = 6
          if(isx.eq.4) kpd9 = 9
!
!         seasonal mean climatology
!
!  albedo
!  there are four albedo fields in this version:
!  two for strong zeneith angle dependent (visible and near ir)
!  and two for weak zeneith angle dependent (vis ans nir)
!
          if (ialb == 0) then
            kpd7=-1
            do k = 1, 4
              call fixrdc(lugb,fnalbc,kpdalb(k),kpd7,kpd9,slmask,
     &                    alb(1,k,nn),len,iret
     &,                   imsk, jmsk, slmskh, gaus,blno, blto
     &,                   outlat, outlon, me)
            enddo
          endif
!
!         monthly mean climatology
!
          mon = mon1
          if (nn .eq. 2) mon = mon2
!cbosu
!cbosu  new snowfree albedo database is monthly.  
          if (ialb == 1) then
            kpd7=-1
            do k = 1, 4
              call fixrdc(lugb,fnalbc,kpdalb(k),kpd7,mon,slmask,
     &                    alb(1,k,nn),len,iret
     &,                   imsk, jmsk, slmskh, gaus,blno, blto
     &,                   outlat, outlon, me)
            enddo
          endif

!     if(lprnt) print *,' mon1=',mon1,' mon2=',mon2
!
!  tsf at the current time t
!
          kpd7=-1
          call fixrdc(lugb,fntsfc,kpdtsf,kpd7,mon,slmask,
     &               tsf(1,nn),len,iret
     &,              imsk, jmsk, slmskh, gaus,blno, blto
     &,              outlat, outlon, me)
!     if(lprnt) print *,' tsf=',tsf(iprnt,nn),' nn=',nn
!
!  tsf...at time t-deltsfc
!
!     fh2 = fh - deltsfc
!     if (fh2 .gt. 0.0) then
!       call fixrd(lugb,fntsfc,kpdtsf,lclim,slmask,
!    &             iy,im,id,ih,fh2,tsfcl2,len,iret
!    &,            imsk, jmsk, slmskh, gaus,blno, blto
!    &,            outlat, outlon, me)
!     else
!       do i=1,len
!         tsfcl2(i) = tsfclm(i)
!       enddo
!     endif
!
!  soil wetness
!
          if(fnwetc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnwetc,kpdwet,kpd7,mon,slmask,
     &                  wet(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
          elseif(fnsmcc(1:8).ne.'        ') then
            if (index(fnsmcc,'global_soilmcpc.1x1.grb') /= 0) then ! the old climo data
              kpd7=-1
              call fixrdc(lugb,fnsmcc,kpdsmc,kpd7,mon,slmask,
     &                    smc(1,lsoil,nn),len,iret
     &,                   imsk, jmsk, slmskh, gaus,blno, blto
     &,                   outlat, outlon, me)
              do l=1,lsoil-1
                do i = 1, len
                 smc(i,l,nn) = smc(i,lsoil,nn)
                enddo
              enddo
            else  ! the new gldas data.  it does not have data defined at landice
                  ! points.  so for efficiency, don't have fixrdc try to
                  ! find a value at landice points as defined by the vet type (vet).
              allocate(slmask_noice(len))
              slmask_noice=1.0
              do i = 1, len
                if (nint(vet(i)) < 1 .or.
     &              nint(vet(i)) == landice_cat) then
                  slmask_noice(i) = 0.0
                endif
              enddo
              do k = 1, lsoil
                if (k==1) kpd7=10     ! 0_10 cm    (pds octs 11 and 12)
                if (k==2) kpd7=2600   ! 10_40 cm
                if (k==3) kpd7=10340  ! 40_100 cm
                if (k==4) kpd7=25800  ! 100_200 cm
                call fixrdc(lugb,fnsmcc,kpdsmc,kpd7,mon,slmask_noice,
     &                      smc(1,k,nn),len,iret
     &,                     imsk, jmsk, slmskh, gaus,blno, blto
     &,                     outlat, outlon, me)
              enddo
              deallocate(slmask_noice)
            endif
          else
            write(6,*) 'climatological soil wetness file not given'
            call abort
          endif
!
!  soil temperature
!
          if(fnstcc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnstcc,kpdstc,kpd7,mon,slmask,
     &                  stc(1,lsoil,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
            do l=1,lsoil-1
              do i = 1, len
               stc(i,l,nn) = stc(i,lsoil,nn)
              enddo
            enddo
          endif
!
!  sea ice
!
          kpd7=-1
          if(fnacnc(1:8).ne.'        ') then
            call fixrdc(lugb,fnacnc,kpdacn,kpd7,mon,slmask,
     &                  acn(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
          elseif(fnaisc(1:8).ne.'        ') then
            call fixrdc(lugb,fnaisc,kpdais,kpd7,mon,slmask,
     &                  ais(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
          else
            write(6,*) 'climatological ice cover file not given'
            call abort
          endif
!
!  snow depth
!
          kpd7=-1
          call fixrdc(lugb,fnsnoc,kpdsno,kpd7,mon,slmask,
     &                sno(1,nn),len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
!
!  snow cover
!
          if(fnscvc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnscvc,kpdscv,kpd7,mon,slmask,
     &                  scv(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
            write(6,*) 'climatological snow cover read in.'
          endif
!
!  surface roughness
!
      if(fnzorc(1:3) == 'sib') then
        if (me == 0) then
          write(6,*) 'roughness length to be set from sib veg type'
        endif
      elseif(fnzorc(1:4) == 'igbp') then
        if (me == 0) then
          write(6,*) 'roughness length to be set from igbp veg type'
        endif
      else
        kpd7=-1
        call fixrdc(lugb,fnzorc,kpdzor,kpd7,mon,slmask,
     &              zor(1,nn),len,iret
     &,             imsk, jmsk, slmskh, gaus,blno, blto
     &,             outlat, outlon, me)
      endif
!
          do i = 1, len
!                           set clouds climatology to zero
            cvclm (i) = 0.
            cvbclm(i) = 0.
            cvtclm(i) = 0.
!
            cnpclm(i) = 0.  !set canopy water content climatology to zero
          enddo
!
!  vegetation cover
!
          if(fnvegc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnvegc,kpdveg,kpd7,mon,slmask,
     &                  veg(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
            if (me .eq. 0) write(6,*) 'climatological vegetation',
     &                                ' cover read in for mon=',mon
          endif

        enddo
!
      mon1s = mon1 ; mon2s = mon2 ; sea1s = sea1 ; sea2s = sea2
!
      if (me .eq. 0) print *,' mon1s=',mon1s,' mon2s=',mon2s
     &,' sea1s=',sea1s,' sea2s=',sea2s
!
        k1  = 1 ; k2 = 2
        m1  = 1 ; m2 = 2
!
        first = .false.
      endif first_time
!
!     to get tsf climatology at the previous call to sfccycle
!
!     if (fh-deltsfc >= 0.0) then
        rjdayh = rjday - deltsfc/24.0
!     else
!       rjdayh = rjday
!     endif
!     if(lprnt) print *,' rjdayh=',rjdayh,' mon1=',mon1,' mon2='
!    &,mon2,' mon1s=',mon1s,' mon2s=',mon2s,' k1=',k1,' k2=',k2
      if (rjdayh .ge. dayhf(mon1)) then
        if (mon2 .eq. 1) mon2 = 13
        wei1x = (dayhf(mon2)-rjdayh)/(dayhf(mon2)-dayhf(mon1))
        wei2x = 1.0 - wei1x
        if (mon2 .eq. 13) mon2 = 1
      else
        rjdayh2 = rjdayh
        if (rjdayh .lt. dayhf(1)) rjdayh2 = rjdayh2 + 365.0
        if (mon1s .eq. mon1) then
          mon1s = mon1 - 1
          if (mon1s .eq. 0) mon1s = 12
          k2  = k1
          k1  = mod(k2,2) + 1
          mon = mon1s
          kpd7=-1
          call fixrdc(lugb,fntsfc,kpdtsf,kpd7,mon,slmask,
     &               tsf(1,k1),len,iret
     &,              imsk, jmsk, slmskh, gaus,blno, blto
     &,              outlat, outlon, me)
        endif
        mon2s = mon1s + 1
!       if (mon2s .eq. 1) mon2s = 13
        wei1x = (dayhf(mon2s)-rjdayh2)/(dayhf(mon2s)-dayhf(mon1s))
        wei2x = 1.0 - wei1x
        if (mon2s .eq. 13) mon2s = 1
        do i=1,len
          tsf2(i) = wei1x * tsf(i,k1) + wei2x * tsf(i,k2)
        enddo
      endif
!
!cbosu new albedo is monthly
      if (sea1 .ne. sea1s) then
         sea1s = sea1
         sea2s = sea2
         m1    = mod(m1,2) + 1
         m2    = mod(m1,2) + 1
!
!     seasonal mean climatology
!
         isx   = sea2/3 + 1
         if (isx .eq. 5) isx = 1
         if(isx.eq.1) kpd9 = 12
         if(isx.eq.2) kpd9 = 3
         if(isx.eq.3) kpd9 = 6
         if(isx.eq.4) kpd9 = 9
!
!  albedo
!  there are four albedo fields in this version:
!  two for strong zeneith angle dependent (visible and near ir)
!  and two for weak zeneith angle dependent (vis ans nir)
!
!cbosu  
        if (ialb == 0) then
           kpd7=-1
           do k = 1, 4
             call fixrdc(lugb,fnalbc,kpdalb(k),kpd7,kpd9,slmask
     &,                 alb(1,k,m2),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
           enddo
        endif

      endif

      if (mon1 .ne. mon1s) then

         mon1s = mon1
         mon2s = mon2
         k1    = mod(k1,2) + 1
         k2    = mod(k1,2) + 1
!
!     monthly mean climatology
!
          mon = mon2
          nn  = k2
!cbosu
          if (ialb == 1) then
            if (me == 0) print*,'bosu 2nd time in clima for month ',
     &                   mon, k1,k2
            kpd7=-1
            do k = 1, 4
              call fixrdc(lugb,fnalbc,kpdalb(k),kpd7,mon,slmask,
     &                    alb(1,k,nn),len,iret
     &,                   imsk, jmsk, slmskh, gaus,blno, blto
     &,                   outlat, outlon, me)
            enddo
          endif
!
!  tsf at the current time t
!
          kpd7=-1
          call fixrdc(lugb,fntsfc,kpdtsf,kpd7,mon,slmask,
     &               tsf(1,nn),len,iret
     &,              imsk, jmsk, slmskh, gaus,blno, blto
     &,              outlat, outlon, me)
!
!  soil wetness
!
          if(fnwetc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnwetc,kpdwet,kpd7,mon,slmask,
     &                  wet(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
          elseif(fnsmcc(1:8).ne.'        ') then
            if (index(fnsmcc,'global_soilmcpc.1x1.grb') /= 0) then ! the old climo data
              kpd7=-1
              call fixrdc(lugb,fnsmcc,kpdsmc,kpd7,mon,slmask,
     &                    smc(1,lsoil,nn),len,iret
     &,                   imsk, jmsk, slmskh, gaus,blno, blto
     &,                   outlat, outlon, me)
              do l=1,lsoil-1
                do i = 1, len
                 smc(i,l,nn) = smc(i,lsoil,nn)
                enddo
              enddo
            else  ! the new gldas data.  it does not have data defined at landice
                  ! points.  so for efficiency, don't have fixrdc try to
                  ! find a value at landice points as defined by the vet type (vet).
              allocate(slmask_noice(len))
              slmask_noice=1.0
              do i = 1, len
                if (nint(vet(i)) < 1 .or.
     &              nint(vet(i)) == landice_cat) then
                  slmask_noice(i) = 0.0
                endif
              enddo
              do k = 1, lsoil
                if (k==1) kpd7=10     ! 0_10 cm   (pds octs 11 and 12)
                if (k==2) kpd7=2600   ! 10_40 cm
                if (k==3) kpd7=10340  ! 40_100 cm
                if (k==4) kpd7=25800  ! 100_200 cm
                call fixrdc(lugb,fnsmcc,kpdsmc,kpd7,mon,slmask_noice,
     &                      smc(1,k,nn),len,iret
     &,                     imsk, jmsk, slmskh, gaus,blno, blto
     &,                     outlat, outlon, me)
              enddo
              deallocate(slmask_noice)
            endif
          else
            write(6,*) 'climatological soil wetness file not given'
            call abort
          endif
!
!  sea ice
!
          kpd7=-1
          if(fnacnc(1:8).ne.'        ') then
            call fixrdc(lugb,fnacnc,kpdacn,kpd7,mon,slmask,
     &                  acn(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
          elseif(fnaisc(1:8).ne.'        ') then
            call fixrdc(lugb,fnaisc,kpdais,kpd7,mon,slmask,
     &                  ais(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
          else
            write(6,*) 'climatological ice cover file not given'
            call abort
          endif
!
!  snow depth
!
          kpd7=-1
          call fixrdc(lugb,fnsnoc,kpdsno,kpd7,mon,slmask,
     &                sno(1,nn),len,iret
     &,               imsk, jmsk, slmskh, gaus,blno, blto
     &,               outlat, outlon, me)
!
!  snow cover
!
          if(fnscvc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnscvc,kpdscv,kpd7,mon,slmask,
     &                  scv(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
            write(6,*) 'climatological snow cover read in.'
          endif
!
!  surface roughness
!
      if(fnzorc(1:3) == 'sib') then
        if (me == 0) then
          write(6,*) 'roughness length to be set from sib veg type'
        endif
      elseif(fnzorc(1:4) == 'igbp') then
        if (me == 0) then
          write(6,*) 'roughness length to be set from igbp veg type'
        endif
      else
        kpd7=-1
        call fixrdc(lugb,fnzorc,kpdzor,kpd7,mon,slmask,
     &              zor(1,nn),len,iret
     &,             imsk, jmsk, slmskh, gaus,blno, blto
     &,             outlat, outlon, me)
      endif
!
!  vegetation cover
!
          if(fnvegc(1:8).ne.'        ') then
            kpd7=-1
            call fixrdc(lugb,fnvegc,kpdveg,kpd7,mon,slmask,
     &                  veg(1,nn),len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
!           if (me .eq. 0) write(6,*) 'climatological vegetation',
!    &                                ' cover read in for mon=',mon
          endif
!
      endif
!
!     now perform the time interpolation
!
! when chosen, set the z0 based on the vegetation type.
! for this option to work, namelist variable fnvetc must be
! set to point at the proper vegetation type file.
      if(fnzorc(1:3) == 'sib') then
        if(fnvetc(1:4) == '   ') then
          if (me==0) write(6,*) "must choose sib veg type climo file"
          call abort
        endif
        zorclm = 0.0
        do i=1,len
          ivtyp=nint(vet(i))
          if (ivtyp >= 1 .and. ivtyp <= 13) then
            zorclm(i) = z0_sib(ivtyp)
          endif
        enddo
      elseif(fnzorc(1:4) == 'igbp') then
        if(fnvetc(1:4) == '   ') then
          if (me==0) write(6,*) "must choose igbp veg type climo file"
          call abort
        endif
        zorclm = 0.0
        do i=1,len
          ivtyp=nint(vet(i))
          if (ivtyp >= 1 .and. ivtyp <= 20) then
            z0_season(1) = z0_igbp_min(ivtyp)
            z0_season(7) = z0_igbp_max(ivtyp)
            if(outlat(i) < 0.0)then
              zorclm(i) = wei1y * z0_season(hyr2) + 
     &                    wei2y *z0_season(hyr1)
             else
              zorclm(i) = wei1y * z0_season(hyr1) + 
     &                    wei2y *z0_season(hyr2)
           endif
          endif
        enddo
      else
        do i=1,len
          zorclm(i) = wei1m * zor(i,k1) + wei2m * zor(i,k2)
        enddo
      endif
!
      do i=1,len
        tsfclm(i) = wei1m * tsf(i,k1) + wei2m * tsf(i,k2)
        snoclm(i) = wei1m * sno(i,k1) + wei2m * sno(i,k2)
        cvclm(i)  = 0.0
        cvbclm(i) = 0.0
        cvtclm(i) = 0.0
        cnpclm(i) = 0.0
        tsfcl2(i) = tsf2(i)
      enddo
!     if(lprnt) print *,' tsfclm=',tsfclm(iprnt),' wei1m=',wei1m
!    &,' wei2m=',wei2m,' tsfk12=',tsf(iprnt,k1),tsf(iprnt,k2)
!
      if (fh .eq. 0.0) then
        do i=1,len
          tsfcl0(i) = tsfclm(i)
        enddo
      endif
      if (rjdayh .ge. dayhf(mon1)) then
        do i=1,len
          tsf2(i)   = wei1x * tsf(i,k1) + wei2x * tsf(i,k2)
          tsfcl2(i) = tsf2(i)
        enddo
      endif
!     if(lprnt) print *,' tsf2=',tsf2(iprnt),' wei1x=',wei1x
!    &,' wei2x=',wei2x,' tsfk12=',tsf(iprnt,k1),tsf(iprnt,k2)
!    &,' mon1s=',mon1s,' mon2s=',mon2s
!    &,' slmask=',slmask(iprnt)
!
      if(fnacnc(1:8).ne.'        ') then
        do i=1,len
          acnclm(i) = wei1m * acn(i,k1) + wei2m * acn(i,k2)
        enddo
      elseif(fnaisc(1:8).ne.'        ') then
        do i=1,len
          aisclm(i) = wei1m * ais(i,k1) + wei2m * ais(i,k2)
        enddo
      endif
!
      if(fnwetc(1:8).ne.'        ') then
        do i=1,len
          wetclm(i) = wei1m * wet(i,k1) + wei2m * wet(i,k2)
        enddo
      elseif(fnsmcc(1:8).ne.'        ') then
        do k=1,lsoil
          do i=1,len
            smcclm(i,k) = wei1m * smc(i,k,k1) + wei2m * smc(i,k,k2)
          enddo
        enddo
      endif
!
      if(fnscvc(1:8).ne.'        ') then
        do i=1,len
          scvclm(i) = wei1m * scv(i,k1) + wei2m * scv(i,k2)
        enddo
      endif
!
      if(fntg3c(1:8).ne.'        ') then
        do i=1,len
          tg3clm(i) =         tg3(i)
        enddo
      elseif(fnstcc(1:8).ne.'        ') then
        do k=1,lsoil
          do i=1,len
            stcclm(i,k) = wei1m * stc(i,k,k1) + wei2m * stc(i,k,k2)
          enddo
        enddo
      endif
!
      if(fnvegc(1:8).ne.'        ') then
        do i=1,len
          vegclm(i) = wei1m * veg(i,k1) + wei2m * veg(i,k2)
        enddo
      endif
!
      if(fnvetc(1:8).ne.'        ') then
        do i=1,len
          vetclm(i) =         vet(i)
        enddo
      endif
!
      if(fnsotc(1:8).ne.'        ') then
        do i=1,len
          sotclm(i) =         sot(i)
        enddo
      endif


!clu ----------------------------------------------------------------------
!
      if(fnvmnc(1:8).ne.'        ') then
        do i=1,len
          vmnclm(i) =         vmn(i)
        enddo
      endif
!
      if(fnvmxc(1:8).ne.'        ') then
        do i=1,len
          vmxclm(i) =         vmx(i)
        enddo
      endif
!
      if(fnslpc(1:8).ne.'        ') then
        do i=1,len
          slpclm(i) =         slp(i)
        enddo
      endif
!
      if(fnabsc(1:8).ne.'        ') then
        do i=1,len
          absclm(i) =         abs(i)
        enddo
      endif
!clu ----------------------------------------------------------------------
!
!cbosu  diagnostic print
      if (me == 0) print*,'monthly albedo weights are ', 
     &             wei1m,' for k', k1, wei2m, ' for k', k2

      if (ialb == 1) then
        do k=1,4
          do i=1,len
            albclm(i,k) = wei1m * alb(i,k,k1) + wei2m * alb(i,k,k2)
          enddo
        enddo
      else
        do k=1,4
          do i=1,len
            albclm(i,k) = wei1s * alb(i,k,m1) + wei2s * alb(i,k,m2)
          enddo
        enddo
      endif
!
      do k=1,2
        do i=1,len
          alfclm(i,k) = alf(i,k)
        enddo
      enddo
!
!  end of climatology reads
!
      return
      end subroutine clima
      subroutine fixrdc(lugb,fngrib,kpds5,kpds7,mon,slmask,
     &                 gdata,len,iret
     &,                imsk, jmsk, slmskh, gaus,blno, blto
     &,                outlat, outlon, me)
      use machine ,      only : kind_io8,kind_io4
      use sfccyc_module, only : mdata
      implicit none
      integer imax,jmax,ijmax,i,j,n,jret,inttyp,iret,imsk,
     &        jmsk,len,lugb,kpds5,mon,lskip,lgrib,ndata,lugi,me,kmami
     &,       jj,w3kindreal,w3kindint
      real (kind=kind_io8) wlon,elon,rnlat,dlat,dlon,rslat,blno,blto
!
!   read in grib climatology files and interpolate to the input
!   grid.  grib files should allow all the necessary parameters
!   to be extracted from the description records.
!
!
      character*500 fngrib
!     character*80 fngrib, asgnstr
!
      real (kind=kind_io8) slmskh(imsk,jmsk)
!
      real (kind=kind_io8) gdata(len), slmask(len)
      real (kind=kind_io8), allocatable :: data(:,:), rslmsk(:,:)
      real (kind=kind_io8), allocatable :: data8(:)
      real (kind=kind_io4), allocatable :: data4(:)
      real (kind=kind_io8), allocatable :: rlngrb(:), rltgrb(:)
!
      logical lmask, yr2kc, gaus, ijordr
      logical, allocatable :: lbms(:)
!
      integer, intent(in) :: kpds7
      integer kpds(1000),kgds(1000)
      integer jpds(1000),jgds(1000), kpds0(1000)
      real (kind=kind_io8) outlat(len), outlon(len)
!
      allocate(data8(1:mdata))
      allocate(lbms(mdata))
!
!     integer imax_sv, jmax_sv, wlon_sv, rnlat_sv, kpds1_sv
!     date imax_sv/0/, jmax_sv/0/, wlon_sv/999.0/, rnlat_sv/999.0/
!    &,    kpds1_sv/-1/
!     save imax_sv, jmax_sv, wlon_sv, rnlat_sv, kpds1_sv
!    &,    rlngrb, rltgrb
!
      iret   = 0
!
      if (me .eq. 0) write(6,*) ' in fixrdc for mon=',mon
     &,' fngrib=',trim(fngrib)
!
      close(lugb)
      call baopenr(lugb,fngrib,iret)
      if (iret .ne. 0) then
        write(6,*) ' error in opening file ',trim(fngrib)
        print *,'error in opening file ',trim(fngrib)
        call abort
      endif
      if (me .eq. 0) write(6,*) ' file ',trim(fngrib),
     &             ' opened. unit=',lugb
!
      lugi = 0
!
      lskip   = -1
      jpds    = -1
      jgds    = -1
      jpds(5) = kpds5
      jpds(7) = kpds7
      kpds    = jpds
      call getgbh(lugb,lugi,lskip,jpds,jgds,lgrib,ndata,
     &            lskip,kpds,kgds,iret)
      if (me .eq. 0) then
      write(6,*) ' first grib record.'
      write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
      write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
      write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
      endif
      yr2kc     = (kpds(8) / 100) .gt. 0
      kpds0     = jpds
      kpds0(4)  = -1
      kpds0(18) = -1
      if(iret.ne.0) then
        write(6,*) ' error in getgbh. iret: ', iret
        if (iret==99) write(6,*) ' field not found.'
        call abort
      endif
!
!   handling climatology file
!
      lskip   = -1
      n       = 0
      jpds    = kpds0
      jpds(9) = mon
      if(jpds(9).eq.13) jpds(9) = 1
      call w3kind(w3kindreal,w3kindint)
      if (w3kindreal==8) then
        call getgb(lugb,lugi,mdata,lskip,jpds,jgds,ndata,lskip,
     &             kpds,kgds,lbms,data8,jret)
      else if (w3kindreal==4) then
        allocate(data4(1:mdata))
        call getgb(lugb,lugi,mdata,lskip,jpds,jgds,ndata,lskip,
     &             kpds,kgds,lbms,data4,jret)
        data8 = real(data4, kind=kind_io8)
        deallocate(data4)
      endif
      if (me .eq. 0) write(6,*) ' input grib file dates=',
     &              (kpds(i),i=8,11)
      if(jret.eq.0) then
        if(ndata.eq.0) then
          write(6,*) ' error in getgb'
          write(6,*) ' kpds=',kpds
          write(6,*) ' kgds=',kgds
          call abort
        endif
        imax=kgds(2)
        jmax=kgds(3)
        ijmax=imax*jmax
        allocate (data(imax,jmax))
        do j=1,jmax
          jj = (j-1)*imax
          do i=1,imax
            data(i,j) = data8(jj+i)
          enddo
        enddo
        if (me .eq. 0) write(6,*) 'imax,jmax,ijmax=',imax,jmax,ijmax
      else
        write(6,*) ' error in getgb - jret=', jret
        call abort
      endif
!
      if (me .eq. 0) then
      write(6,*) ' maxmin of input as is'
      kmami=1
      call maxmin(data(1,1),ijmax,kmami)
      endif
!
      call getarea(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr,me)
      if (me .eq. 0) then
      write(6,*) 'imax,jmax,ijmax,dlon,dlat,ijordr,wlon,rnlat='
      write(6,*)  imax,jmax,ijmax,dlon,dlat,ijordr,wlon,rnlat
      endif
      call subst(data,imax,jmax,dlon,dlat,ijordr)
!
!   first get slmask over input grid
!
        allocate (rlngrb(imax), rltgrb(jmax))
        allocate (rslmsk(imax,jmax))

        call setrmsk(kpds5,slmskh,imsk,jmsk,wlon,rnlat,
     &               data,imax,jmax,rlngrb,rltgrb,lmask,rslmsk
     &,                  gaus,blno, blto, kgds(1), kpds(4), lbms)
!       write(6,*) ' kpds5=',kpds5,' lmask=',lmask
!
                         inttyp = 0
        if(kpds5.eq.225) inttyp = 1
        if(kpds5.eq.230) inttyp = 1
        if(kpds5.eq.236) inttyp = 1  
        if(kpds5.eq.224) inttyp = 1  
        if (me .eq. 0) then
        if(inttyp.eq.1) print *, ' nearest grid point used'
     &,   ' kpds5=',kpds5, ' lmask = ',lmask
        endif
!
        call la2ga(data,imax,jmax,rlngrb,rltgrb,wlon,rnlat,inttyp,
     &             gdata,len,lmask,rslmsk,slmask
     &,            outlat, outlon,me)
!
        deallocate (rlngrb, stat=iret)
        deallocate (rltgrb, stat=iret)
        deallocate (data, stat=iret)
        deallocate (rslmsk, stat=iret)
      call baclose(lugb,iret)
!
      deallocate(data8)
      deallocate(lbms)
      return
      end subroutine fixrdc
      subroutine fixrda(lugb,fngrib,kpds5,slmask,
     &                  iy,im,id,ih,fh,gdata,len,iret
     &,                 imsk, jmsk, slmskh, gaus,blno, blto
     &,                 outlat, outlon, me)
      use machine      , only : kind_io8,kind_io4
      use sfccyc_module, only : mdata
      implicit none
      integer nrepmx,nvalid,imo,iyr,idy,jret,ihr,nrept,lskip,lugi,
     &        lgrib,j,ndata,i,inttyp,jmax,imax,ijmax,ij,jday,len,iret,
     &        jmsk,imsk,ih,kpds5,lugb,iy,id,im,jh,jd,jdoy,jdow,jm,me,
     &        monend,jy,iy4,kmami,iret2,jj,w3kindreal,w3kindint
      real (kind=kind_io8) rnlat,rslat,wlon,elon,dlon,dlat,fh,blno,
     &                     rjday,blto
!
!   read in grib climatology/analysis files and interpolate to the input
!   dates and the grid.  grib files should allow all the necessary parameters
!   to be extracted from the description records.
!
!  nrepmx:  max number of days for going back date search
!  nvalid:  analysis later than (current date - nvalid) is regarded as
!           valid for current analysis
!
      parameter(nrepmx=15, nvalid=4)
!
      character*500 fngrib
!     character*80 fngrib, asgnstr
!
      real (kind=kind_io8) slmskh(imsk,jmsk)
!
      real (kind=kind_io8) gdata(len), slmask(len)
      real (kind=kind_io8), allocatable :: data(:,:),rslmsk(:,:)
      real (kind=kind_io8), allocatable :: data8(:)
      real (kind=kind_io4), allocatable :: data4(:)
      real (kind=kind_io8), allocatable :: rlngrb(:), rltgrb(:)
!
      logical lmask, yr2kc, gaus, ijordr
      logical*1  lbms(mdata)
!
      integer kpds(1000),kgds(1000)
      integer jpds(1000),jgds(1000), kpds0(1000)
      real (kind=kind_io8) outlat(len), outlon(len)
!
! dayhf : julian day of the middle of each month
!
      real (kind=kind_io8) dayhf(13)
      data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0,
     &           196.5,227.5,258.0,288.5,319.0,349.5,380.5/
!
! mjday : number of days in a month
!
      integer mjday(12)
      data mjday/31,28,31,30,31,30,31,31,30,31,30,31/
!
      real (kind=kind_io8) fha(5)
      real(4) fha4(5)
      integer ida(8),jda(8)
!
      allocate(data8(1:mdata))
      iret   = 0
      monend = 9999
!
!  compute jy,jm,jd,jh of forecast and the day of the year
!
      iy4=iy
      if(iy.lt.101) iy4=1900+iy4
      fha=0
      ida=0
      jda=0
      fha(2)=nint(fh)
      ida(1)=iy
      ida(2)=im
      ida(3)=id
      ida(5)=ih
      call w3kind(w3kindreal,w3kindint)
      if(w3kindreal==4) then
        fha4=fha
        call w3movdat(fha4,ida,jda)
      else
        call w3movdat(fha,ida,jda)
      endif
      jy=jda(1)
      jm=jda(2)
      jd=jda(3)
      jh=jda(5)
!     if (me .eq. 0) write(6,*) ' forecast jy,jm,jd,jh,rjday=',
!    &               jy,jm,jd,jh,rjday
      jdow = 0
      jdoy = 0
      jday = 0
      call w3doxdat(jda,jdow,jdoy,jday)
      rjday=jdoy+jda(5)/24.
      if(rjday.lt.dayhf(1)) rjday=rjday+365.

      if (me .eq. 0) write(6,*) ' forecast jy,jm,jd,jh,rjday=',
     &               jy,jm,jd,jh,rjday
!
      if (me .eq. 0) then
      write(6,*) 'forecast jy,jm,jd,jh=',jy,jm,jd,jh
!
      write(6,*) ' '
      write(6,*) '************************************************'
      endif
!
      close(lugb)
      call baopenr(lugb,fngrib,iret)
      if (iret .ne. 0) then
        write(6,*) ' error in opening file ',trim(fngrib)
        print *,'error in opening file ',trim(fngrib)
        call abort
      endif
      if (me .eq. 0) write(6,*) ' file ',trim(fngrib),
     &             ' opened. unit=',lugb
!
      lugi = 0
!
      lskip=-1
      jpds=-1
      jgds=-1
      jpds(5)=kpds5
      kpds = jpds
      call getgbh(lugb,lugi,lskip,jpds,jgds,lgrib,ndata,
     &            lskip,kpds,kgds,iret)
      if (me .eq. 0) then
      write(6,*) ' first grib record.'
      write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
      write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
      write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
      endif
      yr2kc = (kpds(8) / 100) .gt. 0
      kpds0=jpds
      kpds0(4)=-1
      kpds0(18)=-1
      if(iret.ne.0) then
        write(6,*) ' error in getgbh. iret: ', iret
        if(iret==99) write(6,*) ' field not found.'
        call abort
      endif
!
!  handling analysis file
!
!  find record for the given hour/day/month/year
!
      nrept=0
      jpds=kpds0
      lskip = -1
      iyr=jy
      if(iyr.le.100) iyr=2050-mod(2050-iyr,100)
      imo=jm
      idy=jd
      ihr=jh
!     year 2000 compatible data
      if (yr2kc) then
         jpds(8) = iyr
      else
         jpds(8) = mod(iyr,1900)
      endif
   50 continue
      jpds( 8)=mod(iyr-1,100)+1
      jpds( 9)=imo
      jpds(10)=idy
!     jpds(11)=ihr
      jpds(21)=(iyr-1)/100+1
      call w3kind(w3kindreal,w3kindint)
      if (w3kindreal == 8) then
        call getgb(lugb,lugi,mdata,lskip,jpds,jgds,ndata,lskip,
     &             kpds,kgds,lbms,data8,jret)
      elseif (w3kindreal == 4) then
        allocate (data4(1:mdata))
        call getgb(lugb,lugi,mdata,lskip,jpds,jgds,ndata,lskip,
     &             kpds,kgds,lbms,data4,jret)
        data8 = real(data4, kind=kind_io8)
        deallocate(data4)
      endif
      if (me .eq. 0) write(6,*) ' input grib file dates=',
     &              (kpds(i),i=8,11)
      if(jret.eq.0) then
        if(ndata.eq.0) then
          write(6,*) ' error in getgb'
          write(6,*) ' kpds=',kpds
          write(6,*) ' kgds=',kgds
          call abort
        endif
        imax=kgds(2)
        jmax=kgds(3)
        ijmax=imax*jmax
        allocate (data(imax,jmax))
        do j=1,jmax
          jj = (j-1)*imax
          do i=1,imax
            data(i,j) = data8(jj+i)
          enddo
        enddo
      else
        if(nrept.eq.0) then
          if (me .eq. 0) then
          write(6,*) ' no matching dates found.  start searching',
     &               ' nearest matching dates (going back).'
          endif
        endif
!
!  no matching ih found. search nearest hour
!
        if(ihr.eq.6) then
          ihr=0
          go to 50
        elseif(ihr.eq.12) then
          ihr=0
          go to 50
        elseif(ihr.eq.18) then
          ihr=12
          go to 50
        elseif(ihr.eq.0.or.ihr.eq.-1) then
          idy=idy-1
          if(idy.eq.0) then
            imo=imo-1
            if(imo.eq.0) then
              iyr=iyr-1
              if(iyr.lt.0) iyr=99
              imo=12
            endif
            idy=31
            if(imo.eq.4.or.imo.eq.6.or.imo.eq.9.or.imo.eq.11) idy=30
            if(imo.eq.2) then
              if(mod(iyr,4).eq.0) then
                idy=29
              else
                idy=28
              endif
            endif
          endif
          ihr=-1
          if (me .eq. 0) write(6,*) ' decremented dates=',
     &                              iyr,imo,idy,ihr
          nrept=nrept+1
          if(nrept.gt.nvalid) iret=-1
          if(nrept.gt.nrepmx) then
            if (me .eq. 0) then
              write(6,*) ' <warning:cycl> searching range exceeded.'
     &,                  ' may be wrong grib file given'
              write(6,*) ' <warning:cycl> fngrib=',trim(fngrib)
              write(6,*) ' <warning:cycl> terminating search and',
     &                   ' and setting gdata to -999'
              write(6,*) ' range max=',nrepmx
            endif
!           imax=kgds(2)
!           jmax=kgds(3)
!           ijmax=imax*jmax
!           do ij=1,ijmax
!             data(ij)=0.
!           enddo
            go to 100
          endif
          go to 50
        else
          if (me .eq. 0) then
            write(6,*) ' search of analysis for ihr=',ihr,' failed.'
            write(6,*) ' kpds=',kpds
            write(6,*) ' iyr,imo,idy,ihr=',iyr,imo,idy,ihr
          endif
          go to 100
        endif
      endif
!
   80 continue
      if (me .eq. 0) then
      write(6,*) ' maxmin of input as is'
      kmami=1
      call maxmin(data(1,1),ijmax,kmami)
      endif
!
      call getarea(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr,me)
      if (me .eq. 0) then
      write(6,*) 'imax,jmax,ijmax,dlon,dlat,ijordr,wlon,rnlat='
      write(6,*)  imax,jmax,ijmax,dlon,dlat,ijordr,wlon,rnlat
      endif
      call subst(data,imax,jmax,dlon,dlat,ijordr)
!
!   first get slmask over input grid
!
        allocate (rlngrb(imax), rltgrb(jmax))
        allocate (rslmsk(imax,jmax))
        call setrmsk(kpds5,slmskh,imsk,jmsk,wlon,rnlat,
     &               data,imax,jmax,rlngrb,rltgrb,lmask,rslmsk
!    &               data,imax,jmax,abs(dlon),abs(dlat),lmask,rslmsk
!cggg     &,                  gaus,blno, blto, kgds(1))
     &,                  gaus,blno, blto, kgds(1), kpds(4), lbms)

!       write(6,*) ' kpds5=',kpds5,' lmask=',lmask
!
                         inttyp = 0
        if(kpds5.eq.225) inttyp = 1
        if(kpds5.eq.230) inttyp = 1
        if(kpds5.eq.66)  inttyp = 1
        if(inttyp.eq.1) print *, ' nearest grid point used'
!
        call la2ga(data,imax,jmax,rlngrb,rltgrb,wlon,rnlat,inttyp,
     &             gdata,len,lmask,rslmsk,slmask
     &,            outlat, outlon, me)
!
      deallocate (rlngrb, stat=iret)
      deallocate (rltgrb, stat=iret)
      deallocate (data, stat=iret)
      deallocate (rslmsk, stat=iret)
      call baclose(lugb,iret2)
!     write(6,*) ' '
      deallocate(data8)
      return
!
  100 continue
      iret=1
      do i=1,len
        gdata(i) = -999.
      enddo
!
      call baclose(lugb,iret2)
!
      deallocate(data8)
      return
      end subroutine fixrda
      subroutine snodpth2(glacir,snwmax,snoanl, len, me)
      use machine , only : kind_io8,kind_io4
      implicit none
      integer i,me,len
      real (kind=kind_io8) snwmax
!
      real (kind=kind_io8) snoanl(len), glacir(len)
!
      if (me .eq. 0) write(6,*) 'snodpth2'
!
      do i=1,len
!
!  if glacial points has snow in climatology, set sno to snomax
!
        if(glacir(i).ne.0..and.snoanl(i).lt.snwmax*0.5) then
            snoanl(i) = snwmax + snoanl(i)
        endif
!
      enddo
      return
      end
