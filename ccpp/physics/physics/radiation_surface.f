!>  \file radiation_surface.f
!!  This file contains routines that set up surface albedo for SW
!!  radiation and surface emissivity for LW radiation.

!  ==========================================================  !!!!!
!            'module_radiation_surface' description            !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!    this module sets up surface albedo for sw radiation and surface   !
!    emissivity for lw radiation.                                      !
!                                                                      !
!                                                                      !
!    in the module, the externally callabe subroutines are :           !
!                                                                      !
!      'sfc_init'   -- initialization radiation surface data           !
!         inputs:                                                      !
!           ( me )                                                     !
!         outputs:                                                     !
!           (none)                                                     !
!                                                                      !
!      'setalb'     -- set up four-component surface albedoes          !
!         inputs:                                                      !
!           (slmsk,snowf,sncovr,snoalb,zorlf,coszf,tsknf,tairf,hprif,  !
!            alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc            !
!            IMAX)                                                     !
!         outputs:                                                     !
!           (sfcalb)                                                   !
!                                                                      !
!      'setemis'    -- set up surface emissivity for lw radiation      !
!         inputs:                                                      !
!           (xlon,xlat,slmsk,snowf,sncovr,zorlf,tsknf,tairf,hprif,     !
!            IMAX)                                                     !
!         outputs:                                                     !
!           (sfcemis)                                                  !
!                                                                      !
!    external modules referenced:                                      !
!                                                                      !
!       'module machine'             in 'machine.f'                    !
!       'module physcons'            in 'physcons.f'                   !
!       'module module_iounitdef'    in 'iounitdef.f'                  !
!                                                                      !
!                                                                      !
!    program history log:                                              !
!           1995   y.t. hou     - created albaer.f (include albedo     !
!                                 and aerosols calculations)           !
!      nov  1997   y.t. hou     - modified snow albedo                 !
!      jan  1998   y.t. hou     - included grumbine's sea-ice scheme   !
!      feb  1998   h.l. pan     - seasonal interpolation in cycle      !
!      mar  2000   y.t. hou     - modified to use opac aerosol data    !
!      apr  2003   y.t. hou     - seperate albedo and aerosols into    !
!                    two subroutines, rewritten in f90 modulized form  !
!      jan  2005   s. moorthi   - xingren's sea-ice fraction added     !
!      apr  2005   y.t. hou     - revised module structure             !
!      feb  2006   y.t. hou     - add varying surface emissivity,      !
!                    modified sfc albedo structure for modis shceme    !
!      Mar  2006   s. moorthi   - added surface temp over ice fraction !
!      mar  2007   c. marshall & h. wei                                !
!                               - added modis based sfc albedo scheme  !
!      may  2007   y. hou & s. moorthi                                 !
!                               - fix bug in modis albedo over ocean   !
!      aug  2007   h. wei & s. moorthi                                 !
!                               - fix bug in modis albedo over sea-ice !
!      aug  2007   y. hou       - fix bug in emissivity over ocean in  !
!                                 the modis scheme option              !
!      dec  2008   f. yang      - modified zenith angle dependence on  !
!                                 surface albedo over land. (2008 jamc)!
!      aug  2012   y. hou       - minor modification in initialization !
!                                 subr 'sfc_init'.                     !
!      nov  2012   y. hou       - modified control parameters through  !
!                    module 'physparam'.                               !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!


!> \defgroup module_radiation_surface RRTMG Surface Module
!! @{
!> This module sets up surface albedo for sw radiation and surface
!! emissivity for lw radiation.
!!\version NCEP-Radiation_surface   v5.1  Nov 2012
!!
!! In the module, the externally callabe subroutines are : 
!! + sfc_init(): initialization radiation surface data     
!! + setalb(): set up four-component surface albedoes    
!! + setemis(): set up surface emissivity for lw radiation
!========================================!
      module module_radiation_surface    !
!........................................!
!
      use physparam,         only : ialbflg, iemsflg, semis_file,       &
     &                              kind_phys
      use physcons,          only : con_t0c, con_ttp, con_pi, con_tice
      use module_iounitdef,  only : NIRADSF
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGSFC='NCEP-Radiation_surface   v5.1  Nov 2012 '
!    &   VTAGSFC='NCEP-Radiation_surface   v5.0  Aug 2012 '

!  ---  constant parameters
      integer, parameter, public :: NF_ALBD = 4     !< number of surface albedo components
      integer, parameter, public :: IMXEMS = 360    !< number of longtitude points in global emis-type map
      integer, parameter, public :: JMXEMS = 180    !< number of latitude points in global emis-type map
      real (kind=kind_phys), parameter :: f_zero = 0.0
      real (kind=kind_phys), parameter :: f_one  = 1.0
      real (kind=kind_phys), parameter :: rad2dg= 180.0 / con_pi
      integer, allocatable  ::  idxems(:,:)         !< global surface emissivity index array
      integer :: iemslw = 0                         !< global surface emissivity control flag set up in 'sfc_init'
!
      public  sfc_init, setalb, setemis

! =================
      contains
! =================

!> \ingroup module_radiation_surface
!> This subroutine is the initialization program for surface radiation
!! related quantities (albedo, emissivity, etc.)
!!\param me       print control flag
!>\section gen_sfc_init sfc_init General Algorithm
!! @{
!-----------------------------------
      subroutine sfc_init                                               &
     &     ( me )!  ---  inputs:
!  ---  outputs: ( none )

!  ===================================================================  !
!                                                                       !
!  this program is the initialization program for surface radiation     !
!  related quantities (albedo, emissivity, etc.)                        !
!                                                                       !
! usage:         call sfc_init                                          !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!      me           - print control flag                                !
!                                                                       !
!  outputs: (none) to module variables only                             !
!                                                                       !
!  external module variables:                                           !
!     ialbflg       - control flag for surface albedo schemes           !
!                     =0: climatology, based on surface veg types       !
!                     =1:                                               !
!     iemsflg       - control flag for sfc emissivity schemes (ab:2-dig)!
!                     a:=0 set sfc air/ground t same for lw radiation   !
!                       =1 set sfc air/ground t diff for lw radiation   !
!                     b:=0 use fixed sfc emissivity=1.0 (black-body)    !
!                       =1 use varying climtology sfc emiss (veg based) !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: me

!  ---  outputs: ( none )

!  ---  locals:
      integer    :: i, k
!     integer    :: ia, ja
      logical    :: file_exist
      character  :: cline*80
!
!===> ...  begin here
!
      if ( me == 0 ) print *, VTAGSFC   ! print out version tag

!> - Initialization of surface albedo section
!! \n physparam::ialbflg
!!  - = 0: using climatology surface albedo scheme for SW
!!  - = 1: using MODIS based land surface albedo for SW

      if ( ialbflg == 0 ) then

        if ( me == 0 ) then
          print *,' - Using climatology surface albedo scheme for sw'
        endif

      else if ( ialbflg == 1 ) then

        if ( me == 0 ) then
          print *,' - Using MODIS based land surface albedo for sw'
        endif

      else
        print *,' !! ERROR in Albedo Scheme Setting, IALB=',ialbflg
        stop
      endif    ! end if_ialbflg_block

!> - Initialization of surface emissivity section
!! \n physparam::iemsflg
!!  - = 0: fixed SFC emissivity at 1.0
!!  - = 1: input SFC emissivity type map from "semis_file"

      iemslw = mod(iemsflg, 10)          ! emissivity control
      if ( iemslw == 0 ) then            ! fixed sfc emis at 1.0

        if ( me == 0 ) then
          print *,' - Using Fixed Surface Emissivity = 1.0 for lw'
        endif

      elseif ( iemslw == 1 ) then        ! input sfc emiss type map

!  ---  allocate data space
        if ( .not. allocated(idxems) ) then
          allocate ( idxems(IMXEMS,JMXEMS)    )
        endif

!  ---  check to see if requested emissivity data file existed

        inquire (file=semis_file, exist=file_exist)

        if ( .not. file_exist ) then
          if ( me == 0 ) then
            print *,' - Using Varying Surface Emissivity for lw'
            print *,'   Requested data file "',semis_file,'" not found!'
            print *,'   Change to fixed surface emissivity = 1.0 !'
          endif

          iemslw = 0
        else
          close(NIRADSF)
          open (NIRADSF,file=semis_file,form='formatted',status='old')
          rewind NIRADSF

          read (NIRADSF,12) cline
  12      format(a80)

          read (NIRADSF,14) idxems
  14      format(80i1)

          if ( me == 0 ) then
            print *,' - Using Varying Surface Emissivity for lw'
            print *,'   Opened data file: ',semis_file
            print *, cline
!check      print *,' CHECK: Sample emissivity index data'
!           ia = IMXEMS / 5
!           ja = JMXEMS / 5
!           print *, idxems(1:IMXEMS:ia,1:JMXEMS:ja)
          endif

          close(NIRADSF)
        endif    ! end if_file_exist_block

      else
        print *,' !! ERROR in Emissivity Scheme Setting, IEMS=',iemsflg
        stop
      endif   ! end if_iemslw_block

!
      return
!...................................
      end subroutine sfc_init
!-----------------------------------
!! @}

!> \ingroup module_radiation_surface
!> This subroutine computes four components of surface albedos (i.e.,
!! vis-nir, direct-diffused) according to control flag ialbflg.
!! \n 1) climatological surface albedo scheme (\cite briegleb_1992)
!! \n 2) MODIS retrieval based scheme from Boston univ.
!!\param slmsk      (IMAX), sea(0),land(1),ice(2) mask on fcst model grid
!!\param snowf      (IMAX), snow depth water equivalent in mm
!!\param sncovr     (IMAX), snow cover over land
!!\param snoalb     (IMAX), maximum snow albedo over land (for deep snow)
!!\param zorlf      (IMAX), surface roughness in cm
!!\param coszf      (IMAX), cosin of solar zenith angle
!!\param tsknf      (IMAX), ground surface temperature in K
!!\param tairf      (IMAX), lowest model layer air temperature in K
!!\param hprif      (IMAX), topographic sdv in m
!!\n ---  for ialbflg=0 climtological albedo scheme  ---
!!\param alvsf      (IMAX), 60 degree vis albedo with strong cosz dependency
!!\param alnsf      (IMAX), 60 degree nir albedo with strong cosz dependency
!!\param alvwf      (IMAX), 60 degree vis albedo with weak cosz dependency
!!\param alnwf      (IMAX), 60 degree nir albedo with weak cosz dependency
!!\n ---  for ialbflg=1 MODIS based land albedo scheme ---
!!\param alvsf      (IMAX), visible black sky albedo at zenith 60 degree
!!\param alnsf      (IMAX), near-ir black sky albedo at zenith 60 degree
!!\param alvwf      (IMAX), visible white sky albedo
!!\param alnwf      (IMAX), near-ir white sky albedo
!!\param facsf      (IMAX), fractional coverage with strong cosz dependency
!!\param facwf      (IMAX), fractional coverage with weak cosz dependency
!!\param fice       (IMAX), sea-ice fraction
!!\param tisfc      (IMAX), sea-ice surface temperature
!!\param IMAX       array horizontal dimension
!!\param sfcalb     (IMAX,NF_ALBD), mean sfc albedo
!!\n                    ( :, 1) -     near ir direct beam albedo
!!\n                    ( :, 2) -     near ir diffused albedo
!!\n                    ( :, 3) -     uv+vis direct beam albedo
!!\n                    ( :, 4) -     uv+vis diffused albedo
!>\section general_setalb setalb General Algorithm
!! @{
!-----------------------------------
      subroutine setalb                                                 &
     &     ( slmsk,snowf,sncovr,snoalb,zorlf,coszf,tsknf,tairf,hprif,   & !  ---  inputs:
     &       alvsf,alnsf,alvwf,alnwf,facsf,facwf,fice,tisfc,            &
     &       IMAX,                                                      &
     &       sfcalb                                                     & !  ---  outputs:
     &     )

!  ===================================================================  !
!                                                                       !
!  this program computes four components of surface albedos (i.e.       !
!  vis-nir, direct-diffused) according to controflag ialbflg.           !
!   1) climatological surface albedo scheme (briegleb 1992)             !
!   2) modis retrieval based scheme from boston univ.                   !
!                                                                       !
!                                                                       !
! usage:         call setalb                                            !
!                                                                       !
! subprograms called:  none                                             !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     slmsk (IMAX)  - sea(0),land(1),ice(2) mask on fcst model grid     !
!     snowf (IMAX)  - snow depth water equivalent in mm                 !
!     sncovr(IMAX)  - ialgflg=0: not used                               !
!                     ialgflg=1: snow cover over land in fraction       !
!     snoalb(IMAX)  - ialbflg=0: not used                               !
!                     ialgflg=1: max snow albedo over land in fraction  !
!     zorlf (IMAX)  - surface roughness in cm                           !
!     coszf (IMAX)  - cosin of solar zenith angle                       !
!     tsknf (IMAX)  - ground surface temperature in k                   !
!     tairf (IMAX)  - lowest model layer air temperature in k           !
!     hprif (IMAX)  - topographic sdv in m                              !
!           ---  for ialbflg=0 climtological albedo scheme  ---         !
!     alvsf (IMAX)  - 60 degree vis albedo with strong cosz dependency  !
!     alnsf (IMAX)  - 60 degree nir albedo with strong cosz dependency  !
!     alvwf (IMAX)  - 60 degree vis albedo with weak cosz dependency    !
!     alnwf (IMAX)  - 60 degree nir albedo with weak cosz dependency    !
!           ---  for ialbflg=1 modis based land albedo scheme ---       !
!     alvsf (IMAX)  - visible black sky albedo at zenith 60 degree      !
!     alnsf (IMAX)  - near-ir black sky albedo at zenith 60 degree      !
!     alvwf (IMAX)  - visible white sky albedo                          !
!     alnwf (IMAX)  - near-ir white sky albedo                          !
!                                                                       !
!     facsf (IMAX)  - fractional coverage with strong cosz dependency   !
!     facwf (IMAX)  - fractional coverage with weak cosz dependency     !
!     fice  (IMAX)  - sea-ice fraction                                  !
!     tisfc (IMAX)  - sea-ice surface temperature                       !
!     IMAX          - array horizontal dimension                        !
!                                                                       !
!  outputs:                                                             !
!     sfcalb(IMAX,NF_ALBD)                                              !
!           ( :, 1) -     near ir direct beam albedo                    !
!           ( :, 2) -     near ir diffused albedo                       !
!           ( :, 3) -     uv+vis direct beam albedo                     !
!           ( :, 4) -     uv+vis diffused albedo                        !
!                                                                       !
!  module internal control variables:                                   !
!     ialbflg       - =0 use the default climatology surface albedo     !
!                     =1 use modis retrieved albedo and input snow cover!
!                        for land areas                                 !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer, intent(in) :: IMAX

      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       slmsk, snowf, zorlf, coszf, tsknf, tairf, hprif,           &
     &       alvsf, alnsf, alvwf, alnwf, facsf, facwf, fice, tisfc,     &
     &       sncovr, snoalb

!  ---  outputs
      real (kind=kind_phys), dimension(IMAX,NF_ALBD), intent(out) ::    &
     &       sfcalb
!     real (kind=kind_phys), dimension(:,:), intent(out) :: sfcalb

!  ---  locals:
      real (kind=kind_phys) :: asnvb, asnnb, asnvd, asnnd, asevb        &
     &,     asenb, asevd, asend, fsno,  fsea,  rfcs,  rfcw,  flnd       &
     &,     asnow, argh,  hrgh,  fsno0, fsno1, flnd0, fsea0, csnow      &
     &,     a1, a2, b1, b2, b3, ab1bm, ab2bm

      real (kind=kind_phys) ffw, dtgd

      integer :: i, k

!
!===> ...  begin here
!

!> -# If use climatological albedo scheme:
      if ( ialbflg == 0 ) then   ! use climatological albedo scheme

        do i = 1, IMAX

!>    - Modified snow albedo scheme - units convert to m (originally
!!      snowf in mm; zorlf in cm)

         asnow = 0.02*snowf(i)
         argh  = min(0.50, max(.025, 0.01*zorlf(i)))
         hrgh  = min(f_one, max(0.20, 1.0577-1.1538e-3*hprif(i) ) )
         fsno0 = asnow / (argh + asnow) * hrgh
         if (nint(slmsk(i))==0 .and. tsknf(i)>con_tice) fsno0 = f_zero
         fsno1 = f_one - fsno0
         flnd0 = min(f_one, facsf(i)+facwf(i))
         fsea0 = max(f_zero, f_one-flnd0)
         fsno  = fsno0
         fsea  = fsea0 * fsno1
         flnd  = flnd0 * fsno1

!>    - Calculate diffused sea surface albedo

         if (tsknf(i) >= 271.5) then
            asevd = 0.06
            asend = 0.06
         elseif (tsknf(i) < 271.1) then
            asevd = 0.70
            asend = 0.65
         else
            a1 = (tsknf(i) - 271.1)**2
            asevd = 0.7 - 4.0*a1
            asend = 0.65 - 3.6875*a1
         endif

!>    - Calculate diffused snow albedo.

         if (nint(slmsk(i)) == 2) then
            ffw   = f_one - fice(i)
            if (ffw < f_one) then
               dtgd = max(f_zero, min(5.0, (con_ttp-tisfc(i)) ))
               b1   = 0.03 * dtgd
            else
               b1 = f_zero
            endif

            b3   = 0.06 * ffw
            asnvd = (0.70 + b1) * fice(i) + b3
            asnnd = (0.60 + b1) * fice(i) + b3
            asevd = 0.70        * fice(i) + b3
            asend = 0.60        * fice(i) + b3
         else
            asnvd = 0.90
            asnnd = 0.75
         endif

!>    - Calculate direct snow albedo.

         if (coszf(i) < 0.5) then
            csnow = 0.5 * (3.0 / (f_one+4.0*coszf(i)) - f_one)
            asnvb = min( 0.98, asnvd+(1.0-asnvd)*csnow )
            asnnb = min( 0.98, asnnd+(1.0-asnnd)*csnow )
         else
            asnvb = asnvd
            asnnb = asnnd
         endif

!>    - Calculate direct sea surface albedo.

         if (coszf(i) > 0.0001) then
!           rfcs = 1.4 / (f_one + 0.8*coszf(i))
!           rfcw = 1.3 / (f_one + 0.6*coszf(i))
            rfcs = 2.14 / (f_one + 1.48*coszf(i))
            rfcw = rfcs

            if (tsknf(i) >= con_t0c) then
              asevb = max(asevd, 0.026/(coszf(i)**1.7+0.065)            &
     &              + 0.15 * (coszf(i)-0.1) * (coszf(i)-0.5)            &
     &              * (coszf(i)-f_one))
              asenb = asevb
            else
              asevb = asevd
              asenb = asend
            endif
         else
            rfcs  = f_one
            rfcw  = f_one
            asevb = asevd
            asenb = asend
         endif

         a1   = alvsf(i) * facsf(i)
         b1   = alvwf(i) * facwf(i)
         a2   = alnsf(i) * facsf(i)
         b2   = alnwf(i) * facwf(i)
         ab1bm = a1*rfcs + b1*rfcw
         ab2bm = a2*rfcs + b2*rfcw
         sfcalb(i,1) = min(0.99, ab2bm) *flnd + asenb*fsea + asnnb*fsno
         sfcalb(i,2) = (a2 + b2) * 0.96 *flnd + asend*fsea + asnnd*fsno
         sfcalb(i,3) = min(0.99, ab1bm) *flnd + asevb*fsea + asnvb*fsno
         sfcalb(i,4) = (a1 + b1) * 0.96 *flnd + asevd*fsea + asnvd*fsno

        enddo    ! end_do_i_loop

!> -# If use modis based albedo for land area:
      else

        do i = 1, IMAX

!>    - Calculate snow cover input directly for land model, no
!!      conversion needed.

         fsno0 = sncovr(i)

         if (nint(slmsk(i))==0 .and. tsknf(i)>con_tice) fsno0 = f_zero

         if (nint(slmsk(i)) == 2) then
           asnow = 0.02*snowf(i)
           argh  = min(0.50, max(.025, 0.01*zorlf(i)))
           hrgh  = min(f_one, max(0.20, 1.0577-1.1538e-3*hprif(i) ) )
           fsno0 = asnow / (argh + asnow) * hrgh
         endif

         fsno1 = f_one - fsno0
         flnd0 = min(f_one, facsf(i)+facwf(i))
         fsea0 = max(f_zero, f_one-flnd0)
         fsno  = fsno0
         fsea  = fsea0 * fsno1
         flnd  = flnd0 * fsno1

!>    - Calculate diffused sea surface albedo.

         if (tsknf(i) >= 271.5) then
            asevd = 0.06
            asend = 0.06
         elseif (tsknf(i) < 271.1) then
            asevd = 0.70
            asend = 0.65
         else
            a1 = (tsknf(i) - 271.1)**2
            asevd = 0.7 - 4.0*a1
            asend = 0.65 - 3.6875*a1
         endif

!>    - Calculate diffused snow albedo, land area use input max snow
!!      albedo.

         if (nint(slmsk(i)) == 2) then
            ffw   = f_one - fice(i)
            if (ffw < f_one) then
               dtgd = max(f_zero, min(5.0, (con_ttp-tisfc(i)) ))
               b1   = 0.03 * dtgd
            else
               b1 = f_zero
            endif

            b3   = 0.06 * ffw
            asnvd = (0.70 + b1) * fice(i) + b3
            asnnd = (0.60 + b1) * fice(i) + b3
            asevd = 0.70        * fice(i) + b3
            asend = 0.60        * fice(i) + b3
         else
            asnvd = snoalb(i)
            asnnd = snoalb(i)
         endif

!>    - Calculate direct snow albedo.

         if (nint(slmsk(i)) == 2) then
           if (coszf(i) < 0.5) then
              csnow = 0.5 * (3.0 / (f_one+4.0*coszf(i)) - f_one)
              asnvb = min( 0.98, asnvd+(f_one-asnvd)*csnow )
              asnnb = min( 0.98, asnnd+(f_one-asnnd)*csnow )
           else
             asnvb = asnvd
             asnnb = asnnd
           endif
         else
           asnvb = snoalb(i)
           asnnb = snoalb(i)
         endif

!>    - Calculate direct sea surface albedo, use fanglin's zenith angle
!!      treatment.

         if (coszf(i) > 0.0001) then

!           rfcs = 1.89 - 3.34*coszf(i) + 4.13*coszf(i)*coszf(i)        &
!    &           - 2.02*coszf(i)*coszf(i)*coszf(i)
            rfcs = 1.775/(1.0+1.55*coszf(i))

            if (tsknf(i) >= con_t0c) then
              asevb = max(asevd, 0.026/(coszf(i)**1.7+0.065)            &
     &              + 0.15 * (coszf(i)-0.1) * (coszf(i)-0.5)            &
     &              * (coszf(i)-f_one))
              asenb = asevb
            else
              asevb = asevd
              asenb = asend
            endif
         else
            rfcs  = f_one
            asevb = asevd
            asenb = asend
         endif

         ab1bm = min(0.99, alnsf(i)*rfcs)
         ab2bm = min(0.99, alvsf(i)*rfcs)
         sfcalb(i,1) = ab1bm   *flnd + asenb*fsea + asnnb*fsno
         sfcalb(i,2) = alnwf(i)     *flnd + asend*fsea + asnnd*fsno
         sfcalb(i,3) = ab2bm   *flnd + asevb*fsea + asnvb*fsno
         sfcalb(i,4) = alvwf(i)     *flnd + asevd*fsea + asnvd*fsno

        enddo    ! end_do_i_loop

      endif   ! end if_ialbflg
!
      return
!...................................
      end subroutine setalb
!-----------------------------------
!! @}

!> \ingroup module_radiation_surface
!> This subroutine computes surface emissivity for LW radiation.
!!\param xlon      (IMAX), longitude in radiance, ok for both 0->2pi
!!                  or -pi -> +pi ranges
!!\param xlat      (IMAX), latitude  in radiance, default to pi/2 ->
!!                  -pi/2 range, otherwise see in-line comment
!!\param slmsk     (IMAX), sea(0),land(1),ice(2) mask on fcst model grid
!!\param snowf     (IMAX), snow depth water equivalent in mm
!!\param sncovr    (IMAX), snow cover over land
!!\param zorlf     (IMAX), surface roughness in cm
!!\param tsknf     (IMAX), ground surface temperature in K
!!\param tairf     (IMAX), lowest model layer air temperature in K
!!\param hprif     (IMAX), topographic standard deviation in m
!!\param IMAX       array horizontal dimension
!!\param sfcemis  (IMAX), surface emissivity
!>\section general_setemis setemis General Algorithm
!! @{
!-----------------------------------
      subroutine setemis                                                &
     &     ( xlon,xlat,slmsk,snowf,sncovr,zorlf,tsknf,tairf,hprif,      &  !  ---  inputs:
     &       IMAX,                                                      &
     &       sfcemis                                                    &  !  ---  outputs:
     &     )

!  ===================================================================  !
!                                                                       !
!  this program computes surface emissivity for lw radiation.           !
!                                                                       !
!  usage:         call setemis                                          !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     xlon  (IMAX)  - longitude in radiance, ok for both 0->2pi or      !
!                     -pi -> +pi ranges                                 !
!     xlat  (IMAX)  - latitude  in radiance, default to pi/2 -> -pi/2   !
!                     range, otherwise see in-line comment              !
!     slmsk (IMAX)  - sea(0),land(1),ice(2) mask on fcst model grid     !
!     snowf (IMAX)  - snow depth water equivalent in mm                 !
!     sncovr(IMAX)  - ialbflg=1: snow cover over land in fraction       !
!     zorlf (IMAX)  - surface roughness in cm                           !
!     tsknf (IMAX)  - ground surface temperature in k                   !
!     tairf (IMAX)  - lowest model layer air temperature in k           !
!     hprif (IMAX)  - topographic sdv in m                              !
!     IMAX          - array horizontal dimension                        !
!                                                                       !
!  outputs:                                                             !
!     sfcemis(IMAX) - surface emissivity                                !
!                                                                       !
!  -------------------------------------------------------------------  !
!                                                                       !
!  surface type definations:                                            !
!     1. open water                   2. grass/wood/shrub land          !
!     3. tundra/bare soil             4. sandy desert                   !
!     5. rocky desert                 6. forest                         !
!     7. ice                          8. snow                           !
!                                                                       !
!  input index data lon from 0 towards east, lat from n to s            !
!                                                                       !
!  ====================    end of description    =====================  !
!
      implicit none

!  ---  inputs
      integer, intent(in) :: IMAX

      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       xlon,xlat, slmsk, snowf,sncovr, zorlf, tsknf, tairf, hprif

!  ---  outputs
      real (kind=kind_phys), dimension(:), intent(out) :: sfcemis

!  ---  locals:
      integer :: i, i1, i2, j1, j2, idx

      real (kind=kind_phys) :: dltg, hdlt, tmp1, tmp2,                  &
     &      asnow, argh, hrgh, fsno, fsno0, fsno1

!  ---  reference emiss value for diff surface emiss index
!       1-open water, 2-grass/shrub land, 3-bare soil, tundra,
!       4-sandy desert, 5-rocky desert, 6-forest, 7-ice, 8-snow

      real (kind=kind_phys) ::  emsref(8)
      data  emsref / 0.97, 0.95, 0.94, 0.90, 0.93, 0.96, 0.96, 0.99 /

!
!===> ...  begin here
!
!> -# Set sfcemis default to 1.0 or by surface type and condition.
      if ( iemslw == 0 ) then        ! sfc emiss default to 1.0

        sfcemis(:) = f_one
        return

      else                           ! emiss set by sfc type and condition

        dltg = 360.0 / float(IMXEMS)
        hdlt = 0.5 * dltg

!  --- ...  mapping input data onto model grid
!           note: this is a simple mapping method, an upgrade is needed if
!           the model grid is much corcer than the 1-deg data resolution

        lab_do_IMAX : do i = 1, IMAX

          if ( nint(slmsk(i)) == 0 ) then          ! sea point

            sfcemis(i) = emsref(1)

          else if ( nint(slmsk(i)) == 2 ) then     ! sea-ice

            sfcemis(i) = emsref(7)

          else                                     ! land

!  ---  map grid in longitude direction
            i2 = 1
            j2 = 1

            tmp1 = xlon(i) * rad2dg
            if (tmp1 < f_zero) tmp1 = tmp1 + 360.0

            lab_do_IMXEMS : do i1 = 1, IMXEMS
              tmp2 = dltg * (i1 - 1) + hdlt

              if (abs(tmp1-tmp2) <= hdlt) then
               i2 = i1
                exit lab_do_IMXEMS
              endif
            enddo  lab_do_IMXEMS

!  ---  map grid in latitude direction
            tmp1 = xlat(i) * rad2dg           ! if xlat in pi/2 -> -pi/2 range
!           tmp1 = 90.0 - xlat(i)*rad2dg      ! if xlat in 0 -> pi range

            lab_do_JMXEMS : do j1 = 1, JMXEMS
              tmp2 = 90.0 - dltg * (j1 - 1)

              if (abs(tmp1-tmp2) <= hdlt) then
                j2 = j1
                exit lab_do_JMXEMS
              endif
            enddo  lab_do_JMXEMS


            idx = max( 2, idxems(i2,j2) )
            if ( idx >= 7 ) idx = 2
            sfcemis(i) = emsref(idx)

          endif   ! end if_slmsk_block

!> -# Check for snow covered area.

          if ( ialbflg==1 .and. nint(slmsk(i))==1 ) then ! input land area snow cover

            fsno0 = sncovr(i)
            fsno1 = f_one - fsno0
            sfcemis(i) = sfcemis(i)*fsno1 + emsref(8)*fsno0

          else                                           ! compute snow cover from snow depth
            if ( snowf(i) > f_zero ) then
              asnow = 0.02*snowf(i)
              argh  = min(0.50, max(.025, 0.01*zorlf(i)))
              hrgh  = min(f_one, max(0.20, 1.0577-1.1538e-3*hprif(i) ) )
              fsno0 = asnow / (argh + asnow) * hrgh
              if (nint(slmsk(i)) == 0 .and. tsknf(i) > 271.2)           &
     &                               fsno0=f_zero
              fsno1 = f_one - fsno0
              sfcemis(i) = sfcemis(i)*fsno1 + emsref(8)*fsno0
            endif

          endif                                          ! end if_ialbflg

        enddo  lab_do_IMAX

      endif   ! end if_iemslw_block

!chk  print *,' In setemis, iemsflg, sfcemis =',iemsflg,sfcemis

!
      return
!...................................
      end subroutine setemis
!! @}
!-----------------------------------

!
!.........................................!
      end module module_radiation_surface !
!=========================================!
!> @}
