!>  \file radiation_gases.f
!!  This file contains routines that set up ozone climatological
!!  profiles and other constant gas profiles, such as co2, ch4, n2o,
!!  o2, and those of cfc gases.  All data are entered as mixing ratio
!!  by volume, except ozone which is mass mixing ratio (g/g).

!  ==========================================================  !!!!!
!              'module_radiation_gases'  description           !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!   set up ozone climatological profiles and other constant gas        !
!   profiles, such as co2, ch4, n2o, o2, and those of cfc gases.  All  !
!   data are entered as mixing ratio by volume, except ozone which is  !
!   mass mixing ratio (g/g).                                           !
!                                                                      !
!   in the module, the externally callabe subroutines are :            !
!                                                                      !
!      'gas_init'   -- initialization                                  !
!         input:                                                       !
!           ( me )                                                     !
!         output:                                                      !
!           ( none )                                                   !
!                                                                      !
!      'gas_update' -- read in data and update with time               !
!         input:                                                       !
!           ( iyear, imon, iday, ihour, loz1st, ldoco2, me )           !
!         output:                                                      !
!           ( none )                                                   !
!                                                                      !
!      'getozn'     -- setup climatological ozone profile              !
!         input:                                                       !
!           ( prslk,xlat,                                              !
!             IMAX, LM )                                               !
!         output:                                                      !
!           ( o3mmr )                                                  !
!                                                                      !
!      'getgases'   -- setup constant gas profiles for LW and SW       !
!         input:                                                       !
!           ( plvl, xlon, xlat,                                        !
!             IMAX, LMAX )                                             !
!         output:                                                      !
!           ( gasdat )                                                 !
!                                                                      !
!   external modules referenced:                                       !
!       'module machine'                    in 'machine.f'             !
!       'module funcphys'                   in 'funcphys.f'            !
!       'module physcons'                   in 'physcons.f             !
!       'module module_iounitdef'           in 'iounitdef.f'           !
!                                                                      !
!   unit used for radiative active gases:                              !
!      ozone : mass mixing ratio                     (g/g)             !
!      co2   : volume mixing ratio                   (p/p)             !
!      n2o   : volume mixing ratio                   (p/p)             !
!      ch4   : volume mixing ratio                   (p/p)             !
!      o2    : volume mixing ratio                   (p/p)             !
!      co    : volume mixing ratio                   (p/p)             !
!      cfc11 : volume mixing ratio                   (p/p)             !
!      cfc12 : volume mixing ratio                   (p/p)             !
!      cfc22 : volume mixing ratio                   (p/p)             !
!      ccl4  : volume mixing ratio                   (p/p)             !
!      cfc113: volume mixing ratio                   (p/p)             !
!                                                                      !
!                                                                      !
!   program history:                                                   !
!     may 2003 - y-t hou     create rad_module.f that collectively     !
!                  combines several radiation computation supporting   !
!                  programs into fortran 90 module structure (gases    !
!                  and aerosols, etc.)                                 !
!     apr 2004 - y-t hou     modified to add astronomy and surface     !
!                  module components.                                  !
!     feb 2005 - y-t hou     rewrite the component modules into        !
!                  separate individule modules for thier corresponding !
!                  tasks. here as radiation_gases.f                    !
!     mar 2006 - y-t hou     add initialization subroutine to co2 and  !
!                  other gases. historical 2-d co2 data are added.     !
!     sep 2008 - y-t hou     add parameter ictm to control the input   !
!                  data time at the model initial condition.           !
!     oct 2008 - y-t hou     modify the initialization code to add the !
!                  option of superimposing climatology seasonal cycle  !
!                  to the initial condition data (currently co2 only)  !
!     nov 2008 - y-t hou     fix bugs in superimposing climatology     !
!                  seasonal cycle calculations                         !
!     aug 2011 - y-t hou     fix a bug in subr getgases doing vertical !
!                  co2 mapping. (for iflip=0 case, not affact opr).    !
!     aug 2012 - y-t hou     modified subr getozn.  moved the if-first !
!                  block to subr gas_init to ensure threading safe in  !
!                  climatology ozone applications. (not affect gfs)    !
!                  also changed the initialization subr into two parts:!
!                  'gas_init' is called at the start of run to set up  !
!                  module parameters; and 'gas_update' is called within!
!                  the time loop to check and update data sets. defined!
!                  the climatology ozone parameters k1oz,k2oz,facoz as !
!                  module variables and are set in subr 'gas_update'   !
!     nov 2012 - y-t hou     modified control parameters thru module   !
!                  'physparam'.                                        !
!     jan 2013 - z. janjic/y. hou   modified ilon (longitude index)    !
!                  computing formula in subroutine getgases to work    !
!                  properly for models with either of  0->360 or       !
!                  -180->180 zonal grid directions.                    !
!                                                                      !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!


!> \defgroup module_radiation_gases RRTMG Gases Module
!> This module sets up ozone climatological profiles and other constant
!! gas profiles, such as co2, ch4, n2o, o2, and those of cfc gases. All
!! data are entered as mixing ratio by volume, except ozone which is
!! mass mixing ratio (g/g).
!!\version NCEP-Radiation_gases     v5.1  Nov 2012
!========================================!
      module module_radiation_gases      !
!........................................!
!
      use physparam,         only : ico2flg, ictmflg, ioznflg, ivflip,  &
     &                              co2dat_file, co2gbl_file,           &
     &                              co2usr_file, co2cyc_file,           &
     &                              kind_phys, kind_io4
      use funcphys,          only : fpkapx
      use physcons,          only : con_pi
      use ozne_def,          only : JMR => latsozc, LOZ => levozc,      &
     &                              blte => blatc, dlte=> dphiozc,      &
     &                              timeozc => timeozc
      use module_iounitdef,  only : NIO3CLM, NICO2CN
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGGAS='NCEP-Radiation_gases     v5.1  Nov 2012 '
!    &   VTAGGAS='NCEP-Radiation_gases     v5.0  Aug 2012 '

      integer, parameter, public :: NF_VGAS = 10   !< number of gas species
      integer, parameter         :: IMXCO2  = 24   !< input CO2 data longitude points
      integer, parameter         :: JMXCO2  = 12   !< input CO2 data latitude points
      integer, parameter         :: MINYEAR = 1957 !< earlist year 2D CO2 data available

      real (kind=kind_phys), parameter :: resco2=15.0            !< horizontal resolution in degree
      real (kind=kind_phys), parameter :: raddeg=180.0/con_pi    !< rad->deg conversion
      real (kind=kind_phys), parameter :: prsco2=788.0           !< pressure limitation for 2D CO2 (mb)
      real (kind=kind_phys), parameter :: hfpi  =0.5*con_pi      !< half of pi

      real (kind=kind_phys), parameter :: co2vmr_def = 350.0e-6  !< parameter constant for CO2 volume mixing ratio
      real (kind=kind_phys), parameter :: n2ovmr_def = 0.31e-6   !< parameter constant for N2O volume mixing ratio
      real (kind=kind_phys), parameter :: ch4vmr_def = 1.50e-6   !< parameter constant for CH4 volume mixing ratio
      real (kind=kind_phys), parameter :: o2vmr_def  = 0.209     !< parameter constant for O2  volume mixing ratio
      real (kind=kind_phys), parameter :: covmr_def  = 1.50e-8   !< parameter constant for CO  colume mixing ratio
! aer 2003 value
      real (kind=kind_phys), parameter :: f11vmr_def = 3.520e-10
! aer 2003 value
      real (kind=kind_phys), parameter :: f12vmr_def = 6.358e-10
! aer 2003 value
      real (kind=kind_phys), parameter :: f22vmr_def = 1.500e-10
! aer 2003 value
      real (kind=kind_phys), parameter :: cl4vmr_def = 1.397e-10
! gfdl 1999 value
      real (kind=kind_phys), parameter :: f113vmr_def= 8.2000e-11

!  ---  ozone seasonal climatology parameters defined in module ozne_def
!   - 4x5 ozone data parameter
!     integer, parameter :: JMR=45, LOZ=17
!     real (kind=kind_phys), parameter :: blte=-86.0, dlte=4.0
!   - geos ozone data
!     integer, parameter :: JMR=18, LOZ=17
!     real (kind=kind_phys), parameter :: blte=-85.0, dlte=10.0

!  ---  module variables to be set in subroutin gas_init and/or gas_update

! variables for climatology ozone (ioznflg = 0)

      real (kind=kind_phys), allocatable :: pkstr(:), o3r(:,:,:)
      integer :: k1oz = 0,  k2oz = 0
      real (kind=kind_phys) :: facoz = 0.0

!  arrays for co2 2-d monthly data and global mean values from observed data

      real (kind=kind_phys), allocatable :: co2vmr_sav(:,:,:)
      real (kind=kind_phys), allocatable :: co2cyc_sav(:,:,:)

      real (kind=kind_phys) :: co2_glb = co2vmr_def
      real (kind=kind_phys) :: gco2cyc(12)
      data gco2cyc(:) / 12*0.0 /

      integer :: kyrsav  = 0
      integer :: kmonsav = 1

!  ---  public interfaces

      public  gas_init, gas_update, getgases, getozn


! =================
      contains
! =================

!> \ingroup module_radiation_gases
!> This subroutine sets up ozone, co2, etc. parameters. If climatology
!! ozone then read in monthly ozone data.
!!\param me         print message control flag
!>\section gas_init_gen gas_init General Algorithm
!! @{
!-----------------------------------
      subroutine gas_init                                               &
     &     ( me )!  ---  inputs:
!  ---  outputs: ( none )

!  ===================================================================  !
!                                                                       !
!  gas_init sets up ozone, co2, etc. parameters.  if climatology ozone  !
!  then read in monthly ozone data.                                     !
!                                                                       !
!  inputs:                                               dimemsion      !
!     me      - print message control flag                  1           !
!                                                                       !
!  outputs: (to the module variables)                                   !
!    ( none )                                                           !
!                                                                       !
!  external module variables:  (in physparam)                           !
!     ico2flg    - co2 data source control flag                         !
!                   =0: use prescribed co2 global mean value            !
!                   =1: use input global mean co2 value (co2_glb)       !
!                   =2: use input 2-d monthly co2 value (co2vmr_sav)    !
!     ictmflg    - =yyyy#, data ic time/date control flag               !
!                  =   -2: same as 0, but superimpose seasonal cycle    !
!                          from climatology data set.                   !
!                  =   -1: use user provided external data for the fcst !
!                          time, no extrapolation.                      !
!                  =    0: use data at initial cond time, if not existed!
!                          then use latest, without extrapolation.      !
!                  =    1: use data at the forecast time, if not existed!
!                          then use latest and extrapolate to fcst time.!
!                  =yyyy0: use yyyy data for the forecast time, no      !
!                          further data extrapolation.                  !
!                  =yyyy1: use yyyy data for the fcst. if needed, do    !
!                          extrapolation to match the fcst time.        !
!     ioznflg    - ozone data control flag                              !
!                   =0: use climatological ozone profile                !
!                   >0: use interactive ozone profile                   !
!     ivflip     - vertical profile indexing flag                       !
!     co2usr_file- external co2 user defined data table                 !
!     co2cyc_file- external co2 climotology monthly cycle data table    !
!                                                                       !
!  internal module variables:                                           !
!     pkstr, o3r - arrays for climatology ozone data                    !
!                                                                       !
!  usage:    call gas_init                                              !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: me

!  ---  output: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(IMXCO2,JMXCO2) :: co2dat
      real (kind=kind_phys) :: co2g1, co2g2
      real (kind=kind_phys) :: pstr(LOZ)
      real (kind=kind_io4)  :: o3clim4(JMR,LOZ,12), pstr4(LOZ)

      integer    :: imond(12), ilat(JMR,12)
      integer    :: i, j, k, iyr, imo
      logical    :: file_exist, lextpl
      character  :: cline*100, cform*8
      data  cform  / '(24f7.2)' /       !! data format in IMXCO2*f7.2
!
!===>  ...  begin here
!
      if ( me == 0 ) print *, VTAGGAS    ! print out version tag

      kyrsav  = 0
      kmonsav = 1

!  --- ...  climatology ozone data section

      if ( ioznflg > 0 ) then
        if ( me == 0 ) then
          print *,' - Using interactive ozone distribution'
        endif
      else
        if ( timeozc /= 12 ) then
          print *,' - Using climatology ozone distribution'
          print *,' timeozc=',timeozc, ' is not monthly mean',          &
     &            ' - job aborting in subroutin gas_init!!!'
          stop
        endif

        allocate (pkstr(LOZ), o3r(JMR,LOZ,12))
        rewind NIO3CLM

        if ( LOZ == 17 ) then       ! For the operational ozone climatology
          do k = 1, LOZ
            read (NIO3CLM,15) pstr4(k)
   15       format(f10.3)
          enddo

          do imo = 1, 12
            do j = 1, JMR
              read (NIO3CLM,16) imond(imo), ilat(j,imo),                &
     &                          (o3clim4(j,k,imo),k=1,10)
   16         format(i2,i4,10f6.2)
              read (NIO3CLM,20) (o3clim4(j,k,imo),k=11,LOZ)
   20         format(6x,10f6.2)
            enddo
          enddo
        else                      ! For newer ozone climatology
          read (NIO3CLM)
          do k = 1, LOZ
            read (NIO3CLM) pstr4(k)
          enddo

          do imo = 1, 12
            do k = 1, LOZ
              read (NIO3CLM) (o3clim4(j,k,imo),j=1,JMR)
            enddo
          enddo
        endif   ! end if_LOZ_block
!
        do imo = 1, 12
          do k = 1, LOZ
            do j = 1, JMR
              o3r(j,k,imo) = o3clim4(j,k,imo) * 1.655e-6
            enddo
          enddo
        enddo

        do k = 1, LOZ
          pstr(k) = pstr4(k)
        enddo

        if ( me == 0 ) then
          print *,' - Using climatology ozone distribution'
          print *,'   Found ozone data for levels pstr=',               &
     &            (pstr(k),k=1,LOZ)
!         print *,' O3=',(o3r(15,k,1),k=1,LOZ)
        endif

        do k = 1, LOZ
          pkstr(k) = fpkapx(pstr(k)*100.0)
        enddo
      endif   ! end if_ioznflg_block

!  --- ...  co2 data section

      co2_glb = co2vmr_def

      lab_ico2 : if ( ico2flg == 0 ) then

        if ( me == 0 ) then
          print *,' - Using prescribed co2 global mean value=',         &
     &              co2vmr_def
        endif

      else  lab_ico2

        lab_ictm : if ( ictmflg == -1 ) then      ! input user provided data

          inquire (file=co2usr_file, exist=file_exist)
          if ( .not. file_exist ) then
            print *,'   Can not find user CO2 data file: ',co2usr_file, &
     &              ' - Stopped in subroutine gas_init !!'
            stop
          else
            close (NICO2CN)
            open(NICO2CN,file=co2usr_file,form='formatted',status='old')
            rewind NICO2CN
            read (NICO2CN, 25) iyr, cline, co2g1, co2g2
  25        format(i4,a94,f7.2,16x,f5.2)
            co2_glb = co2g1 * 1.0e-6

            if ( ico2flg == 1 ) then
              if ( me == 0 ) then
                print *,' - Using co2 global annual mean value from',   &
     &                  ' user provided data set:',co2usr_file
                print *, iyr,cline(1:94),co2g1,'  GROWTH RATE =', co2g2
              endif
            elseif ( ico2flg == 2 ) then
              allocate ( co2vmr_sav(IMXCO2,JMXCO2,12) )

              do imo = 1, 12
                read (NICO2CN,cform) co2dat
!check          print cform, co2dat

                do j = 1, JMXCO2
                  do i = 1, IMXCO2
                    co2vmr_sav(i,j,imo) = co2dat(i,j) * 1.0e-6
                  enddo
                enddo
              enddo

              if ( me == 0 ) then
                print *,' - Using co2 monthly 2-d data from user',      &
     &                ' provided data set:',co2usr_file
                print *, iyr,cline(1:94),co2g1,'  GROWTH RATE =', co2g2

                print *,' CHECK: Sample of selected months of CO2 data'
                do imo = 1, 12, 3
                  print *,'        Month =',imo
                  print *, (co2vmr_sav(1,j,imo),j=1,jmxco2)
                enddo
              endif
            else
              print *,' ICO2=',ico2flg,' is not a valid selection',     &
     &                ' - Stoped in subroutine gas_init!!!'
              stop
            endif    ! endif_ico2flg_block

            close (NICO2CN)
          endif    ! endif_file_exist_block

        else   lab_ictm                           ! input from observed data

          if ( ico2flg == 1 ) then
            if ( me == 0 ) then
              print *,' - Using observed co2 global annual mean value'
            endiF
          elseif ( ico2flg == 2 ) then
            allocate ( co2vmr_sav(IMXCO2,JMXCO2,12) )

            if ( me == 0 ) then
              print *,' - Using observed co2 monthly 2-d data'
            endif
          else
            print *,' ICO2=',ico2flg,' is not a valid selection',       &
     &              ' - Stoped in subroutine gas_init!!!'
            stop
          endif

          if ( ictmflg == -2 ) then
            inquire (file=co2cyc_file, exist=file_exist)
            if ( .not. file_exist ) then
              if ( me == 0 ) then
                print *,'   Can not find seasonal cycle CO2 data: ',    &
     &               co2cyc_file,' - Stopped in subroutine gas_init !!'
              endif
              stop
            else
              allocate( co2cyc_sav(IMXCO2,JMXCO2,12) )

!  --- ...  read in co2 2-d seasonal cycle data
              close (NICO2CN)
              open (NICO2CN,file=co2cyc_file,form='formatted',          &
     &              status='old')
              rewind NICO2CN
              read (NICO2CN, 35) cline, co2g1, co2g2
  35          format(a98,f7.2,16x,f5.2)
              read (NICO2CN,cform) co2dat        ! skip annual mean part

              if ( me == 0 ) then
                print *,' - Superimpose seasonal cycle to mean CO2 data'
                print *,'   Opened CO2 climatology seasonal cycle data',&
     &                  ' file: ',co2cyc_file
!check          print *, cline(1:98), co2g1, co2g2
              endif

              do imo = 1, 12
                read (NICO2CN,45) cline, gco2cyc(imo)
  45            format(a58,f7.2)
!check          print *, cline(1:58),gco2cyc(imo)
                gco2cyc(imo) = gco2cyc(imo) * 1.0e-6

                read (NICO2CN,cform) co2dat
!check          print cform, co2dat
                do j = 1, JMXCO2
                  do i = 1, IMXCO2
                    co2cyc_sav(i,j,imo) = co2dat(i,j) * 1.0e-6
                  enddo
                enddo
              enddo

              close (NICO2CN)
            endif   ! endif_file_exist_block
          endif

        endif   lab_ictm
      endif   lab_ico2

      return
!
!...................................
      end subroutine gas_init
!! @}
!-----------------------------------

!> \ingroup module_radiation_gases
!> This subroutine reads in 2-d monthly co2 data set for a specified
!! year. Data are in a 15 degree lat/lon horizontal resolution.
!!\param iyear      year of the requested data for fcst
!!\param imon       month of the year
!!\param iday       day of the month
!!\param ihour      hour of the day
!!\param loz1st     clim ozone 1st time update control flag
!!\param ldoco2     co2 update control flag
!!\param me         print message control flag
!>\section gen_gas_update gas_update General Algorithm
!! @{
!-----------------------------------
      subroutine gas_update                                             &
     &     ( iyear, imon, iday, ihour, loz1st, ldoco2, me )!  ---  inputs
!  ---  outputs: ( none )

!  ===================================================================  !
!                                                                       !
!  gas_update reads in 2-d monthly co2 data set for a specified year.   !
!  data are in a 15 degree lat/lon horizontal resolution.               !
!                                                                       !
!  inputs:                                               dimemsion      !
!     iyear   - year of the requested data for fcst         1           !
!     imon    - month of the year                           1           !
!     iday    - day of the month                            1           !
!     ihour   - hour of the day                             1           !
!     loz1st  - clim ozone 1st time update control flag     1           !
!     ldoco2  - co2 update control flag                     1           !
!     me      - print message control flag                  1           !
!                                                                       !
!  outputs: (to the module variables)                                   !
!    ( none )                                                           !
!                                                                       !
!  external module variables:  (in physparam)                           !
!     ico2flg    - co2 data source control flag                         !
!                   =0: use prescribed co2 global mean value            !
!                   =1: use input global mean co2 value (co2_glb)       !
!                   =2: use input 2-d monthly co2 value (co2vmr_sav)    !
!     ictmflg    - =yyyy#, data ic time/date control flag               !
!                  =   -2: same as 0, but superimpose seasonal cycle    !
!                          from climatology data set.                   !
!                  =   -1: use user provided external data for the fcst !
!                          time, no extrapolation.                      !
!                  =    0: use data at initial cond time, if not existed!
!                          then use latest, without extrapolation.      !
!                  =    1: use data at the forecast time, if not existed!
!                          then use latest and extrapolate to fcst time.!
!                  =yyyy0: use yyyy data for the forecast time, no      !
!                          further data extrapolation.                  !
!                  =yyyy1: use yyyy data for the fcst. if needed, do    !
!                          extrapolation to match the fcst time.        !
!     ioznflg    - ozone data control flag                              !
!                   =0: use climatological ozone profile                !
!                   >0: use interactive ozone profile                   !
!     ivflip     - vertical profile indexing flag                       !
!     co2dat_file- external co2 2d monthly obsv data table              !
!     co2gbl_file- external co2 global annual mean data table           !
!                                                                       !
!  internal module variables:                                           !
!     co2vmr_sav - monthly co2 volume mixing ratio     IMXCO2*JMXCO2*12 !
!     co2cyc_sav - monthly cycle co2 vol mixing ratio  IMXCO2*JMXCO2*12 !
!     co2_glb    - global annual mean co2 mixing ratio                  !
!     gco2cyc    - global monthly mean co2 variation       12           !
!     k1oz,k2oz,facoz                                                   !
!                - climatology ozone parameters             1           !
!                                                                       !
!  usage:    call gas_update                                            !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: iyear, imon, iday, ihour, me

      logical, intent(in) :: loz1st, ldoco2

!  ---  output: ( none )

!  ---  locals:
      real (kind=kind_phys), dimension(IMXCO2,JMXCO2) :: co2dat, co2ann
      real (kind=kind_phys) :: co2g1, co2g2, rate

      integer    :: i, id, j, l, iyr, imo, iyr1, iyr2, jyr, idyr
      integer, save :: mdays(13), midmon=15, midm=15, midp=45
!  ---  number of days in a month
      data mdays / 31,28,31,30,31,30,31,31,30,31,30,31,30 /

      logical    :: file_exist, lextpl, change
      character  :: cline*100, cform*8, cfile1*26
      data  cform  / '(24f7.2)' /       !! data format in IMXCO2*f7.2
!
!===>  ...  begin here
!
!> - Ozone data section

      if ( ioznflg == 0 ) then
        midmon = mdays(imon)/2 + 1
        change = loz1st .or. ( (iday==midmon) .and. (ihour==0) )
!
        if ( change ) then
          if ( iday < midmon ) then
            k1oz = mod(imon+10, 12) + 1
            midm = mdays(k1oz)/2 + 1
            k2oz = imon
            midp = mdays(k1oz) + midmon
          else
            k1oz = imon
            midm = midmon
            k2oz = mod(imon, 12) + 1
            midp = mdays(k2oz)/2 + 1 + mdays(k1oz)
          endif
        endif
!
        if (iday < midmon) then
         id = iday + mdays(k1oz)
        else
         id = iday
        endif

        facoz = float(id - midm) / float(midp - midm)
      endif

!> - co2 data section

      if ( ico2flg == 0 ) return    ! use prescribed global mean co2 data
      if ( ictmflg ==-1 ) return    ! use user provided co2 data
      if ( .not. ldoco2 ) return    ! no need to update co2 data

      if ( ictmflg < 0 ) then        ! use user provided external data
        lextpl = .false.                   ! no time extrapolation
        idyr   = iyear                     ! use the model year
      else                           ! use historically observed data
        lextpl = ( mod(ictmflg,10) == 1 )  ! flag for data extrapolation
        idyr   = ictmflg / 10              ! year of data source used
        if ( idyr == 0 ) idyr = iyear      ! not specified, use model year
      endif

!  --- ...  auto select co2 2-d data table for required year

      kmonsav = imon
      if ( kyrsav == iyear ) return
      kyrsav = iyear
      iyr    = iyear

!  --- ...  for data earlier than MINYEAR (1957), the data are in
!           the form of semi-yearly global mean values.  otherwise,
!           data are monthly mean in horizontal 2-d map.

      Lab_if_idyr : if ( idyr < MINYEAR .and. ictmflg > 0 ) then

        if ( me == 0 ) then
          print *,'   Requested CO2 data year',iyear,' earlier than',   &
     &            MINYEAR
          print *,'   Which is the earliest monthly observation',       &
     &            ' data available.'
          print *,'   Thus, historical global mean data is used'
        endif

!  --- ... check to see if requested co2 data file existed

        inquire (file=co2gbl_file, exist=file_exist)
        if ( .not. file_exist ) then
          print *,'   Requested co2 data file "',co2gbl_file,           &
     &            '" not found - Stopped in subroutine gas_update!!'
          stop
        else
          close(NICO2CN)
          open (NICO2CN,file=co2gbl_file,form='formatted',status='old')
          rewind NICO2CN

          read (NICO2CN, 24) iyr1, iyr2, cline
  24      format(i4,4x,i4,a48)

          if ( me == 0 ) then
            print *,'   Opened co2 data file: ',co2gbl_file
!check      print *, iyr1, iyr2, cline(1:48)
          endif

          if ( idyr < iyr1 ) then
            iyr = iyr1
!check      if ( me == 0 ) then
!             print *,'   Using earlist available co2 data, year=',iyr1
!check      endif
          endif

          i = iyr2
          Lab_dowhile1 : do while ( i >= iyr1 )
!           read (NICO2CN,26) jyr, co2g1, co2g2
! 26        format(i4,4x,2f7.2)
            read (NICO2CN, *) jyr, co2g1, co2g2

            if ( i == iyr .and. iyr == jyr ) then
              co2_glb = (co2g1+co2g2) * 0.5e-6
              if ( ico2flg == 2 ) then
                do j = 1, JMXCO2
                  do i = 1, IMXCO2
                    co2vmr_sav(i,j,1:6)  = co2g1 * 1.0e-6
                    co2vmr_sav(i,j,7:12) = co2g2 * 1.0e-6
                  enddo
                enddo
              endif

              if ( me == 0 ) print *,'   Co2 data for year',iyear,      &
     &                               co2_glb
              exit Lab_dowhile1
            else
!check        if ( me == 0 ) print *,'   Skip co2 data for year',i
              i = i - 1
            endif
          enddo  Lab_dowhile1

          close ( NICO2CN )
        endif   ! end if_file_exist_block

      else  Lab_if_idyr

!  --- ...  set up input data file name

        cfile1 = co2dat_file
        write(cfile1(19:22),34) idyr
  34    format(i4.4)

!  --- ... check to see if requested co2 data file existed

        inquire (file=cfile1, exist=file_exist)
        if ( .not. file_exist ) then

          Lab_if_ictm : if ( ictmflg  > 10 ) then    ! specified year of data not found
            if ( me == 0 ) then
              print *,'   Specified co2 data for year',idyr,            &
     &               ' not found !!  Need to change namelist ICTM !!'
              print *,'   *** Stopped in subroutine gas_update !!'
            endif
            stop
          else Lab_if_ictm                        ! looking for latest available data
            if ( me == 0 ) then
              print *,'   Requested co2 data for year',idyr,            &
     &              ' not found, check for other available data set'
            endif

            Lab_dowhile2 : do while ( iyr >= MINYEAR )
              iyr = iyr - 1
              write(cfile1(19:22),34) iyr

              inquire (file=cfile1, exist=file_exist)
              if ( me == 0 ) then
                print *,' Looking for CO2 file ',cfile1
              endif

              if ( file_exist ) then
                exit Lab_dowhile2
              endif
            enddo   Lab_dowhile2

            if ( .not. file_exist ) then
              if ( me == 0 ) then
                print *,'   Can not find co2 data source file'
                print *,'   *** Stopped in subroutine gas_update !!'
              endif
              stop
            endif
          endif  Lab_if_ictm
        endif   ! end if_file_exist_block

!  --- ...  read in co2 2-d data for the requested month

        close(NICO2CN)
        open (NICO2CN,file=cfile1,form='formatted',status='old')
        rewind NICO2CN
        read (NICO2CN, 36) iyr, cline, co2g1, co2g2
  36    format(i4,a94,f7.2,16x,f5.2)

        if ( me == 0 ) then
          print *,'   Opened co2 data file: ',cfile1
          print *, iyr, cline(1:94), co2g1,'  GROWTH RATE =', co2g2
        endif

!  --- ...  add growth rate if needed
        if ( lextpl ) then
!         rate = co2g2 * (iyear - iyr)   ! rate from early year
!         rate = 1.60  * (iyear - iyr)   ! avg rate over long period
          rate = 2.00  * (iyear - iyr)   ! avg rate for recent period
        else
          rate = 0.0
        endif

        co2_glb = (co2g1 + rate) * 1.0e-6
        if ( me == 0 ) then
          print *,'   Global annual mean CO2 data for year',            &
     &              iyear, co2_glb
        endif

        if ( ictmflg == -2 ) then     ! need to calc ic time annual mean first

          if ( ico2flg == 1 ) then
            if ( me==0 ) then
              print *,' CHECK: Monthly deviations of climatology ',     &
     &                'to be superimposed on global annual mean'
              print *, gco2cyc
            endif
          elseif ( ico2flg == 2 ) then
            co2ann(:,:) = 0.0

            do imo = 1, 12
              read (NICO2CN,cform) co2dat
!check        print cform, co2dat

              do j = 1, JMXCO2
                do i = 1, IMXCO2
                  co2ann(i,j) = co2ann(i,j) + co2dat(i,j)
                enddo
              enddo
            enddo

            do j = 1, JMXCO2
              do i = 1, IMXCO2
                co2ann(i,j) = co2ann(i,j) * 1.0e-6 / float(12)
              enddo
            enddo

            do imo = 1, 12
              do j = 1, JMXCO2
                do i = 1, IMXCO2
                  co2vmr_sav(i,j,imo) = co2ann(i,j)+co2cyc_sav(i,j,imo)
                enddo
              enddo
            enddo

            if ( me==0 ) then
              print *,' CHECK: Sample of 2-d annual mean of CO2 ',      &
     &                'data used for year:',iyear
              print *, co2ann(1,:)
              print *,' CHECK: AFTER adding seasonal cycle, Sample ',   &
     &                'of selected months of CO2 data for year:',iyear
              do imo = 1, 12, 3
                print *,'        Month =',imo
                print *, co2vmr_sav(1,:,imo)
              enddo
            endif
          endif   ! endif_icl2flg_block

        else                  ! no need to calc ic time annual mean first

          if ( ico2flg == 2 ) then      ! directly save monthly data
            do imo = 1, 12
              read (NICO2CN,cform) co2dat
!check        print cform, co2dat

              do j = 1, JMXCO2
                do i = 1, IMXCO2
                  co2vmr_sav(i,j,imo) = (co2dat(i,j) + rate) * 1.0e-6
                enddo
              enddo
            enddo

            if ( me == 0 ) then
              print *,' CHECK: Sample of selected months of CO2 ',      &
     &                'data used for year:',iyear
              do imo = 1, 12, 3
                print *,'        Month =',imo
                print *, co2vmr_sav(1,:,imo)
              enddo
            endif
          endif   ! endif_ico2flg_block

          do imo = 1, 12
            gco2cyc(imo) = 0.0
          enddo
        endif   ! endif_ictmflg_block
        close ( NICO2CN )

      endif  Lab_if_idyr

      return
!
!...................................
      end subroutine gas_update
!-----------------------------------
!! @}

!> \ingroup module_radiation_gases
!> This subroutine sets up global distribution of radiation absorbing
!! gases in volume mixing ratio. Currently only co2 has the options
!! from observed values, all other gases are asigned to the
!! climatological values.
!!\param plvl       (IMAX,LMAX+1), pressure at model layer interfaces (mb)
!!\param xlon       (IMAX), grid longitude in radians, ok both 0->2pi
!!                  or -pi -> +pi arrangements
!!\param xlat       (IMAX), grid latitude in radians, default range to
!!                  pi/2 -> -pi/2, otherwise see in-line comment
!!\param IMAX, LMAX      horizontal/vertical dimensions for output data
!!\param gasdat     (IMAX,LMAX,NF_VGAS) - gases volume mixing ratioes
!!\n                    (:,:,1)           - co2
!!\n                    (:,:,2)           - n2o
!!\n                    (:,:,3)           - ch4
!!\n                    (:,:,4)           - o2
!!\n                    (:,:,5)           - co
!!\n                    (:,:,6)           - cfc11
!!\n                    (:,:,7)           - cfc12
!!\n                    (:,:,8)           - cfc22
!!\n                    (:,:,9)           - ccl4
!!\n                    (:,:,10)          - cfc113
!>\section getgases_gen getgases General Algorithm
!! @{
!-----------------------------------
      subroutine getgases                                               &
     &     ( plvl, xlon, xlat,                                          & ! ---  inputs
     &       IMAX, LMAX,                                                &
     &       gasdat                                                     & ! ---  outputs
     &      )

!  ===================================================================  !
!                                                                       !
!  getgases set up global distribution of radiation absorbing  gases    !
!  in volume mixing ratio.  currently only co2 has the options from     !
!  observed values, all other gases are asigned to the climatological   !
!  values.                                                              !
!                                                                       !
!  inputs:                                                              !
!     plvl(IMAX,LMAX+1)- pressure at model layer interfaces (mb)        !
!     xlon(IMAX)       - grid longitude in radians, ok both 0->2pi or   !
!                        -pi -> +pi arrangements                        !
!     xlat(IMAX)       - grid latitude in radians, default range to     !
!                        pi/2 -> -pi/2, otherwise see in-line comment   !
!     IMAX, LMAX       - horiz, vert dimensions for output data         !
!                                                                       !
!  outputs:                                                             !
!     gasdat(IMAX,LMAX,NF_VGAS) - gases volume mixing ratioes           !
!               (:,:,1)           - co2                                 !
!               (:,:,2)           - n2o                                 !
!               (:,:,3)           - ch4                                 !
!               (:,:,4)           - o2                                  !
!               (:,:,5)           - co                                  !
!               (:,:,6)           - cfc11                               !
!               (:,:,7)           - cfc12                               !
!               (:,:,8)           - cfc22                               !
!               (:,:,9)           - ccl4                                !
!               (:,:,10)          - cfc113                              !
!                                                                       !
!  external module variables:  (in physparam)                           !
!     ico2flg    - co2 data source control flag                         !
!                   =0: use prescribed co2 global mean value            !
!                   =1: use input global mean co2 value (co2_glb)       !
!                   =2: use input 2-d monthly co2 value (co2vmr_sav)    !
!     ivflip     - vertical profile indexing flag                       !
!                                                                       !
!  internal module variables used:                                      !
!     co2vmr_sav - saved monthly co2 concentration from sub gas_update  !
!     co2_glb    - saved global annual mean co2 value from  gas_update  !
!     gco2cyc    - saved global seasonal variation of co2 climatology   !
!                  in 12-month form                                     !
!  ** note: for lower atmos co2vmr_sav may have clim monthly deviations !
!           superimposed on init-cond co2 value, while co2_glb only     !
!           contains the global mean value, thus needs to add the       !
!           monthly dglobal mean deviation gco2cyc at upper atmos. for  !
!           ictmflg/=-2, this value will be zero.                       !
!                                                                       !
!  usage:    call getgases                                              !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  input:
      integer,  intent(in)  :: IMAX, LMAX
      real (kind=kind_phys), intent(in) :: plvl(:,:), xlon(:), xlat(:)

!  ---  output:
      real (kind=kind_phys), intent(out) :: gasdat(:,:,:)

!  ---  local:
      integer :: i, k, ilat, ilon

      real (kind=kind_phys) :: xlon1, xlat1, tmp

!===>  ...  begin here

!  --- ...  assign default values

      do k = 1, LMAX
      do i = 1, IMAX
        gasdat(i,k,1) = co2vmr_def
        gasdat(i,k,2) = n2ovmr_def
        gasdat(i,k,3) = ch4vmr_def
        gasdat(i,k,4) = o2vmr_def
        gasdat(i,k,5) = covmr_def
        gasdat(i,k,6) = f11vmr_def
        gasdat(i,k,7) = f12vmr_def
        gasdat(i,k,8) = f22vmr_def
        gasdat(i,k,9) = cl4vmr_def
        gasdat(i,k,10)= f113vmr_def
      enddo
      enddo

!  --- ...  co2 section

      if ( ico2flg == 1 ) then
!  ---  use obs co2 global annual mean value only

        do k = 1, LMAX
          do i = 1, IMAX
            gasdat(i,k,1) = co2_glb + gco2cyc(kmonsav)
          enddo
        enddo

      elseif ( ico2flg == 2 ) then
!  ---  use obs co2 monthly data with 2-d variation at lower atmos
!       otherwise use global mean value

        tmp = raddeg / resco2
        do i = 1, IMAX
          xlon1 = xlon(i)
          if ( xlon1 < 0.0 ) xlon1 = xlon1 + con_pi  ! if xlon in -pi->pi, convert to 0->2pi
          xlat1 = hfpi - xlat(i)                     ! if xlat in pi/2 -> -pi/2 range
!note     xlat1 = xlat(i)                            ! if xlat in 0 -> pi range

          ilon = min( IMXCO2, int( xlon1*tmp + 1 ))
          ilat = min( JMXCO2, int( xlat1*tmp + 1 ))

          if ( ivflip == 0 ) then         ! index from toa to sfc
            do k = 1, LMAX
              if ( plvl(i,k) >= prsco2 ) then
                gasdat(i,k,1) = co2vmr_sav(ilon,ilat,kmonsav)
              else
                gasdat(i,k,1) = co2_glb + gco2cyc(kmonsav)
              endif
            enddo
          else                            ! index from sfc to toa
            do k = 1, LMAX
              if ( plvl(i,k+1) >= prsco2 ) then
                gasdat(i,k,1) = co2vmr_sav(ilon,ilat,kmonsav)
              else
                gasdat(i,k,1) = co2_glb + gco2cyc(kmonsav)
              endif
            enddo
          endif
        enddo
      endif

!
      return
!...................................
      end subroutine getgases
!! @}
!-----------------------------------

!> \ingroup module_radiation_gases
!> This subroutine sets up climatological ozone profile for radiation
!! calculation. This code is originally written by Shrinivas Moorthi.
!!\param prslk       (IMAX,LM), exner function = \f$(p/p0)^{rocp}\f$
!!\param xlat        (IMAX), latitude in radians, default to pi/2 ->
!!                    -pi/2 range, otherwise see in-line comment
!!\param IMAX, LM    horizontal and vertical dimensions
!!\param o3mmr       (IMAX,LM), output ozone profile in mass mixing
!!                   ratio (g/g)
!>\section getozn_gen getozn General Algorithm
!! @{
!-----------------------------------
      subroutine getozn                                                 &
     &     ( prslk,xlat,                                                &                    !  ---  inputs
     &       IMAX, LM,                                                  &
     &       o3mmr                                                      &                    !  ---  outputs
     &     )

!  ===================================================================  !
!                                                                       !
!  getozn sets up climatological ozone profile for radiation calculation!
!                                                                       !
!  this code is originally written By Shrinivas Moorthi                 !
!                                                                       !
!  inputs:                                                              !
!     prslk (IMAX,LM)  - exner function = (p/p0)**rocp                  !
!     xlat  (IMAX)     - latitude in radians, default to pi/2 -> -pi/2  !
!                        range, otherwise see in-line comment           !
!     IMAX, LM         - horizontal and vertical dimensions             !
!                                                                       !
!  outputs:                                                             !
!     o3mmr (IMAX,LM)  - output ozone profile in mass mixing ratio (g/g)!
!                                                                       !
!  module variables:                                                    !
!     k1oz, k2oz       - ozone data interpolation indices               !
!     facoz            - ozone data interpolation factor                !
!     ivflip           - control flag for direction of vertical index   !
!                                                                       !
!  usage:    call getozn                                                !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer,  intent(in) :: IMAX, LM

      real (kind=kind_phys), intent(in) :: prslk(:,:), xlat(:)

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: o3mmr(:,:)

!  ---  locals:
      real (kind=kind_phys) :: o3i(IMAX,LOZ), wk1(IMAX), deglat, elte,  &
     &                         tem, tem1, tem2, tem3, tem4, temp
      integer :: i, j, k, l, j1, j2, ll
!
!===> ...  begin here
!
      elte = blte + (JMR-1)*dlte

      do i = 1, IMAX
        deglat = xlat(i) * raddeg        ! if xlat in pi/2 -> -pi/2 range
!       deglat = 90.0 - xlat(i)*raddeg   ! if xlat in 0 -> pi range

        if (deglat > blte .and. deglat < elte) then
          tem1 = (deglat - blte) / dlte + 1
          j1   = tem1
          j2   = j1 + 1
          tem1 = tem1 - j1
        elseif (deglat <= blte) then
          j1   = 1
          j2   = 1
          tem1 = 1.0
        elseif (deglat >= elte) then
          j1   = JMR
          j2   = JMR
          tem1 = 1.0
        endif

        tem2 = 1.0 - tem1
        do j = 1, LOZ
          tem3     = tem2*o3r(j1,j,k1oz) + tem1*o3r(j2,j,k1oz)
          tem4     = tem2*o3r(j1,j,k2oz) + tem1*o3r(j2,j,k2oz)
          o3i(i,j) = tem4*facoz          + tem3*(1.0 - facoz)
        enddo
      enddo

      do l = 1, LM
        ll = l
        if (ivflip == 1) ll = LM -l + 1

        do i = 1, IMAX
          wk1(i) = prslk(i,ll)
        enddo

        do k = 1, LOZ-1
          temp = 1.0 / (pkstr(k+1) - pkstr(k))

          do i = 1, IMAX
            if (wk1(i) > pkstr(k) .and. wk1(i) <= pkstr(k+1)) then
              tem       = (pkstr(k+1) - wk1(i)) * temp
              o3mmr(I,ll) = tem * o3i(i,k) + (1.0 - tem) * o3i(i,k+1)
            endif
          enddo
        enddo

        do i = 1, IMAX
          if (wk1(i) > pkstr(LOZ)) o3mmr(i,ll) = o3i(i,LOZ)
          if (wk1(i) < pkstr(1))   o3mmr(i,ll) = o3i(i,1)
        enddo
      enddo
!
      return
!...................................
      end subroutine getozn
!! @}
!-----------------------------------

!
!........................................!
      end module module_radiation_gases  !
!========================================!
