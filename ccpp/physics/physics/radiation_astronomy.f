!>  \file radiation_astronomy.f
!!  This file sets up astronomical quantities for solar radiation
!!  calculations.

!  ==========================================================  !!!!!
!          'module_radiation_astronomy'  description           !!!!!
!  ==========================================================  !!!!!
!                                                                      !
!   set up astronomy quantities for solar radiation calculations.      !
!                                                                      !
!   in module 'module_radiation_astronomy', externally accessable      !
!   subroutines are listed below:                                      !
!                                                                      !
!      'sol_init'   -- initialization                                  !
!         input:                                                       !
!           ( me )                                                     !
!         output:                                                      !
!           ( none )                                                   !
!                                                                      !
!      'sol_update' -- update astronomy related quantities             !
!         input:                                                       !
!           ( jdate,kyear,deltsw,deltim,lsol_chg, me )                 !
!         output:                                                      !
!           ( slag,sdec,cdec,solcon )                                  !
!                                                                      !
!      'coszmn'     -- compute cosin of zenith angles                  !
!         input:                                                       !
!           ( xlon,sinlat,coslat,solhr,IM, me )                        !
!         output:                                                      !
!           ( coszen,coszdg )                                          !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!       'module physparam'                  in 'physparam.f'           !
!       'module physcons'                   in 'physcons.f'            !
!                                                                      !
!   program history log:                                               !
!   - a collection of programs to track solar-earth position           !
!     may  1977  ---  ray orzol    (gfdl) created program compjd to    !
!          computes julian day and fraction from year,month,dayand,time!
!     jun  1977  ---  robert white (gfdl) created program cdate to     !
!          computes calendar month, day, year from julian day value.   !
!     jul  1977  ---  robert white (gfdl) created program solar to     !
!          computes radius vector, declination and right ascension of  !
!          sun, equation of time, hour angle, fractional daylight, and !
!          latitudinal mean zenith angle.                              !
!     fall 1988  ---  hualu pan,  updated to limit the iterations in   !
!          newton method and also ccr reduced to avoid non-convergence.!
!     jul  1989  ---  kenneth campana  modified subr solar and created !
!          subr zenith for computations of effective mean cosz and     !
!          daylight fraction.                                          !
!     oct  1990  ---  yu-tai hou      created subr coszmn to replace   !
!          the latitudinal mean cosz by time mean cosz at grid points. !
!     may  1998  ---  mark iredell    y2k compliance                   !
!     dec  2003  ---  yu-tai hou      combined compjd and fcstim and   !
!          rewritten programs in fortran 90 compatable modular form.   !
!     feb  2006  ---  yu-tai hou      add 11-yr solar constant cycle   !
!     mar  2009  ---  yu-tai hou      modified solinit for climate     !
!          hindcast situation responding to ic time.                   !
!     aug  2012  ---  yu-tai hou      modified coszmn to allows sw     !
!          radiation calling interval less than 1 hr limit and linked  !
!          model time step with numb of cosz evaluations. also changed !
!          the initialization subroutine 'solinit' into two parts:     !
!          'sol_init' is called at the start of run to set up module   !
!          parameters; and 'sol_update' is called within the time      !
!          loop to check and update data sets.                         !
!     nov  2012  ---  yu-tai hou      modified control parameters thru !
!          model 'physparam'.                                          !
!     jan  2013  ---  yu-tai hou      modified to include new solar    !
!          constant tables (noaa_a0, noaa_an, cmip_an, cmip_mn)        !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



! \ingroup RRTMG
!> \defgroup module_radiation_astronomy RRTMG Astronomy Module
!> \brief This module sets up astronomical quantities for solar radiation
!!  calculations.
!! \version NCEP-Radiation_astronomy v5.2  Jan 2013
!========================================!
      module module_radiation_astronomy  !
!........................................!
!
      use physparam,         only : isolar, solar_file, kind_phys
      use physcons,          only : con_solr, con_solr_old, con_pi
      use module_iounitdef,  only : NIRADSF
!
      implicit   none
!
      private

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGAST='NCEP-Radiation_astronomy v5.2  Jan 2013 '
!    &   VTAGAST='NCEP-Radiation_astronomy v5.1  Nov 2012 '

!>\name Parameter constants
      real (kind=kind_phys), parameter :: degrad = 180.0/con_pi
      real (kind=kind_phys), parameter :: tpi    = 2.0 * con_pi
      real (kind=kind_phys), parameter :: hpi    = 0.5 * con_pi
      real (kind=kind_phys), parameter :: f12    = 12.0
      real (kind=kind_phys), parameter :: f3600  = 3600.0
      real (kind=kind_phys), parameter :: czlimt = 0.0001      ! ~ cos(89.99427)
!     real (kind=kind_phys), parameter :: pid12  = con_pi/f12  ! angle per hour
      real (kind=kind_phys), parameter :: pid12  = (2.0*asin(1.0))/f12

! \name Module variable (to be set in module_radiation_astronomy::sol_init):
      real (kind=kind_phys), public    :: solc0 = con_solr
      integer   :: isolflg = 10
      character(26) :: solar_fname = ' '

! \name Module variables (to be set in module_radiation_astronomy::sol_update)

! equation of time
      real (kind=kind_phys) :: sollag=0.0
! sine of the solar declination angle
      real (kind=kind_phys) :: sindec=0.0
! cosine of the solar declination angle
      real (kind=kind_phys) :: cosdec=0.0
! solar angle increment per interation of cosz calc
      real (kind=kind_phys) :: anginc=0.0
! saved monthly solar constants (isolflg=4 only)
      real (kind=kind_phys) :: smon_sav(12)
      data smon_sav(1:12) / 12*con_solr /

! saved year  of data used
      integer               :: iyr_sav =0
! total number of zenith angle iterations
      integer               :: nstp    =6

      public  sol_init, sol_update, coszmn


! =================
      contains
! =================

!>\ingroup module_radiation_astronomy 
!> This subroutine initializes astronomy process, and set up module
!! constants.
!!\param me         print message control flag
!>\section sol_init_gen sol_init General Algorithm
!! @{
      subroutine sol_init                                               &
     &     ( me ) !  ---  inputs
!  ---  outputs: ( none )

!  ===================================================================  !
!                                                                       !
!  initialize astronomy process, set up module constants.               !
!                                                                       !
!  inputs:                                                              !
!     me      - print message control flag                              !
!                                                                       !
!  outputs:  (to module variable)                                       !
!     ( none )                                                          !
!                                                                       !
!  external module variable: (in physparam)                             !
!   isolar    - = 0: use the old fixed solar constant in "physcon"      !
!               =10: use the new fixed solar constant in "physcon"      !
!               = 1: use noaa ann-mean tsi tbl abs-scale with cyc apprx !
!               = 2: use noaa ann-mean tsi tbl tim-scale with cyc apprx !
!               = 3: use cmip5 ann-mean tsi tbl tim-scale with cyc apprx!
!               = 4: use cmip5 mon-mean tsi tbl tim-scale with cyc apprx!
!   solar_file- external solar constant data table                      !
!                                                                       !
!  internal module variable:                                            !
!   isolflg   - internal solar constant scheme control flag             !
!   solc0     - solar constant  (w/m**2)                                !
!   solar_fname-file name for solar constant table assigned based on    !
!               the scheme control flag, isolflg.                       !
!                                                                       !
!  usage:    call sol_init                                              !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  input:
      integer,  intent(in) :: me

!  ---  output: ( none )

!  ---  local:
      logical :: file_exist
!
!===>  ...  begin here
!
      if ( me == 0 ) print *, VTAGAST    !print out version tag

!  ---  initialization
      isolflg = isolar
      solc0   = con_solr
      solar_fname = solar_file
      iyr_sav = 0
      nstp    = 6

      if ( isolar == 0 ) then
        solc0   = con_solr_old
        if ( me == 0 ) then
          print *,' - Using old fixed solar constant =', solc0
        endif
      elseif ( isolar == 10 ) then
        if ( me == 0 ) then
          print *,' - Using new fixed solar constant =', solc0
        endif
      elseif ( isolar == 1 ) then        ! noaa ann-mean tsi in absolute scale
        solar_fname(15:26) = 'noaa_a0.txt'

        if ( me == 0 ) then
          print *,' - Using NOAA annual mean TSI table in ABS scale',   &
     &            ' with cycle approximation (old values)!'
        endif

        inquire (file=solar_fname, exist=file_exist)
        if ( .not. file_exist ) then
          isolflg = 10

          if ( me == 0 ) then
            print *,'   Requested solar data file "',solar_fname,       &
     &              '" not found!'
            print *,'   Using the default solar constant value =',solc0,&
     &              ' reset control flag isolflg=',isolflg
          endif
        endif
      elseif ( isolar == 2 ) then        ! noaa ann-mean tsi in tim scale
        solar_fname(15:26) = 'noaa_an.txt'

        if ( me == 0 ) then
          print *,' - Using NOAA annual mean TSI table in TIM scale',   &
     &            ' with cycle approximation (new values)!'
        endif

        inquire (file=solar_fname, exist=file_exist)
        if ( .not. file_exist ) then
          isolflg = 10

          if ( me == 0 ) then
            print *,'   Requested solar data file "',solar_fname,       &
     &              '" not found!'
            print *,'   Using the default solar constant value =',solc0,&
     &              ' reset control flag isolflg=',isolflg
          endif
        endif
      elseif ( isolar == 3 ) then        ! cmip5 ann-mean tsi in tim scale
        solar_fname(15:26) = 'cmip_an.txt'

        if ( me == 0 ) then
          print *,' - Using CMIP5 annual mean TSI table in TIM scale',  &
     &            ' with cycle approximation'
        endif

        inquire (file=solar_fname, exist=file_exist)
        if ( .not. file_exist ) then
          isolflg = 10

          if ( me == 0 ) then
            print *,'   Requested solar data file "',solar_fname,       &
     &              '" not found!'
            print *,'   Using the default solar constant value =',solc0,&
     &              ' reset control flag isolflg=',isolflg
          endif
        endif
      elseif ( isolar == 4 ) then        ! cmip5 mon-mean tsi in tim scale
        solar_fname(15:26) = 'cmip_mn.txt'

        if ( me == 0 ) then
          print *,' - Using CMIP5 monthly mean TSI table in TIM scale', &
     &            ' with cycle approximation'
        endif

        inquire (file=solar_fname, exist=file_exist)
        if ( .not. file_exist ) then
          isolflg = 10

          if ( me == 0 ) then
            print *,'   Requested solar data file "',solar_fname,       &
     &              '" not found!'
            print *,'   Using the default solar constant value =',solc0,&
     &              ' reset control flag isolflg=',isolflg
          endif
        endif
      else                               ! selection error
        isolflg = 10

        if ( me == 0 ) then
          print *,' - !!! ERROR in selection of solar constant data',   &
     &            ' source, ISOL =',isolar
          print *,'   Using the default solar constant value =',solc0,  &
     &              ' reset control flag isolflg=',isolflg
        endif
      endif       ! end if_isolar_block
!
      return
!...................................
      end subroutine sol_init
!! @}
!-----------------------------------

!>\ingroup module_radiation_astronomy 
!> This subroutine computes solar parameters at forecast time.
!!\param jdate     ncep absolute date and time at fcst time
!!                 (yr, mon, day, t-zone, hr, min, sec, mil-sec)
!!\param kyear     usually kyear=jdate(1). if not, it is for hindcast
!!                 mode, and it is usually the init cond time and
!!                 serves as the upper limit of data can be used.
!!\param deltsw    time duration in seconds per sw calculation
!!\param deltim    timestep in seconds
!!\param lsol_chg  logical flags for change solar constant
!!\param me        print message control flag
!!\param slag               equation of time in radians
!!\param sdec, cdec         sin and cos of the solar declination angle
!!\param solcon             sun-earth distance adjusted solar constant
!!                           \f$(w/m^2)\f$
!>\section gen_sol_update sol_update General Algorithm
!! @{
!-----------------------------------
      subroutine sol_update                                             &
     &     ( jdate,kyear,deltsw,deltim,lsol_chg, me,                    &     !  ---  inputs
     &       slag, sdec, cdec, solcon                                   &     !  ---  outputs
     &      )

!  ===================================================================  !
!                                                                       !
!  sol_update computes solar parameters at forecast time                !
!                                                                       !
!  inputs:                                                              !
!     jdate(8)- ncep absolute date and time at fcst time                !
!                (yr, mon, day, t-zone, hr, min, sec, mil-sec)          !
!     kyear   - usually kyear=jdate(1). if not, it is for hindcast mode,!
!               and it is usually the init cond time and serves as the  !
!               upper limit of data can be used.                        !
!     deltsw  - time duration in seconds per sw calculation             !
!     deltim  - timestep in seconds                                     !
!     lsol_chg- logical flags for change solar constant                 !
!     me      - print message control flag                              !
!                                                                       !
!  outputs:                                                             !
!    slag          - equation of time in radians                        !
!    sdec, cdec    - sin and cos of the solar declination angle         !
!    solcon        - sun-earth distance adjusted solar constant (w/m2)  !
!                                                                       !
!                                                                       !
!  module variable:                                                     !
!   solc0   - solar constant  (w/m**2) not adjusted by earth-sun dist   !
!   isolflg - solar constant control flag                               !
!             = 0: use the old fixed solar constant                     !
!             =10: use the new fixed solar constant                     !
!             = 1: use noaa ann-mean tsi tbl abs-scale with cycle apprx !
!             = 2: use noaa ann-mean tsi tbl tim-scale with cycle apprx !
!             = 3: use cmip5 ann-mean tsi tbl tim-scale with cycle apprx!
!             = 4: use cmip5 mon-mean tsi tbl tim-scale with cycle apprx!
!   solar_fname-external solar constant data table                      !
!   sindec  - sine of the solar declination angle                       !
!   cosdec  - cosine of the solar declination angle                     !
!   anginc  - solar angle increment per iteration for cosz calc         !
!   nstp    - total number of zenith angle iterations                   !
!   smon_sav- saved monthly solar constants (isolflg=4 only)            !
!   iyr_sav - saved year  of data previously used                       !
!                                                                       !
!  usage:    call sol_update                                            !
!                                                                       !
!  subprograms called:  solar, prtime                                   !
!                                                                       !
!  external functions called: iw3jdn                                    !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  input:
      integer, intent(in) :: jdate(:), kyear, me
      logical, intent(in) :: lsol_chg

      real (kind=kind_phys), intent(in) :: deltsw, deltim

!  ---  output:
      real (kind=kind_phys), intent(out) :: slag, sdec, cdec, solcon

!  ---  locals:
      real (kind=kind_phys), parameter :: hrday = 1.0/24.0    ! frc day/hour
      real (kind=kind_phys), parameter :: minday= 1.0/1440.0  ! frc day/minute
      real (kind=kind_phys), parameter :: secday= 1.0/86400.0 ! frc day/second

      real (kind=kind_phys) :: smean, solc1, dtswh, smon(12)
      real (kind=kind_phys) :: fjd, fjd1, dlt, r1, alp

      integer :: jd, jd1, iyear, imon, iday, ihr, imin, isec
      integer :: iw3jdn
      integer :: i, iyr, iyr1, iyr2, jyr, nn, nswr, icy1, icy2, icy

      logical :: file_exist
      character :: cline*60
!
!===>  ...  begin here
!
!  --- ...  forecast time
      iyear = jdate(1)
      imon  = jdate(2)
      iday  = jdate(3)
      ihr   = jdate(5)
      imin  = jdate(6)
      isec  = jdate(7)

      if ( lsol_chg ) then   ! get solar constant from data table

        if ( iyr_sav == iyear ) then   ! same year, no new reading necessary
          if ( isolflg==4 ) then
            solc0 = smon_sav(imon)
          endif
        else                           ! need to read in new data
          iyr_sav = iyear

!  --- ...  check to see if the solar constant data file existed

          inquire (file=solar_fname, exist=file_exist)
          if ( .not. file_exist ) then
            print *,' !!! ERROR! Can not find solar constant file!!!'
            stop
          else
            iyr = iyear

            close(NIRADSF)
            open (NIRADSF,file=solar_fname,form='formatted',            &
     &                    status='old')
            rewind NIRADSF

            read (NIRADSF, * ) iyr1,iyr2,icy1,icy2,smean,cline(1:60)
!           read (NIRADSF, 24) iyr1,iyr2,icy1,icy2,smean,cline
! 24        format(4i5,f8.2,a60)

            if ( me == 0 ) then
              print *,'  Updating solar constant with cycle approx'
              print *,'   Opened solar constant data file: ',solar_fname
!check        print *, iyr1, iyr2, icy1, icy2, smean, cline
            endif

!  --- ...  check if there is a upper year limit put on the data table

!           if ( iyear /= kyear ) then
!             icy = icy1 - iyr1 + 1    ! range of the earlest cycle in data table
!             if ( kyear-iyr1 < icy ) then ! need data range at least icy years
                                           ! to perform cycle approximation
!               if ( me == 0 ) then
!                 print *,'  *** the requested year',iyear,' and upper',&
!    &                  'limit',kyear,' do not fit the range of data ', &
!    &                  'table of iyr1, iyr2 =',iyr1,iyr2
!                 print *,'      USE FIXED SOLAR CONSTANT=',con_solr
!               endif
!               solc0 = con_solr
!               isolflg = 10

!             elseif ( kyear < iyr2 ) then

!  --- ...  because the usage limit put on the historical data table,
!           skip those unused data records at first

!               i = iyr2
!               Lab_dowhile0 : do while ( i > kyear )
!                 read (NIRADSF,26) jyr, solc1
! 26              format(i4,f10.4)
!                 read (NIRADSF,*) jyr, solc1
!                 i = i - 1
!               enddo Lab_dowhile0

!               iyr2 = kyear   ! next record will serve the upper limit

!             endif   ! end if_kyear_block
!           endif   ! end if_iyear_block

!  --- ...  checking the cycle range

            if ( iyr < iyr1 ) then
              icy = icy1 - iyr1 + 1    ! range of the earlest cycle in data table
              Lab_dowhile1 : do while ( iyr < iyr1 )
                iyr = iyr + icy
              enddo Lab_dowhile1

              if ( me == 0 ) then
                print *,'   *** Year',iyear,' out of table range!',     &
     &                  iyr1, iyr2
                print *,'       Using the closest-cycle year (',iyr,')'
              endif
            elseif ( iyr > iyr2 ) then
              icy = iyr2 - icy2 + 1    ! range of the latest cycle in data table
              Lab_dowhile2 : do while ( iyr > iyr2 )
                iyr = iyr - icy
              enddo Lab_dowhile2

              if ( me == 0 ) then
                print *,'   *** Year',iyear,' out of table range!',     &
     &                  iyr1, iyr2
                print *,'       Using the closest-cycle year (',iyr,')'
              endif
            endif

!  --- ...  locate the right record for the year of data

            if ( isolflg < 4 ) then        ! use annual mean data tables
              i = iyr2
              Lab_dowhile3 : do while ( i >= iyr1 )
!               read (NIRADSF,26) jyr, solc1
! 26            format(i4,f10.4)
                read (NIRADSF,*) jyr, solc1

                if ( i == iyr .and. iyr == jyr ) then
                  solc0  = smean + solc1

                  if (me == 0) then
                    print *,' CHECK: Solar constant data used for year',&
     &                       iyr, solc1, solc0
                  endif
                  exit Lab_dowhile3
                else
!check            if(me == 0) print *,'  Skip solar const data for yr',i
                  i = i - 1
                endif
              enddo   Lab_dowhile3
            elseif ( isolflg == 4 ) then   ! use monthly mean data tables
              i = iyr2
              Lab_dowhile4 : do while ( i >= iyr1 )
!               read (NIRADSF,26) jyr, smon(:)
! 26            format(i4,12f10.4)
                read (NIRADSF,*) jyr, smon(1:12)

                if ( i == iyr .and. iyr == jyr ) then
                  do nn = 1, 12
                    smon_sav(nn) = smean + smon(nn)
                  enddo
                  solc0  = smean + smon(imon)

                  if (me == 0) then
                    print *,' CHECK: Solar constant data used for year',&
     &                      iyr,' and month',imon
                  endif
                  exit Lab_dowhile4
                else
!check            if(me == 0) print *,'  Skip solar const data for yr',i
                  i = i - 1
                endif
              enddo   Lab_dowhile4
            endif    ! end if_isolflg_block

            close ( NIRADSF )
          endif      ! end if_file_exist_block

        endif    ! end if_iyr_sav_block
      endif   ! end if_lsol_chg_block

!  --- ...  calculate forecast julian day and fraction of julian day

      jd1 = iw3jdn(iyear,imon,iday)

!  --- ...  unlike in normal applications, where day starts from 0 hr,
!           in astronomy applications, day stats from noon.

      if (ihr < 12) then
        jd1 = jd1 - 1
        fjd1= 0.5 + float(ihr)*hrday + float(imin)*minday               &
     &            + float(isec)*secday
      else
        fjd1= float(ihr - 12)*hrday + float(imin)*minday                &
     &            + float(isec)*secday
      endif

      fjd1  = fjd1 + jd1

      jd  = int(fjd1)
      fjd = fjd1 - jd

!> -# Call solar()
      call solar                                                        &
!  ---  inputs:
     &     ( jd, fjd,                                                   &
!  ---  outputs:
     &       r1, dlt, alp                                               &
     &     )

!  --- ...  calculate sun-earth distance adjustment factor appropriate to date
      solcon = solc0 / (r1*r1)

      slag   = sollag
      sdec   = sindec
      cdec   = cosdec

!  --- ...  diagnostic print out

      if (me == 0) then

!> -# Call prtime()
        call prtime                                                     &
!  ---  inputs:
     &     ( jd, fjd, dlt, alp, r1, solcon )
!  ---  outputs: ( none )

      endif

!  --- ...  setting up calculation parameters used by subr coszmn

      nswr  = nint(deltsw / deltim)         ! number of mdl t-step per sw call
      dtswh = deltsw / f3600                ! time length in hours

      if ( deltsw >= f3600 ) then           ! for longer sw call interval
        nn   = max(6, min(12, nint(f3600/deltim) ))   ! num of calc per hour
        nstp = nint(dtswh) * nn + 1                   ! num of calc per sw call
      else                                  ! for shorter sw sw call interval
        nstp = max(2, min(20, nswr)) + 1
!       nn   = nint( float(nstp-1)/dtswh )
      endif

      anginc = pid12 * dtswh / float(nstp-1)          ! solar angle inc during each calc step

      if ( me == 0 ) then
        print *,'   for cosz calculations: nswr,deltim,deltsw,dtswh =', &
     &          nswr,deltim,deltsw,dtswh,'  anginc,nstp =',anginc,nstp
      endif

!     if (me == 0) print*,'in sol_update completed sr solar'
!
      return
!...................................
      end subroutine sol_update
!-----------------------------------
!! @}

!>\ingroup module_radiation_astronomy
!> This subroutine computes radius vector, declination and right
!! ascension of sun, and equation of time.
!!\param[in] jd   julian day
!!\param[in] fjd  fraction of the julian day
!!\param[out] r1  earth-sun radius vector
!!\param[out] dlt declination of sun in radians
!!\param[out] alp right ascension of sun in radians
!>\section solar_gen solar General Algorithm
!! @{
!-----------------------------------
      subroutine solar                                                  &
     &     ( jd, fjd,                                                   &       !  ---  inputs
     &       r1, dlt, alp                                               &       !  ---  outputs
     &     )

!  ===================================================================  !
!                                                                       !
!  solar computes radius vector, declination and right ascension of     !
!  sun, and equation of time.                                           !
!                                                                       !
!  inputs:                                                              !
!    jd       - julian day                                              !
!    fjd      - fraction of the julian day                              !
!                                                                       !
!  outputs:                                                             !
!    r1       - earth-sun radius vector                                 !
!    dlt      - declination of sun in radians                           !
!    alp      - right ascension of sun in radians                       !
!                                                                       !
!  module variables:                                                    !
!    sollag   - equation of time in radians                             !
!    sindec   - sine of declination angle                               !
!    cosdec   - cosine of declination angle                             !
!                                                                       !
!  usage:    call solar                                                 !
!                                                                       !
!  external subroutines called: none                                    !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      real (kind=kind_phys), intent(in) :: fjd
      integer,               intent(in) :: jd

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: r1, dlt, alp

!  ---  locals:
      real (kind=kind_phys), parameter :: cyear = 365.25   ! days of year
      real (kind=kind_phys), parameter :: ccr   = 1.3e-6   ! iteration limit
      real (kind=kind_phys), parameter :: tpp   = 1.55     ! days between epoch and
                                                           ! perihelion passage of 1900
      real (kind=kind_phys), parameter :: svt6  = 78.035   ! days between perihelion passage
                                                           ! and march equinox of 1900
      integer,               parameter :: jdor  = 2415020  ! jd of epoch which is january
                                                           ! 0, 1900 at 12 hours ut

      real (kind=kind_phys) :: dat, t1, year, tyear, ec, angin, ador,   &
     &       deleqn, sni, tini, er, qq, e1, ep, cd, eq, date, em,       &
     &       cr, w1, tst, sun

      integer               :: jdoe, iter

!===>  ...  begin here

! --- ...  computes time in julian centuries after epoch

      t1 = float(jd - jdor) / 36525.0

! --- ...  computes length of anomalistic and tropical years (minus 365 days)

      year = 0.25964134e0 + 0.304e-5 * t1
      tyear= 0.24219879E0 - 0.614e-5 * t1

! --- ...  computes orbit eccentricity and angle of earth's inclination from t

      ec   = 0.01675104e0 - (0.418e-4 + 0.126e-6 * t1) * t1
      angin= 23.452294e0 - (0.0130125e0 + 0.164e-5 * t1) * t1

      ador = jdor
      jdoe = ador + (svt6 * cyear) / (year - tyear)

! --- ...  deleqn is updated svt6 for current date

      deleqn= float(jdoe - jd) * (year - tyear) / cyear
      year  = year + 365.0
      sni   = sin( angin / degrad )
      tini  = 1.0 / tan( angin / degrad )
      er    = sqrt( (1.0 + ec) / (1.0 - ec) )
      qq    = deleqn * tpi / year

! --- ...  determine true anomaly at equinox

      e1    = 1.0
      cd    = 1.0
      iter  = 0

      lab_do_1 : do while ( cd > ccr )

        ep    = e1 - (e1 - ec*sin(e1) - qq) / (1.0 - ec*cos(e1))
        cd    = abs(e1 - ep)
        e1    = ep
        iter  = iter + 1

        if (iter > 10) then
          write(6,*) ' ITERATION COUNT FOR LOOP 32 =', iter
          write(6,*) ' E, EP, CD =', e1, ep, cd
          exit lab_do_1
        endif

      enddo  lab_do_1

      eq   = 2.0 * atan( er * tan( 0.5*e1 ) )

! --- ...  date is days since last perihelion passage

      dat  = float(jd - jdor) - tpp + fjd
      date = mod(dat, year)

! --- ...  solve orbit equations by newton's method

      em   = tpi * date / year
      e1   = 1.0
      cr   = 1.0
      iter = 0

      lab_do_2 : do while ( cr > ccr )

        ep   = e1 - (e1 - ec*sin(e1) - em) / (1.0 - ec*cos(e1))
        cr   = abs(e1 - ep)
        e1   = ep
        iter = iter + 1

        if (iter > 10) then
          write(6,*) ' ITERATION COUNT FOR LOOP 31 =', iter
          exit lab_do_2
        endif

      enddo  lab_do_2

      w1   = 2.0 * atan( er * tan( 0.5*e1 ) )

      r1   = 1.0 - ec*cos(e1)

      sindec = sni * sin(w1 - eq)
      cosdec = sqrt( 1.0 - sindec*sindec )

      dlt  = asin( sindec )
      alp  = asin( tan(dlt)*tini )

      tst  = cos( w1 - eq )
      if (tst < 0.0) alp = con_pi - alp
      if (alp < 0.0) alp = alp + tpi

      sun  = tpi * (date - deleqn) / year
      if (sun < 0.0) sun = sun + tpi
      sollag = sun - alp - 0.03255e0
!
      return
!...................................
      end subroutine solar
!! @}
!-----------------------------------

!>\ingroup module_radiation_astronomy
!> This subroutine computes mean cos solar zenith angle over SW calling
!! interval.
!!\param xlon       grids' longitudes in radians, work both on
!!                  zonal, 0->2pi and -pi->+pi arrangements
!!\param sinlat     sine of the corresponding latitudes
!!\param coslat     cosine of the corresponding latitudes
!!\param solhr      time after 00z in hours
!!\param IM         num of grids in horizontal dimension
!!\param me         print message control flag
!!\param coszen     average of cosz for daytime only in sw call
!!                  interval
!!\param coszdg     average of cosz over entire sw call interval
!>\section coszmn_gen coszmn General Algorithm
!! @{
!-----------------------------------
      subroutine coszmn                                                 &
     &     ( xlon,sinlat,coslat,solhr, IM, me,                          &     !  ---  inputs
     &       coszen, coszdg                                             &     !  ---  outputs
     &     )

!  ===================================================================  !
!                                                                       !
!  coszmn computes mean cos solar zenith angle over sw calling interval !
!                                                                       !
!  inputs:                                                              !
!    xlon  (IM)    - grids' longitudes in radians, work both on zonal   !
!                    0->2pi and -pi->+pi arrangements                   !
!    sinlat(IM)    - sine of the corresponding latitudes                !
!    coslat(IM)    - cosine of the corresponding latitudes              !
!    solhr         - time after 00z in hours                            !
!    IM            - num of grids in horizontal dimension               !
!    me            - print message control flag                         !
!                                                                       !
!  outputs:                                                             !
!    coszen(IM)    - average of cosz for daytime only in sw call interval
!    coszdg(IM)    - average of cosz over entire sw call interval       !
!                                                                       !
!  module variables:                                                    !
!    sollag        - equation of time                                   !
!    sindec        - sine of the solar declination angle                !
!    cosdec        - cosine of the solar declination angle              !
!    anginc        - solar angle increment per iteration for cosz calc  !
!    nstp          - total number of zenith angle iterations            !
!                                                                       !
!  usage:    call comzmn                                                !
!                                                                       !
!  external subroutines called: none                                    !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: IM, me

      real (kind=kind_phys), intent(in) :: sinlat(:), coslat(:),        &
     &       xlon(:), solhr

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: coszen(:), coszdg(:)

!  ---  locals:
      real (kind=kind_phys) :: coszn, cns, ss, cc, solang, rstp

      integer :: istsun(IM), i, it, j, lat

!===>  ...  begin here

      solang = pid12 * (solhr - f12)         ! solar angle at present time
      rstp = 1.0 / float(nstp)

      do i = 1, IM
        coszen(i) = 0.0
        istsun(i) = 0
      enddo

      do it = 1, nstp
        cns = solang + float(it-1)*anginc + sollag

        do i = 1, IM
          ss  = sinlat(i) * sindec
          cc  = coslat(i) * cosdec

          coszn = ss + cc * cos(cns + xlon(i))
          coszen(i) = coszen(i) + max(0.0, coszn)
          if (coszn > czlimt) istsun(i) = istsun(i) + 1
        enddo
      enddo

!  --- ...  compute time averages

      do i = 1, IM
        coszdg(i) = coszen(i) * rstp
        if (istsun(i) > 0) coszen(i) = coszen(i) / istsun(i)
      enddo
!
      return
!...................................
      end subroutine coszmn
!! @}
!-----------------------------------

!>\ingroup module_radiation_astronomy
!> This subroutine prints out forecast date, time, and astronomy
!! quantities.
!!\param[in] jd    forecast julian day
!!\param[in] fjd   forecast fraction of julian day
!!\param[in] dlt   declination angle of sun in radians
!!\param[in] alp   right ascension of sun in radians
!!\param[in] r1    earth-sun radius vector in meter
!!\param[in] solc  solar constant in \f$w/m^2\f$
!>\section prtime_gen prtime General Algorithm
!! @{
!-----------------------------------
      subroutine prtime                                                 &
     &     ( jd, fjd, dlt, alp, r1, solc                                &    !  ---  inputs
     &     )                                                            !  ---  outputs: ( none )

!  ===================================================================  !
!                                                                       !
!  prtime prints out forecast date, time, and astronomy quantities.     !
!                                                                       !
!  inputs:                                                              !
!    jd       - forecast julian day                                     !
!    fjd      - forecast fraction of julian day                         !
!    dlt      - declination angle of sun in radians                     !
!    alp      - right ascension of sun in radians                       !
!    r1       - earth-sun radius vector in meter                        !
!    solc     - solar constant in w/m^2                                 !
!                                                                       !
!  outputs:   ( none )                                                  !
!                                                                       !
!  module variables:                                                    !
!    sollag   - equation of time in radians                             !
!                                                                       !
!  usage:    call prtime                                                !
!                                                                       !
!  external subroutines called: w3fs26                                  !
!                                                                       !
!  ===================================================================  !
!
      implicit none

!  ---  inputs:
      integer, intent(in) :: jd

      real (kind=kind_phys), intent(in) :: fjd, dlt, alp, r1, solc

!  ---  outputs: ( none )

!  ---  locals:
      real (kind=kind_phys), parameter :: sixty  = 60.0

      character(LEN=1),     parameter :: sign   = '-'
      character(LEN=1),     parameter :: sigb   = ' '

      character(LEN=1)     :: dsig
      character(LEN=4)     :: month(12)

      data month / 'JAN.','FEB.','MAR.','APR.','MAY ','JUNE',           &
     &             'JULY','AUG.','SEP.','OCT.','NOV ','DEC.' /

      integer               :: iday, imon, iyear, ihr, ltd, ltm,        &
     &                         ihalp, iyy, jda, mfjd, idaywk, idayyr
      real (kind=kind_phys) :: xmin, dltd, dltm, dlts, halp, ymin,      &
     &                         asec, eqt, eqsec

!===>  ...  begin here

!  --- ...  get forecast hour and minute from fraction of julian day

      if (fjd >= 0.5) then
        jda = jd + 1
        mfjd= nint( fjd*1440.0 )
        ihr = mfjd / 60 - 12
        xmin= float(mfjd) - (ihr + 12)*sixty
      else
        jda = jd
        mfjd= nint( fjd*1440.0 )
        ihr = mfjd / 60 + 12
        xmin= float(mfjd) - (ihr - 12)*sixty
      endif

!  --- ...  get forecast year, month, and day from julian day

      call w3fs26(jda, iyear,imon,iday, idaywk,idayyr)

!  -- ...  compute solar parameters

      dltd = degrad * dlt
      ltd  = dltd
      dltm = sixty * (abs(dltd) - abs(float(ltd)))
      ltm  = dltm
      dlts = sixty * (dltm - float(ltm))

      if ((dltd < 0.0) .and. (ltd == 0.0)) then
        dsig = sign
      else
        dsig = sigb
      endif

      halp = 6.0 * alp / hpi
      ihalp= halp
      ymin = abs(halp - float(ihalp)) * sixty
      iyy  = ymin
      asec = (ymin - float(iyy)) * sixty

      eqt  = 228.55735 * sollag
      eqsec= sixty * eqt

      print 101, iday, month(imon), iyear, ihr, xmin, jd, fjd
 101  format('0 FORECAST DATE',9x,i3,a5,i6,' AT',i3,' HRS',f6.2,' MINS'/&
     &       '  JULIAN DAY',12x,i8,2x,'PLUS',f11.6)

      print 102, r1, halp, ihalp, iyy, asec
 102  format('  RADIUS VECTOR',9x,f10.7/'  RIGHT ASCENSION OF SUN',     &
     &       f12.7,' HRS, OR',i4,' HRS',i4,' MINS',f6.1,' SECS')

      print 103, dltd, dsig, ltd, ltm, dlts, eqt, eqsec, sollag, solc
 103  format('  DECLINATION OF THE SUN',f12.7,' DEGS, OR ',a1,i3,       &
     &       ' DEGS',i4,' MINS',f6.1,' SECS'/'  EQUATION OF TIME',6x,   &
     &       f12.7,' MINS, OR',f10.2,' SECS, OR',f9.6,' RADIANS'/       &
     &       '  SOLAR CONSTANT',8X,F12.7,' (DISTANCE AJUSTED)'//)

!
      return
!...................................
      end subroutine prtime
!! @}
!-----------------------------------

!
!...........................................!
      end module module_radiation_astronomy !
!===========================================!
