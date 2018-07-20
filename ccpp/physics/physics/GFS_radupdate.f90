!>\file GFS_radupdate.f90
!! This file calls many update subroutines to check and update radiation required but 
!! time varying data sets and module variables.
       module GFS_radupdate
       contains

!>\defgroup GFS_radupdate_run RRTMG Update
!!\ingroup RRTMG
!>\brief This subroutine calls many update subroutines to check and update radiation
!! required but time varying data sets and module variables.
!>\section gen_radupdate General Algorithm
!> @{
       subroutine GFS_radupdate_run ( idate,jdate,deltsw,deltim,lsswr, me, &
     &       slag,sdec,cdec,solcon, ictmflg, isolar)

! =================   subprogram documentation block  ================ !
!                                                                      !
! subprogram:   radupdate   calls many update subroutines to check and !
!   update radiation required but time varying data sets and module    !
!   variables.                                                         !
!                                                                      !
! usage:        call radupdate                                         !
!                                                                      !
! attributes:                                                          !
!   language:  fortran 90                                              !
!   machine:   ibm sp                                                  !
!                                                                      !
!  ====================  definition of variables ====================  !
!                                                                      !
! input parameters:                                                    !
!   idate(8)      : ncep absolute date and time of initial condition   !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)     !
!   jdate(8)       : ncep absolute date and time at fcst time          !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)     !
!   deltsw         : sw radiation calling frequency in seconds         !
!   deltim         : model timestep in seconds                         !
!   lsswr          : logical flags for sw radiation calculations       !
!   me             : print control flag                                !
!                                                                      !
!  outputs:                                                            !
!   slag           : equation of time in radians                       !
!   sdec, cdec     : sin and cos of the solar declination angle        !
!   solcon         : sun-earth distance adjusted solar constant(w/m2)  !
!                                                                      !
!  external module variables:                                          !
!   isolar   : solar constant cntrl  (in module physparam)             !
!              = 0: use the old fixed solar constant in "physcon"      !
!              =10: use the new fixed solar constant in "physcon"      !
!              = 1: use noaa ann-mean tsi tbl abs-scale with cycleapprx!
!              = 2: use noaa ann-mean tsi tbl tim-scale with cycleapprx!
!              = 3: use cmip5 ann-mean tsi tbl tim-scale with cyclapprx!
!              = 4: use cmip5 mon-mean tsi tbl tim-scale with cyclapprx!
!   ictmflg  : =yyyy#, external data ic time/date control flag         !
!              =   -2: same as 0, but superimpose seasonal cycle       !
!                      from climatology data set.                      !
!              =   -1: use user provided external data for the         !
!                      forecast time, no extrapolation.                !
!              =    0: use data at initial cond time, if not           !
!                      available, use latest, no extrapolation.        !
!              =    1: use data at the forecast time, if not           !
!                      available, use latest and extrapolation.        !
!              =yyyy0: use yyyy data for the forecast time,            !
!                      no further data extrapolation.                  !
!              =yyyy1: use yyyy data for the fcst. if needed, do       !
!                      extrapolation to match the fcst time.           !
!                                                                      !
!  module variables:                                                   !
!   loz1st   : first-time clim ozone data read flag                    !
!                                                                      !
!  subroutines called: sol_update, aer_update, gas_update              !
!                                                                      !
! ===================================================================  !
!
      use radcons
      use module_radiation_astronomy,              only:  sol_update
      use module_radiation_gases,                  only:  gas_update
      use module_radiation_aerosols,               only:  aer_update

      implicit none

!  ---  inputs:
      integer, intent(in) :: idate(:), jdate(:), me
      logical, intent(in) :: lsswr

      real (kind=kind_phys), intent(in) :: deltsw, deltim

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: slag, sdec, cdec, solcon

!  ---  locals:
      integer :: iyear, imon, iday, ihour
      integer :: kyear, kmon, kday, khour

      logical :: lmon_chg       ! month change flag
      logical :: lco2_chg       ! cntrl flag for updating co2 data
      logical :: lsol_chg       ! cntrl flag for updating solar constant
      integer :: ictmflg, isolar
!
!===> ...  begin here
!
!> -# Set up time stamp at fcst time and that for green house gases
!! (currently co2 only)
!  --- ...  time stamp at fcst time

      iyear = jdate(1)
      imon  = jdate(2)
      iday  = jdate(3)
      ihour = jdate(5)

!  --- ...  set up time stamp used for green house gases (** currently co2 only)

      if ( ictmflg==0 .or. ictmflg==-2 ) then  ! get external data at initial condition time
        kyear = idate(1)
        kmon  = idate(2)
        kday  = idate(3)
        khour = idate(5)
      else                           ! get external data at fcst or specified time
        kyear = iyear
        kmon  = imon
        kday  = iday
        khour = ihour
      endif   ! end if_ictmflg_block

      if ( month0 /= imon ) then
        lmon_chg = .true.
        month0   = imon
      else
        lmon_chg = .false.
      endif

!> -# Call module_radiation_astronomy::sol_update(), yearly update, no
!! time interpolation.
      if (lsswr) then

        if ( isolar == 0 .or. isolar == 10 ) then
          lsol_chg = .false.
        elseif ( iyear0 /= iyear ) then
          lsol_chg = .true.
        else
          lsol_chg = ( isolar==4 .and. lmon_chg )
        endif
        iyear0 = iyear

        call sol_update    &
!  ---  inputs:
     &     ( jdate,kyear,deltsw,deltim,lsol_chg, me,   &
!  ---  outputs:
     &       slag,sdec,cdec,solcon                     &
     &     )

      endif  ! end_if_lsswr_block

!> -# Call module_radiation_aerosols::aer_update(), monthly update, no
!! time interpolation
      if ( lmon_chg ) then
        call aer_update ( iyear, imon, me )
      endif

!> -# Call co2 and other gases update routine:
!! module_radiation_gases::gas_update()
      if ( monthd /= kmon ) then
        monthd = kmon
        lco2_chg = .true.
      else
        lco2_chg = .false.
      endif

      call gas_update ( kyear,kmon,kday,khour,loz1st,lco2_chg, me )

      if ( loz1st ) loz1st = .false.

!> -# Call surface update routine (currently not needed)
!     call sfc_update ( iyear, imon, me )

!> -# Call clouds update routine (currently not needed)
!     call cld_update ( iyear, imon, me )
!
      return
!...................................
      end subroutine GFS_radupdate_run
!-----------------------------------
!> @}

      end module GFS_radupdate
