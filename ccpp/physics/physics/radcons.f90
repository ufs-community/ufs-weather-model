!> \file radcons.f90
!! This file contains module radcons


!> \defgroup radcons GFS RRTMG Contants
!! \ingroup RRTMG
!! This module contains some of the most frequently used math and physics 
!! constants for GCM models.
!! @{
!========================================!
          module radcons                !
!........................................!
!
  use machine,      only : kind_phys
!
  implicit none
!
  public

!  ---  version tag and last revision date
      character(40), parameter ::                                       &
     &   VTAGRAD='NCEP-Radiation_driver    v5.2  Jan 2013 '
!    &   VTAGRAD='NCEP-Radiation_driver    v5.1  Nov 2012 '
!    &   VTAGRAD='NCEP-Radiation_driver    v5.0  Aug 2012 '

!>\name Constant values

!> lower limit of saturation vapor pressure (=1.0e-10)
      real (kind=kind_phys) :: QMIN
!> lower limit of specific humidity (=1.0e-7)
      real (kind=kind_phys) :: QME5
!> lower limit of specific humidity (=1.0e-7)
      real (kind=kind_phys) :: QME6
!> EPSQ=1.0e-12
      real (kind=kind_phys) :: EPSQ
!     parameter (QMIN=1.0e-10, QME5=1.0e-5,  QME6=1.0e-6,  EPSQ=1.0e-12)
      parameter (QMIN=1.0e-10, QME5=1.0e-7,  QME6=1.0e-7,  EPSQ=1.0e-12)
!     parameter (QMIN=1.0e-10, QME5=1.0e-20, QME6=1.0e-20, EPSQ=1.0e-12)

!> lower limit of toa pressure value in mb
      real, parameter :: prsmin = 1.0e-6

!> control flag for LW surface temperature at air/ground interface
!! (default=0, the value will be set in subroutine radinit)
      integer :: itsfc  =0

!> new data input control variables (set/reset in subroutines radinit/radupdate):
      integer :: month0=0,   iyear0=0,   monthd=0

!> control flag for the first time of reading climatological ozone data
!! (set/reset in subroutines radinit/radupdate, it is used only if the
!! control parameter ioznflg=0)
      logical :: loz1st =.true.

! DH* THIS MUST GO BUT NEED IT RIGHT NOW TO DEFINE LEXTOP
!> optional extra top layer on top of low ceiling models
!!\n LTP=0: no extra top layer
      integer, parameter :: LTP = 0   ! no extra top layer
!     integer, parameter :: LTP = 1   ! add an extra top layer
! *DH

!> control flag for extra top layer
      logical, parameter :: lextop = (LTP > 0)

!----------------------------
! Module variable definitions
!----------------------------
! DH* CHECK IF THIS IS NEEDED/TRUE?
!CCPP: copy from GFS_driver.F90
  real(kind=kind_phys), parameter :: con_24  =   24.0_kind_phys
  real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
  real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
  real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys
 ! real(kind=kind_phys), parameter :: qmin    =    1.0e-10




!........................................!
      end module radcons                !
!========================================!
!! @}
