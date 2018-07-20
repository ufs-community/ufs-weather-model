!!!!!  ==========================================================  !!!!!
!!!!!             module "module_iounitdef description             !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!     this module defines fortran unit numbers for input/output data   !
!     files for the ncep gfs model.                                    !
!                                                                      !
!      name      type    description                         unit no.  !
!     ---------------------------------------------------------------  !
!     NISIGI  - input,  sigma file 1                            11     !
!     NISIGI2 - input,  sigma file 2                            12     !
!     NISFCI  - input,  surface initial data                    14     !
!                                                                      !
!     NIMTNVR - input,  montain variance file                   24     !
!     NIDTBTH - input,  equivalent potential temperature file   27     !
!     NICO2TR - input,  co2 transm table for gfdl-lw only       15     !
!     NICO2CN - input,  monthly/yearly 2-d co2 data   (shared)  102    !
!     NIO3PRD - input,  ozone production climatology            28     !
!     NIO3LOS - input,  ozone destruction climatology           29     !
!     NIO3CLM - input,  ozone climatology distribution          48     !
!     NINAMSF - input,  namelist for surface file               35     !
!     NISFCYC - input,  surface cycle files                     101    !
!     NIRADSF - input,  radiation surface data files  (shared)  102    !
!     NICLTUN - input,  cloud tuning table                      43     !
!     NIMICPH - input,  micro physics data file                 1      !
!     NIAERCM - input,  aerosols climatology          (shared)  102    !
!                                                                      !
!     NOSIGR1 - output, first time level sigma restart file     51     !
!     NOSIGR2 - output, second time level sigma restart file    52     !
!     NOSFCR  - output, surface restart file                    53     !
!     NOSIGF  - output, sigma file for post process             61     !
!     NOSFCF  - output, surface file for post process           62     !
!     NOFLXF  - output, flux file for post process              63     !
!     NOD3DF  - output, 3-d  file for post process              64     !
!     NOAERF  - output, 2-d  file for post process              65     !
!hchuang code change [+1L]
!     NOG3DF  - output, 3-d  file for GFS-GOCART specific       69     !
!     NOLOGF  - output, log  file                               99     !
!                                                                      !
!     NIOFRAD - in/out, temperary radiation data file (shared)  16     !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!

!========================================!
      module module_iounitdef            !
!........................................!
!
      implicit   none
!
      public

!  --- ...  input units

      integer, parameter :: NISIGI  = 11
      integer, parameter :: NISIGI2 = 12
      integer, parameter :: NISFCI  = 14
      integer, parameter :: NICO2TR = 15
      integer, parameter :: NICO2CN = 102
      integer, parameter :: NIMTNVR = 24
      integer, parameter :: NIDTBTH = 27
      integer, parameter :: NIO3PRD = 28
      integer, parameter :: NIO3LOS = 29
      integer, parameter :: NINAMSF = 35
      integer, parameter :: NICLTUN = 43
      integer, parameter :: NIO3CLM = 48
      integer, parameter :: NIMICPH = 1
      integer, parameter :: NISFCYC = 101
      integer, parameter :: NIAERCM = 102
      integer, parameter :: NIRADSF = 102

!  --- ... output units

      integer, parameter :: NOSIGR1 = 51
      integer, parameter :: NOSIGR2 = 52
      integer, parameter :: NOSFCR  = 53
      integer, parameter :: NOSIGF  = 61
      integer, parameter :: NOSFCF  = 62
      integer, parameter :: NOFLXF  = 63
      integer, parameter :: NOD3DF  = 64
      integer, parameter :: NOAERF  = 65    ! for g2d_fld
!hchuang code change [+1L]
      integer, parameter :: NOG3DF  = 69
      integer, parameter :: NOLOGF  = 99

!  --- ...  in/out units

      integer, parameter :: NIOFRAD = 16

!
!........................................!
      end module module_iounitdef        !
!========================================!
