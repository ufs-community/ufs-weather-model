module physics_diag_layer

!------------------------------------------------------------------------------------------!
!                                                                                          !
! This module populates the IPD_Diag container withe the elements from GFS physics that    !
! are to be output via the write component in the NEMS system.  The IPD_Diag container     !
! contains properties from the GFS_diag_type, GFS_sfcprop_type and 
!                                                                                          !
!------------------------------------------------------------------------------------------!

  use machine,                   only: kind_phys
  use IPD_typedefs,              only: IPD_diag_type
  use physics_abstraction_layer, only: control_type,  statein_type,  &
                                       stateout_type, sfcprop_type,  &
                                       coupling_type, grid_type,     &
                                       tbd_type,      cldprop_type,  &
                                       radtend_type,  intdiag_type,  &
                                       init_type

  public diag_populate

  CONTAINS
!*******************************************************************************************

!----------------------
! GFS_populate_IPD_Diag
!----------------------
  subroutine diag_populate (IPD_Diag, Model, Statein, Stateout, Sfcprop, Coupling,     &
                            Grid, Tbd, Cldprop, Radtend, Diag, Init_parm)
!------------------------------------------------------------------------------------------!
!   IPD_METADATA                                                                           !
!     IPD_Diag%name           [char*32 ]   variable name in source  [char*32]              !
!     IPD_Diag%output_name    [char*32 ]   output name for variable [char*32]              !
!     IPD_Diag%mod_name       [char*32 ]   module name (e.g. physics, radiation, etc)      !
!     IPD_Diag%file_name      [char*32 ]   output file name for variable                   !
!     IPD_Diag%desc           [char*128]   long description of field                       !
!     IPD_Diag%unit           [char*32 ]   units associated with fields                    !
!     IPD_Diag%type_stat_proc [char*32 ]   type of statistic processing:                   !
!                                          average, accumulation, maximal, minimal, etc.   !
!     IPD_Diag%level_type     [char*32 ]   vertical level of the field                     !
!     IPD_Diag%level          [int*4   ]   vertical level(s)                               !
!     IPD_Diag%cnvfac         [real*8  ]   conversion factors to output in specified units !
!     IPD_Diag%zhour          [real*8  ]   forecast hour when bucket was last emptied      !
!     IPD_Diag%fcst_hour      [real*8  ]   current forecast hour (same as fhour)           !
!     IPD_Diag%data(nb)%var2p(:)       [real*8  ]   pointer to 2D data [=> null() for a 3D field]   !
!     IPD_Diag%data(nb)%var3p(:,:)     [real*8  ]   pointer to 3D data [=> null() for a 2D field]   !
!------------------------------------------------------------------------------------------!

      implicit none
!
!  ---  interface variables
    type(IPD_diag_type),        intent(inout) :: IPD_Diag(:)
    type(control_type),         intent(in)    :: Model
    type(statein_type),         intent(in)    :: Statein(:)
    type(stateout_type),        intent(in)    :: Stateout(:)
    type(sfcprop_type),         intent(in)    :: Sfcprop(:)
    type(coupling_type),        intent(in)    :: Coupling(:)
    type(grid_type),            intent(in)    :: Grid(:)
    type(tbd_type),             intent(in)    :: Tbd(:)
    type(cldprop_type),         intent(in)    :: Cldprop(:)
    type(radtend_type),         intent(in)    :: Radtend(:)
    type(intdiag_type),         intent(in)    :: Diag(:)
    type(init_type),            intent(in)    :: Init_parm

    !--- local variabls
    integer :: idx, nblks, nb, num
    real(kind=kind_phys), parameter :: cn_one = 1._kind_phys
    real(kind=kind_phys), parameter :: cn_100 = 100._kind_phys
    real(kind=kind_phys), parameter :: cn_th  = 1000._kind_phys
    real(kind=kind_phys), parameter :: cn_hr  = 3600._kind_phys

    real(kind=kind_phys), pointer :: var1(:) => null()
     
    nblks = size(Init_parm%blksz)

    !--- initialize GFS_diag
    IPD_Diag(:)%name           = ' '
    IPD_Diag(:)%output_name    = ' '
    IPD_Diag(:)%mod_name       = ' '
    IPD_Diag(:)%file_name      = ' '
    IPD_Diag(:)%desc           = ' '
    IPD_Diag(:)%unit           = ' '
    IPD_Diag(:)%type_stat_proc = ' '
    IPD_Diag(:)%level_type     = ' '
    IPD_Diag(:)%level          = 1
    IPD_Diag(:)%cnvfac         = cn_one
    IPD_Diag(:)%zhour          = Model%zhour
    IPD_Diag(:)%fcst_hour      = Model%fhour
    do idx = 1,size(IPD_Diag,1)
      allocate (IPD_Diag(idx)%data(nblks))
    enddo

    idx = 0

    ! IPD DIAG CONTAINER DATA
    !--- FLUXR
    !--- This data array contains 33 fields including radiation flux, 
    !--- cloud cover, pressure and many other fields, suggest to 
    !--- split into 2D fields:
    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr'
    IPD_Diag(idx)%output_name    = ' '
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'flux from radiation'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = ' '
    IPD_Diag(idx)%level          = Model%nfxr
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%fluxr
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr1'
    IPD_Diag(idx)%output_name    = 'ulwrf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Upward long wave radiation flux at top of atmosphere'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,1)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr2'
    IPD_Diag(idx)%output_name    = 'uswrf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Upward solar radiation flux at top of atmosphere'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,2)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr3'
    IPD_Diag(idx)%output_name    = 'uswrf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Upward solar radiation flux at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,3)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr4'
    IPD_Diag(idx)%output_name    = 'dswrf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Downward solar radiation flux at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,4)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr5'
    IPD_Diag(idx)%output_name    = 'tcdc'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Total cloud cover at high cloud layer'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'high_cloud_lyr'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,5)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr6'
    IPD_Diag(idx)%output_name    = 'tcdc'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Total cloud cover at middle cloud layer'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'mid_cloud_lyr'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,6)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr7'
    IPD_Diag(idx)%output_name    = 'tcdc'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Total cloud cover at low cloud layer'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'low_cloud_lyr'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,7)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr8'
    IPD_Diag(idx)%output_name    = 'pres'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Pressure at high cloud top level'
    IPD_Diag(idx)%unit           = 'pa'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtimsw'
    IPD_Diag(idx)%level_type     = 'high_cloud_top_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,8)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr9'
    IPD_Diag(idx)%output_name    = 'pres'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Pressure at middle cloud top level'
    IPD_Diag(idx)%unit           = 'pa'
    IPD_Diag(idx)%type_stat_proc = 'acc_cloud'
    IPD_Diag(idx)%level_type     = 'mid_cloud_top_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,9)
    enddo
    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr10'
    IPD_Diag(idx)%output_name    = 'pres'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Pressure at low cloud top level'
    IPD_Diag(idx)%unit           = 'pa'
    IPD_Diag(idx)%type_stat_proc = 'acc_cloud'
    IPD_Diag(idx)%level_type     = 'low_cloud_top_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,10)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr11'
    IPD_Diag(idx)%output_name    = 'pres'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Pressure at high cloud bottom level'
    IPD_Diag(idx)%unit           = 'pa'
    IPD_Diag(idx)%type_stat_proc = 'acc_cloud'
    IPD_Diag(idx)%level_type     = 'high_cloud_bot_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,11)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr12'
    IPD_Diag(idx)%output_name    = 'pres'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Pressure at middle cloud bottom level'
    IPD_Diag(idx)%unit           = 'pa'
    IPD_Diag(idx)%type_stat_proc = 'acc_cloud'
    IPD_Diag(idx)%level_type     = 'mid_cloud_bot_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,12)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr13'
    IPD_Diag(idx)%output_name    = 'pres'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Pressure at low cloud bot level'
    IPD_Diag(idx)%unit           = 'pa'
    IPD_Diag(idx)%type_stat_proc = 'acc_cloud'
    IPD_Diag(idx)%level_type     = 'low_cloud_bot_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,13)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr14'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Temperature at high cloud top level'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'acc_cloud'
    IPD_Diag(idx)%level_type     = 'high_cloud_top_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,14)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr15'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Temperature at middle cloud top level'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'acc_cloud'
    IPD_Diag(idx)%level_type     = 'mid_cloud_top_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,15)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr16'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Temperature at low cloud top level'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'acc_cloud'
    IPD_Diag(idx)%level_type     = 'low_cloud_top_lvl'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,16)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr17'
    IPD_Diag(idx)%output_name    = 'tcdc'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Total cloud cover at total atmospheric column'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'entire_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,17)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr18'
    IPD_Diag(idx)%output_name    = 'tcdc'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Total cloud cover (precent) at boundary layer cloud layer'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'bound_lyr_cloud_lyr'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,18)
    enddo

    !--- fluxr19 and fluxr20 are replaced with the surface temperature
    !                        adjusted quantities for output

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr21'
    IPD_Diag(idx)%output_name    = 'duvb'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'UV-B downward solar flux (w/m**2) at land sea surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,21)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr22'
    IPD_Diag(idx)%output_name    = 'cduvb'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Clear sky UV-B downward solar flux (w/m**2) at land sea surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,22)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr23'
    IPD_Diag(idx)%output_name    = 'dswrf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Downward solar radiation flux (w/m**2) at top of atmosphere'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,23)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr24'
    IPD_Diag(idx)%output_name    = 'vbdsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Downward sw uv+vis beam radiation flux (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'surfaces'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,24)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr25'
    IPD_Diag(idx)%output_name    = 'vddsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Downward sw uv+vis diffuse radiation flux (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'surfaces'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,25)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr26'
    IPD_Diag(idx)%output_name    = 'nbdsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Downward sw nir beam radiation flux (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'surfaces'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,26)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr27'
    IPD_Diag(idx)%output_name    = 'nddsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Downward sw nir diffuse radiation flux (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_timsw'
    IPD_Diag(idx)%level_type     = 'surfaces'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,27)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr28'
    IPD_Diag(idx)%output_name    = 'csulf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Clear sky upward long wave radiation flux (w/m**2) at top of atmosphere'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,28)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr29'
    IPD_Diag(idx)%output_name    = 'csusf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Clear sky upward solar radiation flux (w/m**2) at top of atmosphere'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtimer'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,29)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr30'
    IPD_Diag(idx)%output_name    = 'csdlf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Clear sky downward long wave radiation flux (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtimer'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,30)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr31'
    IPD_Diag(idx)%output_name    = 'csusf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Clear sky upward solar radiation flux (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtimer'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,31)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr32'
    IPD_Diag(idx)%output_name    = 'csdsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Clear sky downward solar radiation flux (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtimer'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,32)
    enddo

    idx = idx + 1
    IPD_Diag(idx)%name           = 'fluxr33'
    IPD_Diag(idx)%output_name    = 'csulf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Clear sky upward long wave radiation flux (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtimer'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%fluxr(:,33)
    enddo
    !--- done with FLUXR

    !--- dlwsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dlwsfc'
    IPD_Diag(idx)%output_name    = 'dlwrf'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'downward longwave flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dlwsfc
    enddo
!---need to convert to "ave" for output??

    !---ulwsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'ulwsfc'
    IPD_Diag(idx)%output_name    = 'ulwrf'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'upward long wave flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%ulwsfc
    enddo
!--- need to convert to "ave" for output??

    ! COUPLING FIELDS
    !---nirbmdi
    idx = idx + 1
    IPD_Diag(idx)%name           = 'nirbmdi'
    IPD_Diag(idx)%output_name    = 'nbdsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Downward short wave nir beam radiation flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'Inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Coupling(nb)%nirbmdi
    enddo

    !---nirdfdi
    idx = idx + 1
    IPD_Diag(idx)%name           = 'nirdfdi'
    IPD_Diag(idx)%output_name    = 'nddsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Downward short wave nir diffuse radiation flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'Inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Coupling(nb)%nirdfdi
    enddo
!--- need to convert to "ave" for output

    !---visbmdi
    idx = idx + 1
    IPD_Diag(idx)%name           = 'visbmdi'
    IPD_Diag(idx)%output_name    = 'vbdsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'Downward short wave uv+vis beam radiation flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Coupling(nb)%visbmdi
    enddo

    !---visdfdi
    idx = idx + 1
    IPD_Diag(idx)%name           = 'visdfdi'
    IPD_Diag(idx)%output_name    = 'vddsf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'cpl '
    IPD_Diag(idx)%desc           = 'Downward short wave uv+vis diffuse radiation flux [W/m**2] at surface '
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Coupling(nb)%visdfdi
    enddo

    !---nirbmui
    idx = idx + 1
    IPD_Diag(idx)%name           = 'nirbmui'
    IPD_Diag(idx)%output_name    = 'nbusf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'Upward short wave nir beam radiation flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Coupling(nb)%nirbmui
    enddo

    !---nirdfui
    idx = idx + 1
    IPD_Diag(idx)%name           = 'nirdfui'
    IPD_Diag(idx)%output_name    = 'ndusf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'Upward short wave nir beam radiation flux [W/m**2] at surface '
    IPD_Diag(idx)%unit           = ' '
    IPD_Diag(idx)%type_stat_proc = ' '
    IPD_Diag(idx)%level_type     = ' '
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Coupling(nb)%nirdfui
    enddo

    !---visbmui
    idx = idx + 1
    IPD_Diag(idx)%name           = 'visbmui'
    IPD_Diag(idx)%output_name    = 'vbusf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'Upward short wave uv+vis beam radiation flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = ' '
    IPD_Diag(idx)%type_stat_proc = ' '
    IPD_Diag(idx)%level_type     = ' '
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Coupling(nb)%visbmui
    enddo

    !---visdfui
    idx = idx + 1
    IPD_Diag(idx)%name           = 'visdfui'
    IPD_Diag(idx)%output_name    = 'vdusf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'Upward short wave uv+vis diffuse radiation flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2 '
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Coupling(nb)%visdfui
    enddo
    ! END COUPLING FIELDS

    !---topfsw%upfxc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'sw_upfxc'
    IPD_Diag(idx)%output_name    = 'uswrf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'upward solar radiation flux [w/m**2] at top of atmosphere'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
#ifdef __GFORTRAN__
       Diag(nb)%topfsw_upfxc_gnufix(:) = Diag(nb)%topfsw(:)%upfxc
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topfsw_upfxc_gnufix
#else
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topfsw(:)%upfxc
#endif
    enddo

    !---topfsw%dnfxc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'sw_dnfxc'
    IPD_Diag(idx)%output_name    = 'dswrf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'downward solar radiation flux [w/m**2] at top of atmosphere'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
#ifdef __GFORTRAN__
       Diag(nb)%topfsw_dnfxc_gnufix(:) = Diag(nb)%topfsw(:)%dnfxc
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topfsw_dnfxc_gnufix
#else
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topfsw(:)%dnfxc
#endif
    enddo

    !---topfsw%upfx0
    idx = idx + 1
    IPD_Diag(idx)%name           = 'sw_upfx0'
    IPD_Diag(idx)%output_name    = 'csusf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'Clear sky upward solar radiation flux [w/m**2] at top of atmosphere'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
#ifdef __GFORTRAN__
       Diag(nb)%topfsw_upfx0_gnufix(:) = Diag(nb)%topfsw(:)%upfx0
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topfsw_upfx0_gnufix
#else
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topfsw(:)%upfx0
#endif
    enddo

    !---topflw%upfxc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'lw_upfxc'
    IPD_Diag(idx)%output_name    = 'ulwrf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'upward long wave radiation flux [w/m**2] at top of atmosphere '
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
#ifdef __GFORTRAN__
       Diag(nb)%topflw_upfxc_gnufix(:) = Diag(nb)%topflw(:)%upfxc
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topflw_upfxc_gnufix
#else
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topflw(:)%upfxc
#endif
    enddo

!--- clear sky down long wave is missing?

    !---topflw%upfx0
    idx = idx + 1
    IPD_Diag(idx)%name           = 'lw_upfx0'
    IPD_Diag(idx)%output_name    = 'csulf'
    IPD_Diag(idx)%mod_name       = 'radiation'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'Clear sky upward long wave radiation flux [w/m**2] at top of atmosphere '
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'top_of_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
#ifdef __GFORTRAN__
       Diag(nb)%topflw_upfx0_gnufix(:) = Diag(nb)%topflw(:)%upfx0
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topflw_upfx0_gnufix
#else
       IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%topflw(:)%upfx0
#endif
    enddo

    !---srunoff
    idx = idx + 1
    IPD_Diag(idx)%name           = 'srunoff'
    IPD_Diag(idx)%output_name    = 'ssrun'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'storm water runoff [kg/m**2] at surface'
    IPD_Diag(idx)%unit           = 'kg/m**2 '
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = 1.e3
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%srunoff
    enddo

    !---evbsa
    idx = idx + 1
    IPD_Diag(idx)%name           = 'evbsa'
    IPD_Diag(idx)%output_name    = 'EVBSA'
    IPD_Diag(idx)%mod_name       = 'land'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Direct Evaporation [w/m**2] from Bare Soil'
    IPD_Diag(idx)%unit           = 'w/m**2 '
    IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%evbsa
    enddo

    !---evcwa
    idx = idx + 1
    IPD_Diag(idx)%name           = 'evcwa'
    IPD_Diag(idx)%output_name    = 'evcw'
    IPD_Diag(idx)%mod_name       = 'land'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Canopy water evaporation'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
    IPD_Diag(idx)%level_type     = 'surface '
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%evcwa
    enddo

    !---snohfa
    idx = idx + 1
    IPD_Diag(idx)%name           = 'snohfa'
    IPD_Diag(idx)%output_name    = 'snohf'
    IPD_Diag(idx)%mod_name       = 'land'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Snow Phase Change Heat Flux [w/m**2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%snohfa
    enddo

    !---snowca
    idx = idx + 1
    IPD_Diag(idx)%name           = 'snowca'
    IPD_Diag(idx)%output_name    = 'snowc'
    IPD_Diag(idx)%mod_name       = 'land'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Snow cover (fraction) at surface'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%snowca
    enddo

    !---sbsnoa
    idx = idx + 1
    IPD_Diag(idx)%name           = 'sbsnoa'
    IPD_Diag(idx)%output_name    = 'sbsno'
    IPD_Diag(idx)%mod_name       = 'land'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Snow sublimation [w/m^2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
    IPD_Diag(idx)%level_type     = 'surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%sbsnoa
    enddo

    !---transa
    idx = idx + 1
    IPD_Diag(idx)%name           = 'transa'
    IPD_Diag(idx)%output_name    = 'trans'
    IPD_Diag(idx)%mod_name       = 'land'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Transpiration [w/m**2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%transa
    enddo

    !---soilm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'soilm'
    IPD_Diag(idx)%output_name    = 'soilm'
    IPD_Diag(idx)%mod_name       = 'land'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'total column soil moisture content [kg/m**2] at surface'
    IPD_Diag(idx)%unit           = 'kg/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'soil layer'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_th
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%soilm
    enddo

    !---tmpmin
    idx = idx + 1
    IPD_Diag(idx)%name           = 'tmpmin'
    IPD_Diag(idx)%output_name    = 'tmin'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'minimum temperature at 2 m above ground'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'min'
    IPD_Diag(idx)%level_type     = '2 m above grnd'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%tmpmin
    enddo

    !---tmpmax
    idx = idx + 1
    IPD_Diag(idx)%name           = 'tmpmax'
    IPD_Diag(idx)%output_name    = 'tmax'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'maximum temperature at 2 m above ground'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'max'
    IPD_Diag(idx)%level_type     = '2 m above grnd'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%tmpmax
    enddo

    !---dusfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dusfc'
    IPD_Diag(idx)%output_name    = 'uflx'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Zonal momentum flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dusfc
    enddo

    !---dvsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dvsfc'
    IPD_Diag(idx)%output_name    = 'vflx'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'meridional momentum flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dvsfc
    enddo

    !---dtsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dtsfc'
    IPD_Diag(idx)%output_name    = 'shtfl'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'surface sensible heat flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dtsfc
    enddo

    !---dqsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dqsfc'
    IPD_Diag(idx)%output_name    = 'lhtfl '
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'surface latent heat flux [W/m**2]'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dqsfc
    enddo

    !---totprcp
    idx = idx + 1
    IPD_Diag(idx)%name           = 'totprcp'
    IPD_Diag(idx)%output_name    = 'prate'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'precipitation rate [kg/m**2/s] at surface'
    IPD_Diag(idx)%unit           = 'kg/m**2/s'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%totprcp
    enddo

    !---gflux
    idx = idx + 1
    IPD_Diag(idx)%name           = 'gflux'
    IPD_Diag(idx)%output_name    = 'gflux'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'ground heat flux [W/m**2/s] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%gflux
    enddo

    !---dlwsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dlwsfc'
    IPD_Diag(idx)%output_name    = 'dlwrf'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Time accumulated downward long wave radiation flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dlwsfc
    enddo

    !---ulwsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'ulwsfc'
    IPD_Diag(idx)%output_name    = 'ulwrf'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'time accumulated upward lw flux at surface [W/m**2]'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'accumulation'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%ulwsfc
    enddo

    !---suntim
    idx = idx + 1
    IPD_Diag(idx)%name           = 'suntim'
    IPD_Diag(idx)%output_name    = 'sunsd'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Accumulated sunshine duration time [s]'
    IPD_Diag(idx)%unit           = 's'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%suntim
    enddo

    !---runoff
    idx = idx + 1
    IPD_Diag(idx)%name           = 'runoff'
    IPD_Diag(idx)%output_name    = 'RUNOFF'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'total water runoff [kg/m**2] at surface'
    IPD_Diag(idx)%unit           = 'kg/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_th
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%runoff
    enddo

    !---ep
    idx = idx + 1
    IPD_Diag(idx)%name           = 'ep'
    IPD_Diag(idx)%output_name    = 'pevpr'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Potential evaporation rate [w/m**2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%ep
    enddo

    !---cldwrk
    idx = idx + 1
    IPD_Diag(idx)%name           = 'cldwrk'
    IPD_Diag(idx)%output_name    = 'cwork'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'cloud work function (valid only with SAS) [J/kg] at total atmospheric column'
    IPD_Diag(idx)%unit           = 'J/kg'
    IPD_Diag(idx)%type_stat_proc = 'acc '
    IPD_Diag(idx)%level_type     = 'rntire_atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%cldwrk
    enddo

    !---dugwd
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dugwd'
    IPD_Diag(idx)%output_name    = 'u-gwd'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Zonal gravity wave stress [N/m**2] at surface'
    IPD_Diag(idx)%unit           = 'N/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dugwd
    enddo

    !---dvgwd
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dvgwd'
    IPD_Diag(idx)%output_name    = 'v-gwd'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Meridional gravity wave stress [N/m**2] at surface'
    IPD_Diag(idx)%unit           = 'N/m**2'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dvgwd
    enddo

    !---psmean
    idx = idx + 1
    IPD_Diag(idx)%name           = 'psmean'
    IPD_Diag(idx)%output_name    = 'pres'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'time accumulation surface pressure [kPa*s]'
    IPD_Diag(idx)%unit           = 'kPa*s'
    IPD_Diag(idx)%type_stat_proc = 'acc'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%psmean
    enddo

    !---cnvprcp
    idx = idx + 1
    IPD_Diag(idx)%name           = 'cnvprcp'
    IPD_Diag(idx)%output_name    = 'cprat'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'convective precipitation rate [kg/m**2/s] at surface'
    IPD_Diag(idx)%unit           = 'kg/m**2/s'
    IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%cnvprcp
    enddo

    !---spfhmin
    idx = idx + 1
    IPD_Diag(idx)%name           = 'spfhmin'
    IPD_Diag(idx)%output_name    = 'spfhmin'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'minimum specific humidity [kg/kg] at 2 m above ground'
    IPD_Diag(idx)%unit           = 'kg/kg '
    IPD_Diag(idx)%type_stat_proc = 'min'
    IPD_Diag(idx)%level_type     = '2 m above grnd'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%spfhmin
    enddo

    !---spfhmax
    idx = idx + 1
    IPD_Diag(idx)%name           = 'spfhmaxn'
    IPD_Diag(idx)%output_name    = 'spfhmax'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'maximum specific humidity [kg/kg] at 2 m above ground'
    IPD_Diag(idx)%unit           = 'kg/kg'
    IPD_Diag(idx)%type_stat_proc = 'max'
    IPD_Diag(idx)%level_type     = '2 m above grnd'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%spfhmax
    enddo

    !---rain
    idx = idx + 1
    IPD_Diag(idx)%name           = 'rain'
    IPD_Diag(idx)%output_name    = 'APCP'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'instantaneous total precipitation [m] at surface'
    IPD_Diag(idx)%unit           = 'm'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%rain
    enddo

    !---rainc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'rainc'
    IPD_Diag(idx)%output_name    = 'ACPCP'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'instantaneous convective precipitation [m] at surface'
    IPD_Diag(idx)%unit           = 'm'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%rainc
    enddo

    !---ice
    idx = idx + 1
    IPD_Diag(idx)%name           = 'ice'
    IPD_Diag(idx)%output_name    = 'ICE'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'instantaneous ice fall'
    IPD_Diag(idx)%unit           = ' '
    IPD_Diag(idx)%type_stat_proc = ' '
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%ice
    enddo

    !---snow
    idx = idx + 1
    IPD_Diag(idx)%name           = 'snow'
    IPD_Diag(idx)%output_name    = 'SNOW'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'instantaneous snow fall'
    IPD_Diag(idx)%unit           = ' '
    IPD_Diag(idx)%type_stat_proc = ' '
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%snow
    enddo

    !---graupel
    idx = idx + 1
    IPD_Diag(idx)%name           = 'graupel'
    IPD_Diag(idx)%output_name    = 'GRAUPEL'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'instantaneous graupel fall'
    IPD_Diag(idx)%unit           = ' '
    IPD_Diag(idx)%type_stat_proc = ' '
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%graupel
    enddo

    !---totice
    idx = idx + 1
    IPD_Diag(idx)%name           = 'totice'
    IPD_Diag(idx)%output_name    = 'TOTICE'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'surface ice precipitation rate [kg/m**2/s]'
    IPD_Diag(idx)%unit           = 'kg/m**2/s'
    IPD_Diag(idx)%type_stat_proc = 'accumulation'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%totice
    enddo

    !---totsnw
    idx = idx + 1
    IPD_Diag(idx)%name           = 'totsnw'
    IPD_Diag(idx)%output_name    = 'TOTSNW'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'surface snow precipitation rate [kg/m**2/s]'
    IPD_Diag(idx)%unit           = 'kg/m**2/s'
    IPD_Diag(idx)%type_stat_proc = 'accumulation'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%totsnw
    enddo

    !---totgrp
    idx = idx + 1
    IPD_Diag(idx)%name           = 'totgrp'
    IPD_Diag(idx)%output_name    = 'TOTGRP'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'surface graupel precipitation rate [kg/m**2/s]'
    IPD_Diag(idx)%unit           = 'kg/m**2/s'
    IPD_Diag(idx)%type_stat_proc = 'accumulation'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one/cn_hr
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%totgrp
    enddo

    !---u10m
    idx = idx + 1
    IPD_Diag(idx)%name           = 'u10m'
    IPD_Diag(idx)%output_name    = 'ugrd'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'u wind component [m/s] at 10 m abover ground'
    IPD_Diag(idx)%unit           = 'm/s'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '10 m above grnd'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%u10m
    enddo

    !---v10m
    idx = idx + 1
    IPD_Diag(idx)%name           = 'v10m'
    IPD_Diag(idx)%output_name    = 'vgrd'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'v wind component [m/s] at 10 m above ground'
    IPD_Diag(idx)%unit           = 'm/s'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '10 m above ground'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%v10m
    enddo

    !---zlvl
    idx = idx + 1
    IPD_Diag(idx)%name           = 'zlvl'
    IPD_Diag(idx)%output_name    = 'hgt'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'height [m] at model layer 1'
    IPD_Diag(idx)%unit           = 'm'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'model layer'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%zlvl
    enddo

    !---psurf
    idx = idx + 1
    IPD_Diag(idx)%name           = 'psurf'
    IPD_Diag(idx)%output_name    = 'pres'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'surface pressure [Pa]'
    IPD_Diag(idx)%unit           = 'Pa'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%psurf
    enddo

    !---hpbl
    idx = idx + 1
    IPD_Diag(idx)%name           = 'hpbl'
    IPD_Diag(idx)%output_name    = 'HPBL'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'surface planetary boundary layer height [m]'
    IPD_Diag(idx)%unit           = 'm'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%hpbl
    enddo

    !---pwat
    idx = idx + 1
    IPD_Diag(idx)%name           = 'pwat'
    IPD_Diag(idx)%output_name    = 'pwat'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'atmos column precipitable water [kg/m**2]'
    IPD_Diag(idx)%unit           = 'kg/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'entire atmos'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%pwat
    enddo

    !---t1
    idx = idx + 1
    IPD_Diag(idx)%name           = 't1'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Temperature [K] at model layer 1'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'model lyaer'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%t1
    enddo

    !---q1
    idx = idx + 1
    IPD_Diag(idx)%name           = 'q1'
    IPD_Diag(idx)%output_name    = 'spfh'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'specific humidity [kg/kg] sy model layer 1'
    IPD_Diag(idx)%unit           = 'kg/kg'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'model layer'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%q1
    enddo

    !---u1
    idx = idx + 1
    IPD_Diag(idx)%name           = 'u1'
    IPD_Diag(idx)%output_name    = 'ugrd'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'layer 1 zonal wind [m/s] at model ayer 1'
    IPD_Diag(idx)%unit           = 'm/s'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'model layer'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%u1
    enddo

    !---v1
    idx = idx + 1
    IPD_Diag(idx)%name           = 'v1'
    IPD_Diag(idx)%output_name    = 'vgrd'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'meridional wind [m/s] at model ayer 1'
    IPD_Diag(idx)%unit           = 'm/s'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'model layer'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%v1
    enddo

    !---chh
    idx = idx + 1
    IPD_Diag(idx)%name           = 'chh'
    IPD_Diag(idx)%output_name    = 'CHH'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'thermal exchange coefficient'
    IPD_Diag(idx)%unit           = 'kg/m**2/s'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%chh
    enddo

    !---cmm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'cmm'
    IPD_Diag(idx)%output_name    = 'CMM'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'momentum exchange coefficient'
    IPD_Diag(idx)%unit           = 'm/s'
    IPD_Diag(idx)%type_stat_proc = ' '
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%cmm
    enddo

    !---dlwsfci
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dlwsfci'
    IPD_Diag(idx)%output_name    = 'dlwrf'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'instantaneous downward lw flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dlwsfci
    enddo

    !---ulwsfci
    idx = idx + 1
    IPD_Diag(idx)%name           = 'ulwsfci'
    IPD_Diag(idx)%output_name    = 'ulwrf'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'instantaneous upward lw flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%ulwsfci
    enddo

    !---dswsfci
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dswsfci'
    IPD_Diag(idx)%output_name    = 'dswrf'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'instantaneous downward sw flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dswsfci
    enddo

    !---uswsfci
    idx = idx + 1
    IPD_Diag(idx)%name           = 'uswsfci'
    IPD_Diag(idx)%output_name    = 'uswrf'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'instantaneous upward sw flux [W/m**2] at surface'
    IPD_Diag(idx)%unit           = 'W/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%uswsfci
    enddo

    !---dusfci
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dusfci'
    IPD_Diag(idx)%output_name    = 'uflx'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'instantaneous Zonal compt of momentum flux at surface'
    IPD_Diag(idx)%unit           = 'n/m**2 '
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dusfci
    enddo

    !---dvsfci
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dvsfci'
    IPD_Diag(idx)%output_name    = 'vflx'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'instantaneous v component of surface stress'
    IPD_Diag(idx)%unit           = 'n/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dvsfci
    enddo

    !---dtsfci
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dtsfci'
    IPD_Diag(idx)%output_name    = 'shtfl'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'instantaneous surface sensible heat flux [w/m**2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dtsfci
    enddo

    !---dqsfci
    idx = idx + 1
    IPD_Diag(idx)%name           = 'dqsfci'
    IPD_Diag(idx)%output_name    = 'lhtfl'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'instantaneous latent heat flux [w/m**2] at surface'
    IPD_Diag(idx)%unit           = 'w/m**2 '
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%dqsfci
    enddo

    !---gfluxi
    idx = idx + 1
    IPD_Diag(idx)%name           = 'gfluxi'
    IPD_Diag(idx)%output_name    = 'gflux'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'instantaneous ground heat flux [w/m**2]at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%gfluxi
    enddo

    !---epi
    idx = idx + 1
    IPD_Diag(idx)%name           = 'epi'
    IPD_Diag(idx)%output_name    = 'pevpr'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'cpl'
    IPD_Diag(idx)%desc           = 'instantaneous potential evaporation rate (w/m**2) at surface'
    IPD_Diag(idx)%unit           = 'w/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc '
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%epi
    enddo

    !---smcwlt2
    idx = idx + 1
    IPD_Diag(idx)%name           = 'smcwlt2'
    IPD_Diag(idx)%output_name    = 'wilt'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'wilting point (proportion) at surface'
    IPD_Diag(idx)%unit           = 'proportion '
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%smcwlt2
    enddo

    !---smcref2'
    idx = idx + 1
    IPD_Diag(idx)%name           = 'smcref2'
    IPD_Diag(idx)%output_name    = 'SMCREF2'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Field capacity (fraction) at surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%smcref2
    enddo

    !---wet1
    idx = idx + 1
    IPD_Diag(idx)%name           = 'wet1'
    IPD_Diag(idx)%output_name    = 'WET1'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'gocart_cpl '
    IPD_Diag(idx)%desc           = 'normalized top soil layer wetness [frqaction]'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'soil layer'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%wet1
    enddo

    !---sr
    idx = idx + 1
    IPD_Diag(idx)%name           = 'sr'
    IPD_Diag(idx)%output_name    = 'cpofp'
    IPD_Diag(idx)%mod_name       = 'physics'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Frozen precipitation fraction'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Diag(nb)%sr
    enddo

    if (Model%ldiag3d) then
      !---dt3dt
      idx = idx + 1
      IPD_Diag(idx)%name           = 'dt3dt1'
      IPD_Diag(idx)%output_name    = 'LWHR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Long wave radiative heating rate [K/s] at model layers '
      IPD_Diag(idx)%unit           = 'K/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dt3dt(:,:,1)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dt3dt2'
      IPD_Diag(idx)%output_name    = 'SWHR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Short wave radiative heating rate [K/s] at model layers '
      IPD_Diag(idx)%unit           = 'K/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dt3dt(:,:,2)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dt3dt3'
      IPD_Diag(idx)%output_name    = 'VDFHR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Vertical diffusion heating rate [K/s] at model layers '
      IPD_Diag(idx)%unit           = 'K/s'
      IPD_Diag(idx)%type_stat_proc = 'amm_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dt3dt(:,:,3)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dt3dt4'
      IPD_Diag(idx)%output_name    = 'CNVHR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Deep convective heating rate [K/s] at model layers '
      IPD_Diag(idx)%unit           = 'K/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dt3dt(:,:,4)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dt3dt5'
      IPD_Diag(idx)%output_name    = 'SHAHR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Shallow convective heating rate [K/s] at model layers '
      IPD_Diag(idx)%unit           = 'K/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dt3dt(:,:,5)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dt3dt6'
      IPD_Diag(idx)%output_name    = 'LRGHR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Large scale condensate heat rate [K/s] at model layers '
      IPD_Diag(idx)%unit           = 'K/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dt3dt(:,:,6)
      enddo
      
      !---dq3dt
      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt1'
      IPD_Diag(idx)%output_name    = 'VDFMR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Vertical diffusion moistening rate [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,1)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt2'
      IPD_Diag(idx)%output_name    = 'CNVMR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Deep convective moistening rate [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,2)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt3'
      IPD_Diag(idx)%output_name    = 'SHAMR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Shallow convective moistening rate [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,3)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt4'
      IPD_Diag(idx)%output_name    = 'LRGMR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Large scale moistening rate [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,4)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt5'
      IPD_Diag(idx)%output_name    = 'VDFOZ'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Ozone vertical diffusion [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,5)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt6'
      IPD_Diag(idx)%output_name    = 'POZ'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Ozone production [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,6)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt7'
      IPD_Diag(idx)%output_name    = 'TOZ'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Ozone tendency [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,7)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt8'
      IPD_Diag(idx)%output_name    = 'POZT'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Ozone production from temperature term [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,8)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dq3dt9'
      IPD_Diag(idx)%output_name    = 'POZO'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Ozone production from col ozone term [kg/kg/s] at model layers'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dq3dt(:,:,9)
      enddo

      !---du3dt
      idx = idx + 1
      IPD_Diag(idx)%name           = 'du3dt1'
      IPD_Diag(idx)%output_name    = 'VDFUA'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Vertical diffusion zonal acceleration [m/s**2] at model layers'
      IPD_Diag(idx)%unit           = 'm/s**2'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%du3dt(:,:,1)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'du3dt2'
      IPD_Diag(idx)%output_name    = 'GWDU'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Gravity wave drag zonal acceleration [m/s**2] at model layers'
      IPD_Diag(idx)%unit           = 'm/s**2'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%du3dt(:,:,2)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'du3dt3'
      IPD_Diag(idx)%output_name    = 'CNVU'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Convective zonal momentum mixing acceleration [m/s**2] at model layers'
      IPD_Diag(idx)%unit           = 'm/s**2'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%du3dt(:,:,3)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'du3dt4'
      IPD_Diag(idx)%output_name    = 'CNGWDU'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Convective Gravity wave drag zonal momentum mixing acceleration [m/s**2] at model layers'
      IPD_Diag(idx)%unit           = 'm/s**2'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%du3dt(:,:,4)
      enddo

      !---dv3dt
      idx = idx + 1
      IPD_Diag(idx)%name           = 'dv3dt1'
      IPD_Diag(idx)%output_name    = 'VDFVA'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Vertical diffusion meridional acceleration [m/s**2] at model layers'
      IPD_Diag(idx)%unit           = 'm/s**2'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dv3dt(:,:,1)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dv3dt2'
      IPD_Diag(idx)%output_name    = 'GWDV'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Gravity wave drag meridional acceleration [m/s**2] at model layers'
      IPD_Diag(idx)%unit           = 'm/s**2'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dv3dt(:,:,2)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dv3dt3'
      IPD_Diag(idx)%output_name    = 'CNVV'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Convective meridional momentum mixing acceleration [m/s**2] at model layers'
      IPD_Diag(idx)%unit           = 'm/s**2'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dv3dt(:,:,3)
      enddo

      idx = idx + 1
      IPD_Diag(idx)%name           = 'dv3dt4'
      IPD_Diag(idx)%output_name    = 'CNGWDV'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Convective Gravity wave drag meridional acceleration [m/s**2] at model layers'
      IPD_Diag(idx)%unit           = 'm/s**2'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtime'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%dv3dt(:,:,4)
      enddo
    
      idx = idx + 1
      IPD_Diag(idx)%name           = 'cldcov'
      IPD_Diag(idx)%output_name    = 'CDLYR'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'diag3d'
      IPD_Diag(idx)%desc           = 'Amount of non-convective cloud [%] at model layers'
      IPD_Diag(idx)%unit           = '%'
      IPD_Diag(idx)%type_stat_proc = 'acc_rtimsw'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      IPD_Diag(idx)%fcst_hour      = Model%fhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Diag(nb)%cldcov
      enddo
    endif

    if (Model%lgocart) then
      !---dqdtv
      idx = idx + 1
      IPD_Diag(idx)%name           = 'dqdti'
      IPD_Diag(idx)%output_name    = 'DQDTV'
      IPD_Diag(idx)%mod_name       = 'physics'
      IPD_Diag(idx)%file_name      = 'gocart_cpl'
      IPD_Diag(idx)%desc           = 'instantaneous total moisture tendency [kg/kg/s]'
      IPD_Diag(idx)%unit           = 'kg/kg/s'
      IPD_Diag(idx)%type_stat_proc = 'inst'
      IPD_Diag(idx)%level_type     = 'model layer'
      IPD_Diag(idx)%level          = 64
      IPD_Diag(idx)%cnvfac         = cn_one
      IPD_Diag(idx)%zhour          = Model%zhour
      do nb = 1,nblks
        IPD_Diag(idx)%data(nb)%var3p => Coupling(nb)%dqdti
      enddo
    endif

    ! GFS_SFCPROP CONTAINER DATA: with ialb=0: climatological albedo scheme
    !---alnsf 
    idx = idx + 1
    IPD_Diag(idx)%name           = 'alnsf'
    IPD_Diag(idx)%output_name    = 'alnsf'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'mean nir albedo with strong cosz dependency'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = ''
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%alnsf
    enddo

    !---alnwf
    idx = idx + 1
    IPD_Diag(idx)%name           = 'alnwf'
    IPD_Diag(idx)%output_name    = 'ALNWF'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'nir albedo with weak cosz dependency [%] at surface'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%alnwf
    enddo

    !---alvsf
    idx = idx + 1
    IPD_Diag(idx)%name           = 'alvsf'
    IPD_Diag(idx)%output_name    = 'ALVSF'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'uv+vis albedo with strong cosz dependency [%] at surface'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%alvsf
    enddo

    !---alvwf
    idx = idx + 1
    IPD_Diag(idx)%name           = 'alvwf'
    IPD_Diag(idx)%output_name    = 'ALVWF'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'vis albedo with weak cosz dependency'
    IPD_Diag(idx)%unit           = '%'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%alvwf
    enddo

    !---canopy
    idx = idx + 1
    IPD_Diag(idx)%name           = 'canopy'
    IPD_Diag(idx)%output_name    = 'CNWAT'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Canopy water content (kg/m**2)'
    IPD_Diag(idx)%unit           = 'kg/m**2 '
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%canopy
    enddo

    !---f10m
    idx = idx + 1
    IPD_Diag(idx)%name           = 'f10m'
    IPD_Diag(idx)%output_name    = 'F10M'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'ratio of 10-meter wind speed to the lowest model layer wind speed [numeric]'
    IPD_Diag(idx)%unit           = 'numeric'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '10 m above grnd'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%f10m
    enddo

    !---facsf
    idx = idx + 1
    IPD_Diag(idx)%name           = 'facsf'
    IPD_Diag(idx)%output_name    = 'FACSF'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'fractional coverage with strong cosz dependency'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%facsf
    enddo

    !---facwf
    idx = idx + 1
    IPD_Diag(idx)%name           = 'facwf'
    IPD_Diag(idx)%output_name    = 'FACWF'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'fractional coverage with weak cosz dependency [fraction] '
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%facwf
    enddo

    !---ffhh
    idx = idx + 1
    IPD_Diag(idx)%name           = 'ffhh'
    IPD_Diag(idx)%output_name    = 'FFHH'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'stability profile function [numeric] for heat at surface'
    IPD_Diag(idx)%unit           = 'numeric'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%ffhh
    enddo

    !---ffmm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'ffmm'
    IPD_Diag(idx)%output_name    = 'FFMM'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'stability profile function [numeric] for momentum at surface layer'
    IPD_Diag(idx)%unit           = 'numeric'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%ffmm
    enddo

    !---fice
    idx = idx + 1
    IPD_Diag(idx)%name           = 'fice'
    IPD_Diag(idx)%output_name    = 'icec'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'surface ice concentration (ice=1, no ice=0) [fraction]'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%fice
    enddo

    !---hice
    idx = idx + 1
    IPD_Diag(idx)%name           = 'hice'
    IPD_Diag(idx)%output_name    = 'icetk'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'ice thickness [m] at surface'
    IPD_Diag(idx)%unit           = 'm'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%hice
    enddo

    !---snoalb
    idx = idx + 1
    IPD_Diag(idx)%name           = 'snoalb'
    IPD_Diag(idx)%output_name    = 'mxsalb'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'maximum snow albedo in fraction'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%snoalb
    enddo

    !---shdmax
    idx = idx + 1
    IPD_Diag(idx)%name           = 'shdmax'
    IPD_Diag(idx)%output_name    = 'SHDMAX'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'maximum fractional coverage of green vegetation'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'max'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%shdmax
    enddo

    !---shdmin
    idx = idx + 1
    IPD_Diag(idx)%name           = 'shdmin'
    IPD_Diag(idx)%output_name    = 'SHDMIN'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'minimum fractional coverage of green vegetation'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'min'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%shdmin
    enddo

    !---snowd
    idx = idx + 1
    IPD_Diag(idx)%name           = 'snowd'
    IPD_Diag(idx)%output_name    = 'SNOWD'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'surface snow depth [m]'
    IPD_Diag(idx)%unit           = 'm'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%snowd
    enddo

    !---stype
    idx = idx + 1
    IPD_Diag(idx)%name           = 'stype'
    IPD_Diag(idx)%output_name    = 'sotyp'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'soil type at surface'
    IPD_Diag(idx)%unit           = 'numeric'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%stype
    enddo

    !---q2m
    idx = idx + 1
    IPD_Diag(idx)%name           = 'q2m'
    IPD_Diag(idx)%output_name    = 'spfh'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'specific humidity [kg/kg] at 2 m above ground'
    IPD_Diag(idx)%unit           = 'kg/kg'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '2 m above grnd'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%q2m
    enddo

    !---t2m
    idx = idx + 1
    IPD_Diag(idx)%name           = 't2m'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'temperature [K] at 2 m above ground'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '2 m above grnd'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%t2m
    enddo

    !---tsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'tsfc'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'surface temperature [K]'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%tsfc
    enddo

    !---tg3
    idx = idx + 1
    IPD_Diag(idx)%name           = 'tg3'
    IPD_Diag(idx)%output_name    = 'TG3'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'deep soil temperature [K]'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%tg3
    enddo

    !---tisfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'tisfc'
    IPD_Diag(idx)%output_name    = 'TISFC'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'surface temperature over ice fraction [K]'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%tisfc
    enddo

    !---tprcp
    idx = idx + 1
    IPD_Diag(idx)%name           = 'tprcp'
    IPD_Diag(idx)%output_name    = 'TPRCP'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'total precipitation at surface [m]'
    IPD_Diag(idx)%unit           = 'm'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%tprcp
    enddo

    !---vtype
    idx = idx + 1
    IPD_Diag(idx)%name           = 'vtype'
    IPD_Diag(idx)%output_name    = 'VGTYP'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'vegetation type [index]'
    IPD_Diag(idx)%unit           = 'index'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%vtype
    enddo

    !---weasd
    idx = idx + 1
    IPD_Diag(idx)%name           = 'weasd'
    IPD_Diag(idx)%output_name    = 'WEASD'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Water equivalent of accumulated snow depth [kg/m**2] at surface'
    IPD_Diag(idx)%unit           = 'kg/m**2'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%weasd
    enddo

    !---hgtsfc
    idx = idx + 1
    IPD_Diag(idx)%name           = 'hgtsfc'
    IPD_Diag(idx)%output_name    = 'HGT'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'surface geopotential height [gpm]'
    IPD_Diag(idx)%unit           = 'gpm'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%oro
    enddo

    !---slmsk
    idx = idx + 1
    IPD_Diag(idx)%name           = 'slmsk'
    IPD_Diag(idx)%output_name    = 'LAND'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'sea-land-ice mask (0-sea, 1-land, 2-ice)'
    IPD_Diag(idx)%unit           = 'numeric'
    IPD_Diag(idx)%type_stat_proc = ''
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%slmsk
    enddo

    !---zorl
    idx = idx + 1
    IPD_Diag(idx)%name           = 'zorl'
    IPD_Diag(idx)%output_name    = 'sfcr'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'surface roughness [m]'
    IPD_Diag(idx)%unit           = 'm'
    IPD_Diag(idx)%type_stat_proc = ' '
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%zorl
    enddo

    !---vfrac
    idx = idx + 1
    IPD_Diag(idx)%name           = 'vfrac'
    IPD_Diag(idx)%output_name    = 'veg'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'sfc'
    IPD_Diag(idx)%desc           = 'vegetation fraction'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = 'sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%vfrac
    enddo

    !---slc 0-10cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'slc'
    IPD_Diag(idx)%output_name    = 'soill'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'liquid soil moisture content at 0-10cm below land surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '0-10cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%slc(:,1)
    enddo

    !---slc 10-40cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'slc'
    IPD_Diag(idx)%output_name    = 'soill'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'liquid soil moisture content at 10-40cm below land surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '10-40cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%slc(:,2)
    enddo

    !---slc 40-100cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'slc'
    IPD_Diag(idx)%output_name    = 'SLC'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'liquid soil moisture content at 40-100cm below land surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '40-100cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%slc(:,3)
    enddo

    !---slc 100-200cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'slc'
    IPD_Diag(idx)%output_name    = 'soill'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'liquid soil moisture content at 100-200cm below land surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '100-200cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%slc(:,4)
    enddo

    !---smc 0-10cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'smc'
    IPD_Diag(idx)%output_name    = 'soilw'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Volumetric soil moist content (frac) at 0-10cm below land surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '0-10cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%smc(:,1)
    enddo

    !---smc 10-40cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'smc'
    IPD_Diag(idx)%output_name    = 'soilw'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Volumetric soil moist content (frac) at 10-40cm below land surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '10-40cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%smc(:,2)
    enddo

    !---smc 40-100cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'smc'
    IPD_Diag(idx)%output_name    = 'soilw'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Volumetric soil moist content (frac) at 40-100cm below land surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '40-100cm'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%smc(:,3)
    enddo

    !---smc 100-200cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'smc'
    IPD_Diag(idx)%output_name    = 'soilw'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'Volumetric soil moist content (frac) at 100-200cm below land surface'
    IPD_Diag(idx)%unit           = 'fraction'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '100-200cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%smc(:,4)
    enddo

    !---stc 0-10cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'stc'
    IPD_Diag(idx)%output_name    = 'TMP'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = ' '
    IPD_Diag(idx)%desc           = 'soil temperature at 0-10cm below land surface'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '0-10cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%stc(:,1)
    enddo

    !---stc 10-40cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'stc'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'soil temperature at 10-40cm below land surface'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '10-40cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%stc(:,2)
    enddo

    !---stc 40-100cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'stc'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'soil temperature at 40-100cm below land surface'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '40-100cm below land surface'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%stc(:,3)
    enddo

    !---stc 100-200cm
    idx = idx + 1
    IPD_Diag(idx)%name           = 'stc'
    IPD_Diag(idx)%output_name    = 'tmp'
    IPD_Diag(idx)%mod_name       = 'surface'
    IPD_Diag(idx)%file_name      = 'flx'
    IPD_Diag(idx)%desc           = 'soil temperature at 100-200cm below land surface'
    IPD_Diag(idx)%unit           = 'K'
    IPD_Diag(idx)%type_stat_proc = 'inst'
    IPD_Diag(idx)%level_type     = '100-200cm below sfc'
    IPD_Diag(idx)%level          = 1
    IPD_Diag(idx)%cnvfac         = cn_one
    IPD_Diag(idx)%zhour          = Model%zhour
    IPD_Diag(idx)%fcst_hour      = Model%fhour
    do nb = 1,nblks
      IPD_Diag(idx)%data(nb)%var2p => Sfcprop(nb)%stc(:,4)
    enddo

    if (idx > size(IPD_Diag)) then
      print *,'GFS_populate_IPD_Diag: increase size declaration of IPD_Diag'
      stop
    endif

  end subroutine diag_populate

end module physics_diag_layer

