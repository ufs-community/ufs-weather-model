module GFS_initialize_scm

  implicit none

  private

!----------------
! Public entities
!----------------
  public  GFS_initialize_scm_init, GFS_initialize_scm_run, GFS_initialize_scm_finalize

  CONTAINS
!*******************************************************************************************

!--------------
! GFS initialze
!--------------

  subroutine GFS_initialize_scm_init()
  end subroutine GFS_initialize_scm_init

  subroutine GFS_initialize_scm_finalize()
  end subroutine GFS_initialize_scm_finalize

!> \section arg_table_GFS_initialize_scm_run Argument Table
!! | local_name           | standard_name                                               | long_name                                                               | units         | rank | type                          |    kind   | intent | optional |
!! |----------------------|-------------------------------------------------------------|-------------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Model                | FV3-GFS_Control_type                                        | Fortran DDT containing FV3-GFS model control parameters                 | DDT           |    0 | GFS_control_type              |           | inout  | F        |
!! | Statein              | FV3-GFS_Statein_type                                        | Fortran DDT containing FV3-GFS prognostic state data in from dycore     | DDT           |    0 | GFS_statein_type              |           | inout  | F        |
!! | Stateout             | FV3-GFS_Stateout_type                                       | Fortran DDT containing FV3-GFS prognostic state to return to dycore     | DDT           |    0 | GFS_stateout_type             |           | inout  | F        |
!! | Sfcprop              | FV3-GFS_Sfcprop_type                                        | Fortran DDT containing FV3-GFS surface fields                           | DDT           |    0 | GFS_sfcprop_type              |           | inout  | F        |
!! | Coupling             | FV3-GFS_Coupling_type                                       | derived type GFS_coupling_type in FV3                                   | DDT           |    0 | GFS_coupling_type             |           | inout  | F        |
!! | Grid                 | FV3-GFS_Grid_type                                           | Fortran DDT containing FV3-GFS grid and interpolation related data      | DDT           |    0 | GFS_grid_type                 |           | inout  | F        |
!! | Tbd                  | FV3-GFS_Tbd_type                                            | derived type GFS_tbd_type in FV3                                        | DDT           |    0 | GFS_tbd_type                  |           | inout  | F        |
!! | Cldprop              | FV3-GFS_Cldprop_type                                        | derived type GFS_cldprop_type in FV3                                    | DDT           |    0 | GFS_cldprop_type              |           | inout  | F        |
!! | Radtend              | FV3-GFS_Radtend_type                                        | derived type GFS_radtend_type in FV3                                    | DDT           |    0 | GFS_radtend_type              |           | inout  | F        |
!! | Diag                 | FV3-GFS_Diag_type                                           | Fortran DDT containing FV3-GFS fields targeted for diagnostic output    | DDT           |    0 | GFS_diag_type                 |           | inout  | F        |
!! | Interstitial         | FV3-GFS_Interstitial_type                                   | derived type GFS_interstitial_type in FV3                               | DDT           |    0 | GFS_interstitial_type         |           | inout  | F        |
!! | Init_parm            | FV3-GFS_Init_type                                           | dervied type GFS_init_type in FV3                                       | DDT           |    0 | GFS_init_type                 |           | in     | F        |
!! | n_ozone_layers       | vertical_dimension_of_ozone_forcing_data_from_host          | number of vertical layers in ozone forcing data coming from host        | count         |    0 | integer                       |           | in     | F        |
!! | n_ozone_lats         | number_of_latitutde_points_in_ozone_forcing_data_from_host  | number of latitude points in ozone forcing data coming from host        | count         |    0 | integer                       |           | in     | F        |
!! | n_ozone_times        | number_of_time_levels_in_ozone_forcing_data_from_host       | number of time levels in ozone forcing data coming from host            | count         |    0 | integer                       |           | in     | F        |
!! | n_ozone_coefficients | number_of_coefficients_in_ozone_forcing_data_from_host      | number of coeffcients in ozone forcing data coming from host            | count         |    0 | integer                       |           | in     | F        |
!! | ozone_lat            | latitude_of_ozone_forcing_data_from_host                    | latitude value of the ozone forcing data coming from host               | degree        |    1 | real                          | kind_phys | in     | F        |
!! | ozone_pres           | natural_log_of_ozone_forcing_data_pressure_levels_from_host | natural logarithm of the pressure levels of the ozone forcing data      | Pa            |    1 | real                          | kind_phys | in     | F        |
!! | ozone_time           | time_levels_in_ozone_forcing_data_from_host                 | time values of the ozone forcing data coming from host                  | day           |    1 | real                          | kind_phys | in     | F        |
!! | ozone_forcing_in     | ozone_forcing_from_host                                     | ozone forcing data from host                                            | various       |    4 | real                          | kind_phys | in     | F        |
!! | errmsg               | error_message                                               | error message for error handling in CCPP                                | none          |    0 | character                     | len=*     | out    | F        |
!! | errflg               | error_flag                                                  | error flag for error handling in CCPP                                   | flag          |    0 | integer                       |           | out    | F        |
!!
  subroutine GFS_initialize_scm_run (Model, Statein, Stateout, Sfcprop,           &
                             Coupling, Grid, Tbd, Cldprop, Radtend, Diag,         &
                             Interstitial, Init_parm, n_ozone_lats,               &
                             n_ozone_layers, n_ozone_times, n_ozone_coefficients, &
                             ozone_lat, ozone_pres, ozone_time, ozone_forcing_in, &
                             errmsg, errflg)

    use machine,             only: kind_phys
    use GFS_typedefs,        only: GFS_init_type,                          &
                                   GFS_statein_type,  GFS_stateout_type,   &
                                   GFS_sfcprop_type,  GFS_coupling_type,   &
                                   GFS_control_type,  GFS_grid_type,       &
                                   GFS_tbd_type,      GFS_cldprop_type,    &
                                   GFS_radtend_type,  GFS_diag_type,       &
                                   GFS_interstitial_type
    use funcphys,            only: gfuncphys
    use module_microphysics, only: gsmconst
    use cldwat2m_micro,      only: ini_micro
    use aer_cloud,           only: aer_cloud_init
    use module_ras,          only: ras_init
    use ozne_def,            only: latsozp, levozp, timeoz, oz_coeff, oz_lat, oz_pres, oz_time, ozplin

    !--- interface variables
    type(GFS_control_type),      intent(inout) :: Model
    type(GFS_statein_type),      intent(inout) :: Statein
    type(GFS_stateout_type),     intent(inout) :: Stateout
    type(GFS_sfcprop_type),      intent(inout) :: Sfcprop
    type(GFS_coupling_type),     intent(inout) :: Coupling
    type(GFS_grid_type),         intent(inout) :: Grid
    type(GFS_tbd_type),          intent(inout) :: Tbd
    type(GFS_cldprop_type),      intent(inout) :: Cldprop
    type(GFS_radtend_type),      intent(inout) :: Radtend
    type(GFS_diag_type),         intent(inout) :: Diag
    type(GFS_interstitial_type), intent(inout) :: Interstitial
    type(GFS_init_type),         intent(in)    :: Init_parm

    integer, intent(in) :: n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times
    real(kind=kind_phys), intent(in) :: ozone_lat(:), ozone_pres(:), ozone_time(:), ozone_forcing_in(:,:,:,:)

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    !
    !     !--- local variables
    real(kind=kind_phys), allocatable :: si(:)
    real(kind=kind_phys), parameter   :: p_ref = 101325.0d0

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

!     !--- set control properties (including namelist read)
    call Model%init (Init_parm%nlunit, Init_parm%fn_nml,           &
                     Init_parm%me, Init_parm%master,               &
                     Init_parm%logunit, Init_parm%isc,             &
                     Init_parm%jsc, Init_parm%nx, Init_parm%ny,    &
                     Init_parm%levs, Init_parm%cnx, Init_parm%cny, &
                     Init_parm%gnx, Init_parm%gny,                 &
                     Init_parm%dt_dycore, Init_parm%dt_phys,       &
                     Init_parm%bdat, Init_parm%cdat,               &
                     Init_parm%tracer_names, Init_parm%blksz)

    !allocate memory for the variables stored in ozne_def and set them
    allocate(oz_lat(n_ozone_lats), oz_pres(n_ozone_layers), oz_time(n_ozone_times+1))
    allocate(ozplin(n_ozone_lats, n_ozone_layers, n_ozone_coefficients, n_ozone_times))
    latsozp = n_ozone_lats
    levozp = n_ozone_layers
    timeoz = n_ozone_times
    oz_coeff = n_ozone_coefficients
    oz_lat = ozone_lat
    oz_pres = ozone_pres
    oz_time = ozone_time
    ozplin = ozone_forcing_in

    call Statein%create(1, Model)
    call Stateout%create(1, Model)
    call Sfcprop%create(1, Model)
    call Coupling%create(1, Model)
    call Grid%create(1, Model)
    call Tbd%create(1, 1, Model)
    call Cldprop%create(1, Model)
    call Radtend%create(1, Model)
    !--- internal representation of diagnostics
    call Diag%create(1, Model)
    !--- internal representation of interstitials for CCPP physics
    call Interstitial%create(1, Model)

!     !--- populate the grid components
    call GFS_grid_populate (Grid, Init_parm%xlon, Init_parm%xlat, Init_parm%area)

    !--- read in and initialize ozone and water
    if (Model%ntoz > 0) then
      call setindxoz (Init_parm%blksz, Grid%xlat_d, Grid%jindx1_o3, &
                        Grid%jindx2_o3, Grid%ddy_o3)
    endif

    if (Model%h2o_phys) then
      call setindxh2o (Init_parm%blksz, Grid%xlat_d, Grid%jindx1_h, &
                         Grid%jindx2_h, Grid%ddy_h)
    endif

!     !--- Call gfuncphys (funcphys.f) to compute all physics function tables.
    call gfuncphys ()
!
    call gsmconst (Model%dtp, Model%me, .TRUE.)
!
!     !--- define sigma level for radiation initialization
!     !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
!     !--- ps is replaced with p0. The value of p0 uses that in http://www.emc.ncep.noaa.gov/officenotes/newernotes/on461.pdf
!     !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
    allocate(si(Model%levr+1))
    si = (Init_parm%ak + Init_parm%bk * p_ref - Init_parm%ak(Model%levr+1)) &
             / (p_ref - Init_parm%ak(Model%levr+1))
    call rad_initialize (si, Model%levr, Model%ictm, Model%isol, &
           Model%ico2, Model%iaer, Model%ialb, Model%iems,       &
           Model%ntcw, Model%num_p3d, Model%npdf3d, Model%ntoz,  &
           Model%iovr_sw, Model%iovr_lw, Model%isubc_sw,         &
           Model%isubc_lw, Model%crick_proof, Model%ccnorm,      &
           Model%norad_precip, Model%idate,Model%iflip, Model%me)
    deallocate (si)
!
!     !--- initialize Morrison-Gettleman microphysics
    if (Model%ncld == 2) then
      call ini_micro (Model%mg_dcs, Model%mg_qcvar, Model%mg_ts_auto_ice)
      call aer_cloud_init ()
    endif
!
!     !--- initialize ras
    if (Model%ras) call ras_init (Model%levs, Model%me)
!
!     !--- initialize soil vegetation
    call set_soilveg(Model%me, Model%isot, Model%ivegsrc, Model%nlunit)
!
!     !--- lsidea initialization
    if (Model%lsidea) then
      print *,' LSIDEA is active but needs to be reworked for FV3 - shutting down'
      stop
      !--- NEED TO get the logic from the old phys/gloopb.f initialization area
    endif
!
!     !--- sncovr may not exist in ICs from chgres.
!     !--- FV3GFS handles this as part of the IC ingest
!     !--- this not is placed here to alert users to the need to study
!     !--- the FV3GFS_io.F90 module

  end subroutine GFS_initialize_scm_run

   !------------------
   ! GFS_grid_populate
   !------------------
  subroutine GFS_grid_populate (Grid, xlon, xlat, area)
    use machine,             only: kind_phys
    use physcons,            only: pi => con_pi
    use GFS_typedefs,        only: GFS_grid_type

    implicit none

    type(GFS_grid_type)              :: Grid
    real(kind=kind_phys), intent(in) :: xlon(:,:)
    real(kind=kind_phys), intent(in) :: xlat(:,:)
    real(kind=kind_phys), intent(in) :: area(:,:)

    !--- local variables
    integer :: n_columns, i

    n_columns = size(Grid%xlon)

    do i=1, n_columns
     Grid%xlon = xlon(i,1)
     Grid%xlat = xlat(i,1)
     Grid%xlat_d(i) = xlat(i,1) * 180.0_kind_phys/pi
     Grid%sinlat(i) = sin(Grid%xlat(i))
     Grid%coslat(i) = sqrt(1.0_kind_phys - Grid%sinlat(i)*Grid%sinlat(i))
     Grid%area(i)   = area(i,1)
     Grid%dx(i)     = sqrt(area(i,1))
    end do

  end subroutine GFS_grid_populate

end module GFS_initialize_scm
