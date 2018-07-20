module GFS_driver

  use machine,                  only: kind_phys
  use GFS_typedefs,             only: GFS_init_type,                          &
                                      GFS_statein_type,  GFS_stateout_type,   &
                                      GFS_sfcprop_type,  GFS_coupling_type,   &
                                      GFS_control_type,  GFS_grid_type,       &
                                      GFS_tbd_type,      GFS_cldprop_type,    &
                                      GFS_radtend_type,  GFS_diag_type,       &
                                      GFS_interstitial_type
#ifndef CCPP
  use module_radiation_driver,  only: GFS_radiation_driver
  use module_physics_driver,    only: GFS_physics_driver
#endif
  use module_radsw_parameters,  only: topfsw_type, sfcfsw_type
  use module_radlw_parameters,  only: topflw_type, sfcflw_type
  use funcphys,                 only: gfuncphys

  implicit none

  private

!--------------------------------------------------------------------------------
! GFS_init_type
!--------------------------------------------------------------------------------
!   This container is the minimum set of data required from the dycore/atmosphere
!   component to allow proper initialization of the GFS physics
!
!   Type is defined in GFS_typedefs.F90
!--------------------------------------------------------------------------------
! type GFS_init_type
!   public
!   integer :: me                                !< my MPI-rank
!   integer :: master                            !< master MPI-rank
!   integer :: isc                               !< starting i-index for this MPI-domain
!   integer :: jsc                               !< starting j-index for this MPI-domain
!   integer :: nx                                !< number of points in i-dir for this MPI rank
!   integer :: ny                                !< number of points in j-dir for this MPI rank
!   integer :: levs                              !< number of vertical levels
!   integer :: cnx                               !< number of points in i-dir for this cubed-sphere face
!                                                !< equal to gnx for lat-lon grids
!   integer :: cny                               !< number of points in j-dir for this cubed-sphere face
!                                                !< equal to gny for lat-lon grids
!   integer :: gnx                               !< number of global points in x-dir (i) along the equator
!   integer :: gny                               !< number of global points in y-dir (j) along any meridian
!   integer :: nlunit                            !< fortran unit number for file opens
!   integer :: logunit                           !< fortran unit number for writing logfile
!   integer :: dt_dycore                         !< dynamics time step in seconds
!   integer :: dt_phys                           !< physics  time step in seconds
!   integer :: bdat(8)                           !< model begin date in GFS format   (same as idat)
!   integer :: cdat(8)                           !< model current date in GFS format (same as jdat)
!   !--- blocking data
!   integer, pointer :: blksz(:)                 !< for explicit data blocking
!                                                !< default blksz(1)=[nx*ny]
!   !--- ak/bk for pressure level calculations
!   integer, pointer :: ak(:)                    !< from surface (k=1) to TOA (k=levs)
!   integer, pointer :: bk(:)                    !< from surface (k=1) to TOA (k=levs)
!   !--- grid metrics
!   real(kind=kind_phys), pointer :: xlon(:,:)   !< column longitude for MPI rank
!   real(kind=kind_phys), pointer :: xlat(:,:)   !< column latitude  for MPI rank
!   real(kind=kind_phys), pointer :: area(:,:)   !< column area for length scale calculations
!
!   character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
!                                                 !< based on name location in array
!   character(len=65) :: fn_nml                  !< namelist filename
! end type GFS_init_type
!--------------------------------------------------------------------------------

!------------------
! Module parameters
!------------------

!----------------------------
! Module variable definitions
!----------------------------
  real(kind=kind_phys), parameter :: con_24  =   24.0_kind_phys
  real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
  real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
  real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys
  real(kind=kind_phys), parameter :: qmin    =    1.0e-10

  integer, allocatable :: blksz(:)

!----------------
! Public entities
!----------------
  public  GFS_initialize              !< GFS initialization routine
  public  GFS_time_vary_step          !< perform operations needed prior radiation or physics
#ifndef CCPP
  public  GFS_radiation_driver        !< radiation_driver (was grrad)
  public  GFS_physics_driver          !< physics_driver (was gbphys)
  public  GFS_stochastic_driver       !< stochastic physics
#endif
  public  GFS_finalize

  CONTAINS
!*******************************************************************************************


!---------------
! GFS initialize
!---------------
  subroutine GFS_initialize (Model, Statein, Stateout, Sfcprop,     &
                             Coupling, Grid, Tbd, Cldprop, Radtend, &
                             Diag, Interstitial, Init_parm)

    use module_microphysics, only: gsmconst
! Not available in FV3v0-CCPP demo
#if 0
    use cldwat2m_micro,      only: ini_micro
    use aer_cloud,           only: aer_cloud_init
#endif
    use module_ras,          only: ras_init
#ifdef OPENMP
    use omp_lib
#endif

    !--- interface variables
    type(GFS_control_type),      intent(inout) :: Model
    type(GFS_statein_type),      intent(inout) :: Statein(:)
    type(GFS_stateout_type),     intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),      intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),     intent(inout) :: Coupling(:)
    type(GFS_grid_type),         intent(inout) :: Grid(:)
    type(GFS_tbd_type),          intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),      intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),      intent(inout) :: Radtend(:)
    type(GFS_diag_type),         intent(inout) :: Diag(:)
    type(GFS_interstitial_type), intent(inout) :: Interstitial(:)
    type(GFS_init_type),         intent(in)    :: Init_parm

    !--- local variables
    integer :: nb
    integer :: nblks
    integer :: blkszmax
    integer :: nt
    integer :: nthreads
    integer :: ntrac
    real(kind=kind_phys), allocatable :: si(:)
    real(kind=kind_phys), parameter   :: p_ref = 101325.0d0


    nblks = size(Init_parm%blksz)
    ntrac = size(Init_parm%tracer_names)
    allocate (blksz(nblks))
    blksz(:) = Init_parm%blksz(:)

    !--- set control properties (including namelist read)
    call Model%init (Init_parm%nlunit, Init_parm%fn_nml,           &
                     Init_parm%me, Init_parm%master,               &
                     Init_parm%logunit, Init_parm%isc,             &
                     Init_parm%jsc, Init_parm%nx, Init_parm%ny,    &
                     Init_parm%levs, Init_parm%cnx, Init_parm%cny, &
                     Init_parm%gnx, Init_parm%gny,                 &
                     Init_parm%dt_dycore, Init_parm%dt_phys,       &
                     Init_parm%bdat, Init_parm%cdat,               &
                     Init_parm%tracer_names, Init_parm%blksz)

    call read_o3data  (Model%ntoz, Model%me, Model%master)
    call read_h2odata (Model%h2o_phys, Model%me, Model%master)

    blkszmax = 0
    do nb = 1,nblks
      call Statein      (nb)%create (Init_parm%blksz(nb), Model)
      call Stateout     (nb)%create (Init_parm%blksz(nb), Model)
      call Sfcprop      (nb)%create (Init_parm%blksz(nb), Model)
      call Coupling     (nb)%create (Init_parm%blksz(nb), Model)
      call Grid         (nb)%create (Init_parm%blksz(nb), Model)
      call Tbd          (nb)%create (Init_parm%blksz(nb), nb, Model)
      call Cldprop      (nb)%create (Init_parm%blksz(nb), Model)
      call Radtend      (nb)%create (Init_parm%blksz(nb), Model)
      !--- internal representation of diagnostics
      call Diag         (nb)%create (Init_parm%blksz(nb), Model)
      !--- maximum blocksize
      blkszmax = max(blkszmax, Init_parm%blksz(nb))
    enddo

#ifdef CCPP

#ifdef OPENMP
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

! Initialize the Interstitial data type in parallel so that
! each thread creates (touches) its Interstitial(nt) first
!$OMP parallel do default (shared) &
!$OMP            schedule (static,1) &
!$OMP            private  (nt)
    do nt=1,nthreads
      call Interstitial (nt)%create (blkszmax, Model)
    enddo
!$OMP end parallel do
#endif

    !--- populate the grid components
    call GFS_grid_populate (Grid, Init_parm%xlon, Init_parm%xlat, Init_parm%area)

    !--- read in and initialize ozone and water
    if (Model%ntoz > 0) then
      do nb = 1, nblks
        call setindxoz (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_o3, &
                        Grid(nb)%jindx2_o3, Grid(nb)%ddy_o3)
      enddo
    endif

    if (Model%h2o_phys) then
      do nb = 1, nblks
        call setindxh2o (Init_parm%blksz(nb), Grid(nb)%xlat_d, Grid(nb)%jindx1_h, &
                         Grid(nb)%jindx2_h, Grid(nb)%ddy_h)
      enddo
    endif

    !--- Call gfuncphys (funcphys.f) to compute all physics function tables.
    call gfuncphys ()

    call gsmconst (Model%dtp, Model%me, .TRUE.)

    !--- define sigma level for radiation initialization
    !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
    !--- ps is replaced with p0. The value of p0 uses that in http://www.emc.ncep.noaa.gov/officenotes/newernotes/on461.pdf
    !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
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

    !--- initialize Morrison-Gettleman microphysics
    if (Model%ncld == 2) then
      write(0,*) "Morrison-Gettleman microphysics not available in FV3v0-CCPP demo"
      stop
#if 0
      call ini_micro (Model%mg_dcs, Model%mg_qcvar, Model%mg_ts_auto_ice)
      call aer_cloud_init ()
#endif
    endif

    !--- initialize ras
    if (Model%ras) call ras_init (Model%levs, Model%me)

    !--- initialize soil vegetation
    call set_soilveg(Model%me, Model%isot, Model%ivegsrc, Model%nlunit)

    !--- lsidea initialization
    if (Model%lsidea) then
      print *,' LSIDEA is active but needs to be reworked for FV3 - shutting down'
      stop
      !--- NEED TO get the logic from the old phys/gloopb.f initialization area
    endif

    !--- sncovr may not exist in ICs from chgres.
    !--- FV3GFS handles this as part of the IC ingest
    !--- this not is placed here to alert users to the need to study
    !--- the FV3GFS_io.F90 module

  end subroutine GFS_initialize


!-------------------------------------------------------------------------
! time_vary_step
!-------------------------------------------------------------------------
!    routine called prior to radiation and physics steps to handle:
!      1) sets up various time/date variables
!      2) sets up various triggers
!      3) defines random seed indices for radiation (in a reproducible way)
!      5) interpolates coefficients for prognostic ozone calculation
!      6) performs surface data cycling via the GFS gcycle routine
!-------------------------------------------------------------------------
  subroutine GFS_time_vary_step (Model, Statein, Stateout, Sfcprop, Coupling, &
                                 Grid, Tbd, Cldprop, Radtend, Diag)

    use GFS_phys_time_vary_1,  only: GFS_phys_time_vary_1_run
    use GFS_phys_time_vary_2,  only: GFS_phys_time_vary_2_run
    use GFS_rad_time_vary,     only: GFS_rad_time_vary_run 
    implicit none

    !--- interface variables
    type(GFS_control_type),   intent(inout) :: Model
    type(GFS_statein_type),   intent(inout) :: Statein(:)
    type(GFS_stateout_type),  intent(inout) :: Stateout(:)
    type(GFS_sfcprop_type),   intent(inout) :: Sfcprop(:)
    type(GFS_coupling_type),  intent(inout) :: Coupling(:)
    type(GFS_grid_type),      intent(inout) :: Grid(:)
    type(GFS_tbd_type),       intent(inout) :: Tbd(:)
    type(GFS_cldprop_type),   intent(inout) :: Cldprop(:)
    type(GFS_radtend_type),   intent(inout) :: Radtend(:)
    type(GFS_diag_type),      intent(inout) :: Diag(:)

    ! CCPP error handling variables (not used)
    character(len=512) :: errmsg
    integer            :: errflg

    call GFS_phys_time_vary_1_run (Model, errmsg, errflg)

    call GFS_rad_time_vary_run (Model, Statein, Tbd, errmsg, errflg)

    call GFS_phys_time_vary_2_run (Grid, Model, Tbd, Sfcprop, Cldprop, Diag, errmsg, errflg)

  end subroutine GFS_time_vary_step


#ifndef CCPP
!-------------------------------------------------------------------------
! GFS stochastic_driver
!-------------------------------------------------------------------------
!    routine called prior to radiation and physics steps to handle:
!      1) sets up various time/date variables
!      2) sets up various triggers
!      3) defines random seed indices for radiation (in a reproducible way)
!      5) interpolates coefficients for prognostic ozone calculation
!      6) performs surface data cycling via the GFS gcycle routine
!-------------------------------------------------------------------------
  subroutine GFS_stochastic_driver (Model, Statein, Stateout, Sfcprop, Coupling, &
                                    Grid, Tbd, Cldprop, Radtend, Diag)

    use GFS_stochastics, only: GFS_stochastics_run
    !use memcheck, only: memcheck_run

    implicit none

    !--- interface variables
    type(GFS_control_type),   intent(in   ) :: Model
    type(GFS_statein_type),   intent(in   ) :: Statein
    type(GFS_stateout_type),  intent(in   ) :: Stateout
    type(GFS_sfcprop_type),   intent(in   ) :: Sfcprop
    type(GFS_coupling_type),  intent(inout) :: Coupling
    type(GFS_grid_type),      intent(in   ) :: Grid
    type(GFS_tbd_type),       intent(in   ) :: Tbd
    type(GFS_cldprop_type),   intent(in   ) :: Cldprop
    type(GFS_radtend_type),   intent(in   ) :: Radtend
    type(GFS_diag_type),      intent(inout) :: Diag

    ! CCPP error handling variables (not used)
    character(len=512) :: errmsg
    integer            :: errflg

    call GFS_stochastics_run(Model, Statein, Stateout, Sfcprop, Coupling, &
                             Grid, Tbd, Cldprop, Radtend, Diag, errmsg, errflg)

    !errmsg = 'end of GFS_stochastic_driver'
    !call memcheck_run(Model%sec, Tbd%blkno, errmsg, errflg)

  end subroutine GFS_stochastic_driver
#endif


!------------------
! GFS_grid_populate
!------------------
  subroutine GFS_grid_populate (Grid, xlon, xlat, area)
    use physcons,                 only: pi => con_pi

    implicit none

    type(GFS_grid_type)              :: Grid(:)
    real(kind=kind_phys), intent(in) :: xlon(:,:)
    real(kind=kind_phys), intent(in) :: xlat(:,:)
    real(kind=kind_phys), intent(in) :: area(:,:)

    !--- local variables
    integer :: nb, ix, blksz, i, j

    blksz = size(Grid(1)%xlon)

    nb = 1
    ix = 0
    do j = 1,size(xlon,2)
      do i = 1,size(xlon,1)
        ix=ix+1
        if (ix .gt. blksz) then
          nb = nb + 1
          ix = 1
        endif
        Grid(nb)%xlon(ix)   = xlon(i,j)
        Grid(nb)%xlat(ix)   = xlat(i,j)
        Grid(nb)%xlat_d(ix) = xlat(i,j) * 180.0_kind_phys/pi
        Grid(nb)%sinlat(ix) = sin(Grid(nb)%xlat(ix))
        Grid(nb)%coslat(ix) = sqrt(1.0_kind_phys - Grid(nb)%sinlat(ix)*Grid(nb)%sinlat(ix))
        Grid(nb)%area(ix)   = area(i,j)
        Grid(nb)%dx(ix)     = sqrt(area(i,j))
      enddo
    enddo

  end subroutine GFS_grid_populate


!-------------
! GFS finalize
!-------------
  subroutine GFS_finalize ()
  end subroutine GFS_finalize

end module GFS_driver
