  module fv_arrays_stub_mod   

#include <fms_platform.h>
  use mpp_domains_mod,       only: domain2d
  use fms_io_mod,            only: restart_file_type
  use time_manager_mod,      only: time_type
  use mpp_mod,               only: mpp_broadcast
  use platform_mod,          only: r8_kind
  public

  integer, public, parameter :: R_GRID = r8_kind

  !Several 'auxiliary' structures are introduced here. These are for
  ! the internal use by certain modules, and although fv_atmos_type
  !  contains one of each of these structures all memory management
  !   is performed by the module in question.


!>@brief The type 'fv_grid_type' is made up of grid-dependent information from fv_grid_tools and fv_grid_utils.
!>@details It should not contain any user options (that goes in a different structure) nor data which
!! is altered outside of those two modules.
  type fv_grid_type
     real(kind=R_GRID), allocatable, dimension(:,:,:) :: grid_64, agrid_64
     real(kind=R_GRID), allocatable, dimension(:,:) :: dx_64, dy_64
     real(kind=R_GRID), allocatable, dimension(:,:) :: dxc_64, dyc_64
     real(kind=R_GRID), allocatable, dimension(:,:) :: dxa_64, dya_64

     real, allocatable, dimension(:,:,:) :: grid, agrid

     real, allocatable, dimension(:,:,:) :: e1,e2
     real, allocatable, dimension(:,:) :: dx, dy
     real, allocatable, dimension(:,:) :: dxc, dyc
     real, allocatable, dimension(:,:) :: dxa, dya
     real, allocatable, dimension(:,:) :: rdx, rdy
     real, allocatable, dimension(:,:) :: rdxc, rdyc
     real, allocatable, dimension(:,:) :: rdxa, rdya

     real(kind=R_GRID), allocatable :: ee1(:,:,:)
     real(kind=R_GRID), allocatable :: ee2(:,:,:)
     real(kind=R_GRID), allocatable :: ec1(:,:,:)
     real(kind=R_GRID), allocatable :: ec2(:,:,:)
     real(kind=R_GRID), allocatable :: ew(:,:,:,:)
     real(kind=R_GRID), allocatable :: es(:,:,:,:)


     !- 3D Super grid to contain all geometrical factors --
     ! the 3rd dimension is 9
     real, allocatable :: sin_sg(:,:,:)
     real, allocatable :: cos_sg(:,:,:)
     !--------------------------------------------------


     integer, dimension(:,:,:), allocatable :: iinta, jinta, iintb, jintb

     !Scalar data

     integer :: npx_g, npy_g, ntiles_g ! global domain

     logical :: g_sum_initialized = .false. !< Not currently used but can be useful
     logical:: sw_corner, se_corner, ne_corner, nw_corner

     real(kind=R_GRID) :: da_min, da_max, da_min_c, da_max_c

     real  :: acapN, acapS

     logical :: latlon = .false.
     logical :: cubed_sphere = .false.
     logical :: have_south_pole = .false.
     logical :: have_north_pole = .false.
     logical :: stretched_grid = .false.

     logical :: square_domain = .false.


     !! Convenience pointers

     integer, pointer :: grid_type !< Which type of grid to use. If 0, the equidistant gnomonic
                                   !< cubed-sphere will be used. If 4, a doubly-periodic
                                   !< f-plane cartesian grid will be used. If 5, a user-defined
                                   !< orthogonal grid will be used. If -1, the grid is read
                                   !< from INPUT/grid_spec.nc. Values 2, 3, 5, 6, and 7 are not
                                   !< supported and will likely not run. The default value is 0.

     logical, pointer :: nested   !< Whether this is a nested grid. .false. by default.
     logical, pointer :: regional !< Is this a (stand-alone) limited area regional domain?
     logical :: bounded_domain !< Is this a regional or nested domain?

  end type fv_grid_type

  type fv_flags_type

     !! FOR EACH VARIABLE IN FV_FLAGS:
     !! 1. Must be defined here:
     !! 2. Must be broadcast in fv_atmos_data
     !! 3. If a namelist entry, a pointer must
     !!    be defined and associated in fv_control
     !! 4. Must NOT appear in fv_current_grid_mod.
     !!    (this module will soon be removed)
     !! 5. Must be referenced through Atm%flagstruct,
     !!    not Atm%, unless a convenience
     !!    pointer is defined

!-----------------------------------------------------------------------
! Grid descriptor file setup
!-----------------------------------------------------------------------
   character(len=16) :: restart_resolution = 'both'
   character(len=80) :: grid_name = 'Gnomonic'
   character(len=120):: grid_file = 'Inline'
  integer      :: grid_type = 0     !< -1: read from file; 0: ED Gnomonic
!                                   !<  0: the "true" equal-distance Gnomonic grid
!                                   !<  1: the traditional equal-distance Gnomonic grid
!                                   !<  2: the equal-angular Gnomonic grid
!                                   !<  3: the lat-lon grid -- to be implemented
!                                   !<  4: double periodic boundary condition on Cartesian grid
!                                   !<  5: a user-defined orthogonal grid for stand alone regional model

!------------------------------------------
! Model Domain parameters
!------------------------------------------
   integer :: npx   !< Number of grid corners in the x-direction on one tile of the domain;
                    !< so one more than the number of grid cells across a tile. On the cubed sphere
                    !< this is one more than the number of cells across a cube face. Must be set.
   integer :: npy   !< Number of grid corners in the y-direction on one tile of the
                    !< domain. This value should be identical to npx on a cubed-sphere grid;
                    !< doubly periodic or nested grids do not have this restriction. Must be set.
   integer :: npz   !< Number of vertical levels. Each choice of npz comes with a
                    !< pre-defined set of hybrid sigma-pressure levels and model top
                    !< (see fv_eta.F90). Must be set.
   integer :: ntiles = 1  !< Number of tiles on the domain. For the cubed sphere, this
                          !< should be 6, one tile for each face of the cubed sphere; normally for
                          !< most other domains (including nested grids) this should be set to 1.
                          !< Must be set.
   integer :: ndims = 2   !< Lat-Lon Dims for Grid in Radians

  !>Convenience pointers
  integer, pointer :: grid_number

  end type fv_flags_type


  type fv_grid_bounds_type

     integer :: is,  ie,  js,  je
     integer :: isd, ied, jsd, jed
     integer :: isc, iec, jsc, jec

     integer :: ng = 3 !default

  end type fv_grid_bounds_type

  type fv_atmos_type

     logical :: allocated = .false.
     logical :: dummy = .false. ! same as grids_on_this_pe(n)
     integer :: grid_number = 1
     character(len=32) :: nml_filename = "input.nml"

     !Timestep-related variables.

     type(time_type) :: Time_init, Time, Run_length, Time_end, Time_step_atmos

     logical :: grid_active = .true. !Always active for now

!-----------------------------------------------------------------------
! Five prognostic state variables for the f-v dynamics
!-----------------------------------------------------------------------
    !! Convenience pointers
    integer, pointer :: npx, npy, npz, ng

     integer, allocatable, dimension(:) :: pelist

     type(fv_grid_bounds_type) :: bd

     type(fv_flags_type) :: flagstruct
     type(domain2D) :: domain

     type(domain2D) :: domain_for_coupler !< domain used in coupled model with halo = 1.

     !global tile and tile_of_mosaic only have a meaning for the CURRENT pe
     integer :: num_contact, npes_per_tile, global_tile, tile_of_mosaic, npes_this_grid
     integer :: layout(2), io_layout(2) = (/ 1,1 /)   !< layout: Processor layout on each tile.
                                                      !< The number of PEs assigned to a domain must equal
                                                      !< layout(1)*layout(2)*ntiles. Must be set.
                                                      !< io_layout: Layout of output files on each tile. 1,1 by default,
                                                      !< which combines all restart and history files on a tile into one file.
                                                      !< For 0,0, every process writes out its own restart and history files.
                                                      !< If not equal to 1,1, you will have to use mppnccombine to combine these
                                                      !< output files prior to post-processing, or if you want to change the
                                                      !< number of PEs. Both entries must divide the respective value in layout.

!!!!!!!!!!!!!!!!
! From fv_grid_tools
!!!!!!!!!!!!!!!!


     real    :: ptop

  type(fv_grid_type) :: gridstruct



!!!!!!!!!!!!!!
! From fv_io !
!!!!!!!!!!!!!!

     !Hold on to coarse-grid global grid, so we don't have to waste processor time getting it again when starting to do grid nesting
     real(kind=R_GRID), allocatable, dimension(:,:,:,:) :: grid_global

  integer :: atmos_axes(4)


  end type fv_atmos_type

contains

!>@brief The subroutine 'allocate_fv_atmos_type' allocates the fv_atmos_type
!>@details It includes an option to define dummy grids that have scalar and
!! small arrays defined as null 3D arrays.
  subroutine allocate_fv_atmos_type(Atm, isd_in, ied_in, jsd_in, jed_in, is_in, ie_in, js_in, je_in, &
       npx_in, npy_in, npz_in, ndims_in,  dummy, alloc_2d, ngrids_in)

    !WARNING: Before calling this routine, be sure to have set up the
    ! proper domain parameters from the namelists (as is done in
    ! fv_control.F90)

    implicit none
    type(fv_atmos_type), intent(INOUT), target :: Atm
    integer, intent(IN) :: isd_in, ied_in, jsd_in, jed_in, is_in, ie_in, js_in, je_in
    integer, intent(IN) :: npx_in, npy_in, npz_in, ndims_in
    logical, intent(IN) :: dummy, alloc_2d
    integer, intent(IN) :: ngrids_in
    integer:: isd, ied, jsd, jed, is, ie, js, je
    integer:: npx, npy, npz, ndims, ng

    !For 2D utility arrays
    integer:: isd_2d, ied_2d, jsd_2d, jed_2d, is_2d, ie_2d, js_2d, je_2d
    integer:: npx_2d, npy_2d, npz_2d, ndims_2d, ng_2d
    integer :: i,j,k, ns, n

    if (Atm%allocated) return

    if (dummy) then
       isd     =  0
       ied=   -1
       jsd=   0
       jed=   -1
       is=   0
       ie=   -1
       js=   0
       je=   -1
       npx=   1
       npy=   1
       npz=   1
       ndims=   1
    else
       isd     =  isd_in
       ied=   ied_in
       jsd=   jsd_in
       jed=   jed_in
       is=   is_in
       ie=   ie_in
       js=   js_in
       je=   je_in
       npx=   npx_in
       npy=   npy_in
       npz=   npz_in
       ndims=   ndims_in
    endif

    if ((.not. dummy) .or. alloc_2d) then
       isd_2d     =  isd_in
       ied_2d=   ied_in
       jsd_2d=   jsd_in
       jed_2d=   jed_in
       is_2d=   is_in
       ie_2d=   ie_in
       js_2d=   js_in
       je_2d=   je_in
       npx_2d=   npx_in
       npy_2d=   npy_in
       npz_2d=   npz_in
       ndims_2d=   ndims_in
    else
       isd_2d     =  0
       ied_2d=   -1
       jsd_2d=   0
       jed_2d=   -1
       is_2d=   0
       ie_2d=   -1
       js_2d=   0
       je_2d=   -1
       npx_2d=   1
       npy_2d=   1
       npz_2d=   npz_in !for ak, bk, which are 1D arrays and thus OK to allocate
       ndims_2d=   1
    endif


    !Convenience pointers
    Atm%npx => Atm%flagstruct%npx
    Atm%npy => Atm%flagstruct%npy
    Atm%npz => Atm%flagstruct%npz

    Atm%ng => Atm%bd%ng
    Atm%flagstruct%ndims = ndims_in


    allocate ( Atm%gridstruct% dx(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct% dx_64(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct%rdx(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct% dy(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct% dy_64(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%rdy(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )

    allocate ( Atm%gridstruct% dxc(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct% dxc_64(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%rdxc(isd_2d:ied_2d+1,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct% dyc(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct% dyc_64(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct%rdyc(isd_2d:ied_2d  ,jsd_2d:jed_2d+1) )

    allocate ( Atm%gridstruct% dxa(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct% dxa_64(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%rdxa(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct% dya(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct% dya_64(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )
    allocate ( Atm%gridstruct%rdya(isd_2d:ied_2d  ,jsd_2d:jed_2d  ) )

    allocate ( Atm%gridstruct%grid (isd_2d:ied_2d+1,jsd_2d:jed_2d+1,1:ndims_2d) )
    allocate ( Atm%gridstruct%grid_64 (isd_2d:ied_2d+1,jsd_2d:jed_2d+1,1:ndims_2d) )
    allocate ( Atm%gridstruct%agrid(isd_2d:ied_2d  ,jsd_2d:jed_2d  ,1:ndims_2d) )
    allocate ( Atm%gridstruct%agrid_64(isd_2d:ied_2d  ,jsd_2d:jed_2d  ,1:ndims_2d) )

    allocate ( Atm%gridstruct%  e1(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )
    allocate ( Atm%gridstruct%  e2(3,isd_2d:ied_2d+1,jsd_2d:jed_2d+1) )

    allocate (Atm%gridstruct%iinta(4, isd_2d:ied_2d ,jsd_2d:jed_2d), &
         Atm%gridstruct%jinta(4, isd_2d:ied_2d ,jsd_2d:jed_2d),  &
         Atm%gridstruct%iintb(4, is_2d:ie_2d+1 ,js_2d:je_2d+1), &
         Atm%gridstruct%jintb(4, is_2d:ie_2d+1 ,js_2d:je_2d+1) )

    !!Convenience pointers
    Atm%gridstruct%grid_type => Atm%flagstruct%grid_type
    Atm%flagstruct%grid_number => Atm%grid_number

    Atm%allocated = .true.
    if (dummy) Atm%dummy = .true.

  end subroutine allocate_fv_atmos_type

end module fv_arrays_stub_mod
