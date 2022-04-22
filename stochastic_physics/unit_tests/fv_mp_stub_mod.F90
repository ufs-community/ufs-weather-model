      module fv_mp_stub_mod
! !USES:
      use fms_mod,         only : fms_init, fms_end
      use mpp_mod,         only : FATAL, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED, WARNING
      use mpp_mod,         only : mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_set_warn_level
      use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync
      use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
      use mpp_mod,         only : mpp_chksum, mpp_broadcast
      use mpp_mod,         only : mpp_send, mpp_recv, mpp_sync_self, EVENT_RECV, mpp_gather
      use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, FOLD_NORTH_EDGE, CGRID_NE
      use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
      use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D, mpp_get_ntile_count
      use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain
      use mpp_domains_mod, only : mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
      use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
      use mpp_domains_mod, only : mpp_check_field, mpp_define_layout 
      use mpp_domains_mod, only : mpp_get_neighbor_pe, mpp_define_mosaic, mpp_define_io_domain
      use mpp_domains_mod, only : NORTH, NORTH_EAST, EAST, SOUTH_EAST
      use mpp_domains_mod, only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST
      use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
      use mpp_domains_mod, only : mpp_group_update_initialized, mpp_do_group_update
      use mpp_domains_mod, only : mpp_create_group_update,mpp_reset_group_update_field
      use mpp_domains_mod, only : group_halo_update_type => mpp_group_update_type
      use mpp_domains_mod, only: nest_domain_type
      use mpp_parameter_mod, only : WUPDATE, EUPDATE, SUPDATE, NUPDATE, XUPDATE, YUPDATE
      use fv_arrays_stub_mod, only: fv_atmos_type, fv_grid_bounds_type
      use fms_io_mod, only: set_domain
      use mpp_mod, only : mpp_get_current_pelist, mpp_set_current_pelist
      use mpp_domains_mod, only : mpp_get_domain_shift
      use ensemble_manager_mod, only : get_ensemble_id

      implicit none
      private

      integer, parameter:: ng    = 3     ! Number of ghost zones required
      integer, parameter :: MAX_NNEST=20, MAX_NTILE=50

#include "mpif.h"
      integer, parameter :: XDir=1
      integer, parameter :: YDir=2
      integer :: commglobal, ierror, npes

      !need tile as a module variable so that some of the mp_ routines below will work
      integer::tile

      integer, allocatable, dimension(:)       :: npes_tile, tile1, tile2
      integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
      integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
      integer, allocatable, dimension(:,:)     :: layout2D, global_indices
      integer :: numthreads, gid, masterproc

      logical :: master

      integer :: this_pe_grid = 0
      integer, EXTERNAL :: omp_get_thread_num, omp_get_num_threads      

      integer :: npes_this_grid

      !! CLEANUP: these are currently here for convenience
      !! Right now calling switch_current_atm sets these to the value on the "current" grid
      !!  (as well as changing the "current" domain) 
      integer :: is, ie, js, je
      integer :: isd, ied, jsd, jed
      integer :: isc, iec, jsc, jec

      integer, allocatable :: grids_master_procs(:)
      integer, dimension(MAX_NNEST) :: tile_fine = 0 !Global index of LAST tile in a mosaic
      type(nest_domain_type) :: global_nest_domain !ONE structure for ALL levels of nesting
      public commglobal
      public mp_start, mp_assign_gid, mp_stop!, npes
      public domain_decomp
      public fill_corners, XDir, YDir
      public switch_current_domain, switch_current_Atm, broadcast_domains
      public setup_master
      public start_group_halo_update, complete_group_halo_update
      public group_halo_update_type, grids_master_procs, tile_fine
      public global_nest_domain, MAX_NNEST, MAX_NTILE, ng

      interface start_group_halo_update
        module procedure start_var_group_update_2d
        module procedure start_var_group_update_3d
        module procedure start_var_group_update_4d
        module procedure start_vector_group_update_2d
        module procedure start_vector_group_update_3d
      end interface start_group_halo_update

      INTERFACE fill_corners
        MODULE PROCEDURE fill_corners_2d_r4
        MODULE PROCEDURE fill_corners_2d_r8
        MODULE PROCEDURE fill_corners_xy_2d_r4
        MODULE PROCEDURE fill_corners_xy_2d_r8
        MODULE PROCEDURE fill_corners_xy_3d_r4
        MODULE PROCEDURE fill_corners_xy_3d_r8
      END INTERFACE

      INTERFACE fill_corners_agrid
        MODULE PROCEDURE fill_corners_agrid_r4
        MODULE PROCEDURE fill_corners_agrid_r8
      END INTERFACE

      INTERFACE fill_corners_cgrid
        MODULE PROCEDURE fill_corners_cgrid_r4
        MODULE PROCEDURE fill_corners_cgrid_r8
      END INTERFACE

      INTERFACE fill_corners_dgrid
        MODULE PROCEDURE fill_corners_dgrid_r4
        MODULE PROCEDURE fill_corners_dgrid_r8
      END INTERFACE

      !! The routines aggregate elements from many processes into one process. 
      ! WARNING only works with one level (ldim == 1)

      integer :: halo_update_type = 1

contains

        subroutine mp_assign_gid

          gid = mpp_pe()
          npes = mpp_npes()

        end subroutine mp_assign_gid

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!>@brief The subroutine 'mp_start' starts SPMD processes
        subroutine mp_start(commID, halo_update_type_in)
          integer, intent(in), optional :: commID
          integer, intent(in), optional :: halo_update_type_in

         integer :: ios
         integer :: unit

         masterproc = mpp_root_pe()
         commglobal = MPI_COMM_WORLD
         if( PRESENT(commID) )then
             commglobal = commID
         end if
         halo_update_type = halo_update_type_in

         numthreads = 1
!$OMP PARALLEL
!$OMP MASTER
!$       numthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL

         if ( mpp_pe()==mpp_root_pe() ) then
            master = .true.
         else
           master = .false.
         endif

         if (mpp_npes() > 1)  call MPI_BARRIER(commglobal, ierror)

      end subroutine mp_start
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

      subroutine setup_master(pelist_local)

        integer, intent(IN) :: pelist_local(:)

        if (ANY(gid == pelist_local)) then
        
           masterproc = pelist_local(1)
           master = (gid == masterproc)

        endif

      end subroutine setup_master


!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!>@brief The subroutine 'mp_stop' stops all SPMD processes
      subroutine mp_stop()

         call MPI_BARRIER(commglobal, ierror)
         if (gid==masterproc) print*, 'Stopping PEs : ', npes
         call fms_end()
        ! call MPI_FINALIZE (ierror)

      end subroutine mp_stop
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!>@brief The subroutine 'domain_decomp' sets up the domain decomposition.
      subroutine domain_decomp(npx,npy,nregions,grid_type,nested,layout,io_layout,bd,tile,square_domain,&
           npes_per_tile,domain,domain_for_coupler,num_contact,pelist)
         integer, intent(IN)  :: npx,npy,grid_type
         integer, intent(INOUT) :: nregions, tile
         logical, intent(IN):: nested
         integer, intent(INOUT) :: layout(2), io_layout(2)

         integer, allocatable :: pe_start(:), pe_end(:)

         integer :: nx,ny,n,num_alloc
         character(len=32) :: type = "unknown"
         logical :: is_symmetry 
         logical :: debug=.false.
         integer, allocatable :: tile_id(:)

         integer i
         integer :: npes_x, npes_y 

         integer, intent(INOUT) :: pelist(:)
         integer, intent(OUT) :: num_contact, npes_per_tile
         logical, intent(OUT) :: square_domain
         type(domain2D), intent(OUT) :: domain, domain_for_coupler
         type(fv_grid_bounds_type), intent(INOUT) :: bd

         nx = npx-1
         ny = npy-1

         npes_x = layout(1)
         npes_y = layout(2)


         call mpp_domains_init(MPP_DOMAIN_TIME)

         select case(nregions)
         case ( 1 )  ! Lat-Lon "cyclic"

            select case (grid_type)
            case (0,1,2) !Gnomonic nested grid
               if (nested) then
                  type = "Cubed-sphere nested grid"
               else
                  type = "Cubed-sphere, single face"
               end if
               nregions = 1
               num_contact = 0
               npes_per_tile = npes_x*npes_y !/nregions !Set up for concurrency
               is_symmetry = .true.
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )

               if ( npes_x == 0 ) then 
                  npes_x = layout(1)
               endif
               if ( npes_y == 0 ) then
                  npes_y = layout(2)
               endif

               if ( npes_x==npes_y .and. (npx-1)==((npx-1)/npes_x)*npes_x )  square_domain = .true.

               if ( (npx/npes_x < ng) .or. (npy/npes_y < ng) ) then
                  write(*,310) npes_x, npes_y, npx/npes_x, npy/npes_y
                 call mp_stop
                 call exit(1)
              endif
           
              layout = (/npes_x,npes_y/)
            case (3)   ! Lat-Lon "cyclic"
               type="Lat-Lon: cyclic"
               nregions = 4
               num_contact = 8
               if( mod(npes,nregions) .NE. 0 ) then
                  call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cyclic mosaic, npes should be multiple of nregions. ' // &
                                       'No test is done for Cyclic mosaic. ' )
                  return
               end if
               npes_per_tile = npes/nregions
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
               layout = (/1,npes_per_tile/) ! force decomp only in Lat-Direction
            case (4)   ! Cartesian, double periodic
               type="Cartesian: double periodic"
               nregions = 1
               num_contact = 2
               npes_per_tile = npes/nregions
               if(npes_x*npes_y == npes_per_tile) then
                  layout = (/npes_x,npes_y/)
               else
                  call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
               endif
            case (5)   ! latlon patch
               type="Lat-Lon: patch"
               nregions = 1
               num_contact = 0
               npes_per_tile = npes/nregions
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            case (6)   ! latlon strip
               type="Lat-Lon: strip"
               nregions = 1
               num_contact = 1
               npes_per_tile = npes/nregions
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            case (7)   ! Cartesian, channel
               type="Cartesian: channel"
               nregions = 1
               num_contact = 1
               npes_per_tile = npes/nregions
               call mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )
            end select

         case ( 6 )  ! Cubed-Sphere
            type="Cubic: cubed-sphere"
            if (nested) then
               call mpp_error(FATAL, 'For a nested grid with grid_type < 3 nregions_domain must equal 1.')
            endif
            nregions = 6
            num_contact = 12
            !--- cubic grid always have six tiles, so npes should be multiple of 6
            npes_per_tile = npes_x*npes_y
            call  mpp_define_layout( (/1,npx-1,1,npy-1/), npes_per_tile, layout )

            if ( npes_x == 0 ) then 
               npes_x = layout(1)
            endif
            if ( npes_y == 0 ) then
               npes_y = layout(2)
            endif

            if ( npes_x==npes_y .and. (npx-1)==((npx-1)/npes_x)*npes_x )  square_domain = .true.

            if ( (npx/npes_x < ng) .or. (npy/npes_y < ng) ) then
               write(*,310) npes_x, npes_y, npx/npes_x, npy/npes_y
 310           format('Invalid layout, NPES_X:',i4.4,'NPES_Y:',i4.4,'ncells_X:',i4.4,'ncells_Y:',i4.4)
               call mp_stop
               call exit(1)
            endif
           
            layout = (/npes_x,npes_y/)
         case default
            call mpp_error(FATAL, 'domain_decomp: no such test: '//type)
         end select

         allocate(layout2D(2,nregions), global_indices(4,nregions), npes_tile(nregions) )
         allocate(pe_start(nregions),pe_end(nregions))
         npes_tile = npes_per_tile
         do n = 1, nregions
            global_indices(:,n) = (/1,npx-1,1,npy-1/)
            layout2D(:,n)         = layout
               pe_start(n) = pelist(1) + (n-1)*layout(1)*layout(2)
            pe_end(n)   = pe_start(n) + layout(1)*layout(2) -1
         end do
         num_alloc=max(1,num_contact)
         allocate(tile1(num_alloc), tile2(num_alloc) )
         allocate(istart1(num_alloc), iend1(num_alloc), jstart1(num_alloc), jend1(num_alloc) )
         allocate(istart2(num_alloc), iend2(num_alloc), jstart2(num_alloc), jend2(num_alloc) )
 
         is_symmetry = .true.
         select case(nregions)
         case ( 1 )

            select case (grid_type)
            case (0,1,2) !Gnomonic nested grid
               !No contacts, don't need to do anything
            case (3)   ! Lat-Lon "cyclic"
               !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
               tile1(1) = 1; tile2(1) = 2
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
               !--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
               tile1(2) = 1; tile2(2) = 3
               istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
               istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
               !--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
               tile1(3) = 1; tile2(3) = 2
               istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
               istart2(3) = nx; iend2(3) = nx; jstart2(3) = 1;  jend2(3) = ny
               !--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH)
               tile1(4) = 1; tile2(4) = 3
               istart1(4) = 1;  iend1(4) = nx; jstart1(4) = ny;  jend1(4) = ny
               istart2(4) = 1;  iend2(4) = nx; jstart2(4) = 1;   jend2(4) = 1
               !--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
               tile1(5) = 2; tile2(5) = 4
               istart1(5) = 1;  iend1(5) = nx; jstart1(5) = 1;  jend1(5) = 1
               istart2(5) = 1;  iend2(5) = nx; jstart2(5) = ny; jend2(5) = ny
               !--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
               tile1(6) = 2; tile2(6) = 4
               istart1(6) = 1;  iend1(6) = nx; jstart1(6) = ny;  jend1(6) = ny
               istart2(6) = 1;  iend2(6) = nx; jstart2(6) = 1;   jend2(6) = 1
               !--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST)
               tile1(7) = 3; tile2(7) = 4
               istart1(7) = nx; iend1(7) = nx; jstart1(7) = 1;  jend1(7) = ny
               istart2(7) = 1;  iend2(7) = 1;  jstart2(7) = 1;  jend2(7) = ny
               !--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
               tile1(8) = 3; tile2(8) = 4
               istart1(8) = 1;  iend1(8) = 1;  jstart1(8) = 1;  jend1(8) = ny
               istart2(8) = nx; iend2(8) = nx; jstart2(8) = 1;  jend2(8) = ny
               is_symmetry = .false.
            case (4)   ! Cartesian, double periodic
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
               !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
               tile1(2) = 1; tile2(2) = 1
               istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
               istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
            case (5)   ! latlon patch

            case (6)   !latlon strip
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            case (7)   ! Cartesian, channel
               !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
               tile1(1) = 1; tile2(1) = 1
               istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
               istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            end select

         case ( 6 )  ! Cubed-Sphere
            !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
            tile1(1) = 1; tile2(1) = 2
            istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
            istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
            !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
            tile1(2) = 1; tile2(2) = 3
            istart1(2) = 1;  iend1(2) = nx; jstart1(2) = ny; jend1(2) = ny
            istart2(2) = 1;  iend2(2) = 1;  jstart2(2) = ny; jend2(2) = 1
            !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
            tile1(3) = 1; tile2(3) = 5
            istart1(3) = 1;  iend1(3) = 1;  jstart1(3) = 1;  jend1(3) = ny
            istart2(3) = nx; iend2(3) = 1;  jstart2(3) = ny; jend2(3) = ny
            !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
            tile1(4) = 1; tile2(4) = 6
            istart1(4) = 1;  iend1(4) = nx; jstart1(4) = 1;  jend1(4) = 1
            istart2(4) = 1;  iend2(4) = nx; jstart2(4) = ny; jend2(4) = ny
            !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
            tile1(5) = 2; tile2(5) = 3
            istart1(5) = 1;  iend1(5) = nx; jstart1(5) = ny; jend1(5) = ny
            istart2(5) = 1;  iend2(5) = nx; jstart2(5) = 1;  jend2(5) = 1
            !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
            tile1(6) = 2; tile2(6) = 4
            istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
            istart2(6) = nx; iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = 1
            !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
            tile1(7) = 2; tile2(7) = 6
            istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;  jend1(7) = 1
            istart2(7) = nx; iend2(7) = nx; jstart2(7) = ny; jend2(7) = 1
            !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
            tile1(8) = 3; tile2(8) = 4
            istart1(8) = nx; iend1(8) = nx; jstart1(8) = 1;  jend1(8) = ny
            istart2(8) = 1;  iend2(8) = 1;  jstart2(8) = 1;  jend2(8) = ny
            !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
            tile1(9) = 3; tile2(9) = 5
            istart1(9) = 1;  iend1(9) = nx; jstart1(9) = ny; jend1(9) = ny
            istart2(9) = 1;  iend2(9) = 1;  jstart2(9) = ny; jend2(9) = 1
            !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
            tile1(10) = 4; tile2(10) = 5
            istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny; jend1(10) = ny
            istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;  jend2(10) = 1
            !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
            tile1(11) = 4; tile2(11) = 6
            istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
            istart2(11) = nx; iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = 1
            !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
            tile1(12) = 5; tile2(12) = 6
            istart1(12) = nx; iend1(12) = nx; jstart1(12) = 1;  jend1(12) = ny
            istart2(12) = 1;  iend2(12) = 1;  jstart2(12) = 1;  jend2(12) = ny
         end select

         if ( ANY(pelist == gid) ) then
            allocate(tile_id(nregions))
            if( nested ) then
               if( nregions .NE. 1 ) then
                  call mpp_error(FATAL, 'domain_decomp: nregions should be 1 for nested region, contact developer')
               endif
               tile_id(1) = 7   ! TODO need update for multiple nests
            else
               do n = 1, nregions
                  tile_id(n) = n
               enddo
            endif
            call mpp_define_mosaic(global_indices, layout2D, domain, nregions, num_contact, tile1, tile2, &
                 istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                 pe_start=pe_start, pe_end=pe_end, symmetry=is_symmetry,              &
                 shalo = ng, nhalo = ng, whalo = ng, ehalo = ng, tile_id=tile_id, name = type)
            call mpp_define_mosaic(global_indices, layout2D, domain_for_coupler, nregions, num_contact, tile1, tile2, &
                 istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,                  &
                 pe_start=pe_start, pe_end=pe_end, symmetry=is_symmetry,                          &
                 shalo = 1, nhalo = 1, whalo = 1, ehalo = 1, tile_id=tile_id, name = type)
            deallocate(tile_id)
            call mpp_define_io_domain(domain, io_layout)
            call mpp_define_io_domain(domain_for_coupler, io_layout)

         endif

       deallocate(pe_start,pe_end)
       deallocate(layout2D, global_indices, npes_tile)
       deallocate(tile1, tile2)
       deallocate(istart1, iend1, jstart1, jend1)
       deallocate(istart2, iend2, jstart2, jend2)

       !--- find the tile number
       tile = (gid-pelist(1))/npes_per_tile+1 
       if (ANY(pelist == gid)) then
          npes_this_grid = npes_per_tile*nregions
          tile = tile
          call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
          call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
          
          bd%is = is
          bd%js = js
          bd%ie = ie
          bd%je = je

          bd%isd = isd
          bd%jsd = jsd
          bd%ied = ied
          bd%jed = jed

          bd%isc = is
          bd%jsc = js
          bd%iec = ie
          bd%jec = je

          if (debug .and. nregions==1) then
             tile=1
             write(*,200) tile, is, ie, js, je
             !   call mp_stop
             !   stop
          endif
200       format(i4.4, ' ', i4.4, ' ', i4.4, ' ', i4.4, ' ', i4.4, ' ')
       else
          
          bd%is = 0
          bd%js = 0
          bd%ie = -1
          bd%je = -1

          bd%isd = 0
          bd%jsd = 0
          bd%ied = -1
          bd%jed = -1

          bd%isc = 0
          bd%jsc = 0
          bd%iec = -1
          bd%jec = -1

       endif

      end subroutine domain_decomp
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

subroutine start_var_group_update_2d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group !< The data type that store information for group update
  real, dimension(:,:),         intent(inout) :: array !< The array which is having its halos points exchanged
  type(domain2D),               intent(inout) :: domain !< contains domain information
  integer,      optional,       intent(in)    :: flags !< Optional integer indicating which directions the data should be sent
  integer,      optional,       intent(in)    :: position  !< An optional argument indicating the position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete !< Optional argument indicating whether the halo updates
                                                          !! should be initiated immediately or wait for second pass_..._start call
  real                                        :: d_type
  logical                                     :: is_complete
! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  else
    call mpp_create_group_update(group, array, domain, flags=flags, position=position, &
             whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete .and. halo_update_type == 1) then 
     call mpp_start_group_update(group, domain, d_type)
  endif

end subroutine start_var_group_update_2d


subroutine start_var_group_update_3d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group !< The data type that store information for group update
  real, dimension(:,:,:),       intent(inout) :: array !< The array which is having its halos points exchanged
  type(domain2D),               intent(inout) :: domain !< contains domain information
  integer,           optional,  intent(in)    :: flags !< Optional integer indicating which directions the data should be sent 
  integer,           optional,  intent(in)    :: position !< An optional argument indicating the position
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete !< Optional argument indicating whether the halo updates
                                                          !! should be initiated immediately or wait for second pass_..._start call
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  else
    call mpp_create_group_update(group, array, domain, flags=flags, position=position, &
          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete .and. halo_update_type == 1 ) then
     call mpp_start_group_update(group, domain, d_type)
  endif

end subroutine start_var_group_update_3d

subroutine start_var_group_update_4d(group, array, domain, flags, position, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group !< The data type that store information for group update
  real, dimension(:,:,:,:),     intent(inout) :: array !< The array which is having its halos points exchanged
  type(domain2D),               intent(inout) :: domain !< contains domain information
  integer,           optional,  intent(in)    :: flags !< Optional integer indicating which directions the data should be sent 
  integer,           optional,  intent(in)    :: position !< An optional argument indicating the position
                                                          !! This is may be CORNER, but is CENTER by default
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete !< Optional argument indicating whether the halo updates
                                                          !! should be initiated immediately or wait for second pass_..._start call
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   array - The array which is having its halos points exchanged.
!  (in)      domain - contains domain information.
!  (in)      flags  - An optional integer indicating which directions the
!                       data should be sent.  
!  (in)      position - An optional argument indicating the position.  This is
!                       may be CORNER, but is CENTER by default.
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  integer :: dirflag

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,array)
  else
    call mpp_create_group_update(group, array, domain, flags=flags, position=position, &
              whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete .and. halo_update_type == 1 ) then
     call mpp_start_group_update(group, domain, d_type)
  endif

end subroutine start_var_group_update_4d



subroutine start_vector_group_update_2d(group, u_cmpt, v_cmpt, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group !< The data type that store information for group update
  real,       dimension(:,:),   intent(inout) :: u_cmpt, v_cmpt !< The nominal zonal (u) and meridional (v)
                                                                !! components of the vector pair that 
                                                                !! is having its halos points exchanged
  type(domain2d),               intent(inout) :: domain !< Contains domain decomposition information
  integer,            optional, intent(in)    :: flags !< Optional integer indicating which directions the data should be sent 
  integer,            optional, intent(in)    :: gridtype !< An optional flag, which may be one of A_GRID, BGRID_NE,
                                                          !! CGRID_NE or DGRID_NE, indicating where the two components of th 
                                                          !! vector are discretized
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete !< Optional argument indicating whether the halo updates
                                                          !! should be initiated immediately or wait for second pass_..._start call
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged. 
!  (in)      domain - Contains domain decomposition information.
!  (in)      flags - An optional integer indicating which directions the
!                        data should be sent. 
!  (in)      gridtype - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      CGRID_NE or DGRID_NE, indicating where the two components of the
!                      vector are discretized. 
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, domain, &
            flags=flags, gridtype=gridtype, &
            whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete .and. halo_update_type == 1 ) then
     call mpp_start_group_update(group, domain, d_type)
  endif

end subroutine start_vector_group_update_2d

subroutine start_vector_group_update_3d(group, u_cmpt, v_cmpt, domain, flags, gridtype, whalo, ehalo, shalo, nhalo, complete)
  type(group_halo_update_type), intent(inout) :: group !< The data type that store information for group update
  real,       dimension(:,:,:), intent(inout) :: u_cmpt, v_cmpt !! The nominal zonal (u) and meridional (v)
                                                                !! components of the vector pair that 
                                                                !! is having its halos points exchanged.
  type(domain2d),               intent(inout) :: domain !< Contains domain decomposition information
  integer,            optional, intent(in)    :: flags !< Optional integer indicating which directions the data should be sent
  integer,            optional, intent(in)    :: gridtype !< An optional flag, which may be one of A_GRID, BGRID_NE,
                                                          !! CGRID_NE or DGRID_NE, indicating where the two components of th 
                                                          !! vector are discretized
  integer,      optional,       intent(in)    :: whalo, ehalo, shalo, nhalo
  logical,      optional,       intent(in)    :: complete !< Optional argument indicating whether the halo updates
                                                          !! should be initiated immediately or wait for second pass_..._start call
  real                                        :: d_type
  logical                                     :: is_complete

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!                    This data will be used in do_group_pass.
!  (inout)   u_cmpt - The nominal zonal (u) component of the vector pair which
!                     is having its halos points exchanged.
!  (inout)   v_cmpt - The nominal meridional (v) component of the vector pair
!                     which is having its halos points exchanged. 
!  (in)      domain - Contains domain decomposition information.
!  (in)      flags - An optional integer indicating which directions the
!                        data should be sent. 
!  (in)      gridtype - An optional flag, which may be one of A_GRID, BGRID_NE,
!                      CGRID_NE or DGRID_NE, indicating where the two components of the
!                      vector are discretized. 
!  (in)      complete - An optional argument indicating whether the halo updates
!                       should be initiated immediately or wait for second 
!                       pass_..._start call.  Omitting complete is the same as 
!                       setting complete to .true.

  if (mpp_group_update_initialized(group)) then
    call mpp_reset_group_update_field(group,u_cmpt, v_cmpt)
  else
    call mpp_create_group_update(group, u_cmpt, v_cmpt, domain, &
            flags=flags, gridtype=gridtype, &
            whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)
  endif

  is_complete = .TRUE.
  if(present(complete)) is_complete = complete
  if(is_complete .and. halo_update_type == 1) then
     call mpp_start_group_update(group, domain, d_type)
  endif

end subroutine start_vector_group_update_3d


subroutine complete_group_halo_update(group, domain)
  type(group_halo_update_type), intent(inout) :: group !< The data type that store information for group update
  type(domain2d),               intent(inout) :: domain !< Contains domain decomposition information
  real                                        :: d_type

! Arguments: 
!  (inout)   group - The data type that store information for group update. 
!  (in)      domain - Contains domain decomposition information.

  if( halo_update_type == 1 ) then
  call mpp_complete_group_update(group, domain, d_type)
  else
    call mpp_do_group_update(group, domain, d_type)
  endif

end subroutine complete_group_halo_update



!Depreciated
subroutine broadcast_domains(Atm,current_pelist,current_npes)
  
  type(fv_atmos_type), intent(INOUT) :: Atm(:)
  integer, intent(IN) :: current_npes
  integer, intent(IN) :: current_pelist(current_npes)

  integer :: n, i
  integer :: ens_root_pe, ensemble_id

  !I think the idea is that each process needs to properly be part of a pelist,
  !the pelist on which the domain is currently defined is ONLY for the pes which have the domain.

  ! This is needed to set the proper pelist for the ensemble.  The pelist
  ! needs to include the non-nested+nested tile for the ensemble.
  ensemble_id = get_ensemble_id()
  ens_root_pe = (ensemble_id-1)*npes

  !Pelist needs to be set to ALL ensemble PEs for broadcast_domain to work
  call mpp_set_current_pelist((/ (i,i=ens_root_pe,npes-1+ens_root_pe) /))
  do n=1,size(Atm)
     call mpp_broadcast_domain(Atm(n)%domain)
     call mpp_broadcast_domain(Atm(n)%domain_for_coupler)
  end do
  call mpp_set_current_pelist(current_pelist)

end subroutine broadcast_domains

!depreciated
subroutine switch_current_domain(new_domain,new_domain_for_coupler)

  type(domain2D), intent(in), target :: new_domain, new_domain_for_coupler
  logical, parameter :: debug = .FALSE.

  !--- find the tile number
  !tile = mpp_pe()/npes_per_tile+1 
  !ntiles = mpp_get_ntile_count(new_domain)
  call mpp_get_compute_domain( new_domain, is,  ie,  js,  je  )
  isc = is ; jsc = js
  iec = ie ; jec = je
  call mpp_get_data_domain   ( new_domain, isd, ied, jsd, jed )
!  if ( npes_x==npes_y .and. (npx-1)==((npx-1)/npes_x)*npes_x )  square_domain = .true.

!  if (debug .AND. (gid==masterproc)) write(*,200) tile, is, ie, js, je
!200 format('New domain: ', i4.4, ' ', i4.4, ' ', i4.4, ' ', i4.4, ' ', i4.4, ' ')

  call set_domain(new_domain)


end subroutine switch_current_domain

!depreciated
subroutine switch_current_Atm(new_Atm, switch_domain)

  type(fv_atmos_type), intent(IN), target :: new_Atm
  logical, intent(IN), optional :: switch_domain
  logical, parameter :: debug = .false.
  logical :: swD


  call mpp_error(FATAL, "switch_current_Atm depreciated. call set_domain instead.")

!!$  if (debug .AND. (gid==masterproc)) print*, 'SWITCHING ATM STRUCTURES', new_Atm%grid_number
!!$  if (present(switch_domain)) then
!!$     swD = switch_domain
!!$  else
!!$     swD = .true.
!!$  end if
!!$  if (swD) call switch_current_domain(new_Atm%domain, new_Atm%domain_for_coupler)

!!$  if (debug .AND. (gid==masterproc)) WRITE(*,'(A, 6I5)') 'NEW GRID DIMENSIONS: ', &
!!$       isd, ied, jsd, jed, new_Atm%npx, new_Atm%npy

end subroutine switch_current_Atm

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     
      subroutine fill_corners_2d_r4(q, npx, npy, FILL, AGRID, BGRID)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: q
         integer, intent(IN):: npx,npy
         integer, intent(IN):: FILL  !< X-Dir or Y-Dir 
         logical, OPTIONAL, intent(IN) :: AGRID, BGRID 
         integer :: i,j

         if (present(BGRID)) then
            if (BGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-j  ,1-i  ) = q(i+1  ,1-j    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-j  ,npy+i) = q(i+1  ,npy+j  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+j,1-i  ) = q(npx-i,1-j    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+j,npy+i) = q(npx-i,npy+j  )  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              end select
            endif
          elseif (present(AGRID)) then
            if (AGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-i    ,1-j    ) = q(1-j    ,i        )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-i    ,npy-1+j) = q(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+i,1-j    ) = q(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+i,npy-1+j) = q(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng        
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                   enddo
                 enddo          
              end select
            endif
          endif

      end subroutine fill_corners_2d_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     
      subroutine fill_corners_2d_r8(q, npx, npy, FILL, AGRID, BGRID)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: q
         integer, intent(IN):: npx,npy
         integer, intent(IN):: FILL  ! <X-Dir or Y-Dir 
         logical, OPTIONAL, intent(IN) :: AGRID, BGRID 
         integer :: i,j

         if (present(BGRID)) then
            if (BGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-j  ,1-i  ) = q(i+1  ,1-j    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-j  ,npy+i) = q(i+1  ,npy+j  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+j,1-i  ) = q(npx-i,1-j    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+j,npy+i) = q(npx-i,npy+j  )  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng
                     if ((is==    1) .and. (js==    1)) q(1-i  ,1-j  ) = q(1-j  ,i+1    )  !SW Corner 
                     if ((is==    1) .and. (je==npy-1)) q(1-i  ,npy+j) = q(1-j  ,npy-i  )  !NW Corner
                     if ((ie==npx-1) .and. (js==    1)) q(npx+i,1-j  ) = q(npx+j,i+1    )  !SE Corner
                     if ((ie==npx-1) .and. (je==npy-1)) q(npx+i,npy+j) = q(npx+j,npy-i  )  !NE Corner
                    enddo
                 enddo
              end select
            endif
          elseif (present(AGRID)) then
            if (AGRID) then
              select case (FILL)
              case (XDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-i    ,1-j    ) = q(1-j    ,i        )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-i    ,npy-1+j) = q(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+i,1-j    ) = q(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+i,npy-1+j) = q(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
              case (YDir)
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo
              case default
                 do j=1,ng
                    do i=1,ng        
                       if ((is==    1) .and. (js==    1)) q(1-j    ,1-i    ) = q(i        ,1-j    )  !SW Corner 
                       if ((is==    1) .and. (je==npy-1)) q(1-j    ,npy-1+i) = q(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) q(npx-1+j,1-i    ) = q(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) q(npx-1+j,npy-1+i) = q(npx-1-i+1,npy-1+j)  !NE Corner
                   enddo
                 enddo          
              end select
            endif
          endif

      end subroutine fill_corners_2d_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     fill_corners_xy_2d_r8
      subroutine fill_corners_xy_2d_r8(x, y, npx, npy, DGRID, AGRID, CGRID, VECTOR)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: x !<(isd:ied  ,jsd:jed+1)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: y !<(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j

         real(kind=8) :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            call fill_corners_dgrid(x, y, npx, npy, mySign)
         elseif (present(CGRID)) then
            call fill_corners_cgrid(x, y, npx, npy, mySign)
         elseif (present(AGRID)) then
            call fill_corners_agrid(x, y, npx, npy, mySign)
         else
            call fill_corners_agrid(x, y, npx, npy, mySign)
         endif

      end subroutine fill_corners_xy_2d_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     fill_corners_xy_2d_r4
      subroutine fill_corners_xy_2d_r4(x, y, npx, npy, DGRID, AGRID, CGRID, VECTOR)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: x !<(isd:ied  ,jsd:jed+1)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: y !<(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j

         real(kind=4) :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            call fill_corners_dgrid(x, y, npx, npy, mySign)
         elseif (present(CGRID)) then
            call fill_corners_cgrid(x, y, npx, npy, mySign)
         elseif (present(AGRID)) then
            call fill_corners_agrid(x, y, npx, npy, mySign)
         else
            call fill_corners_agrid(x, y, npx, npy, mySign)
         endif

      end subroutine fill_corners_xy_2d_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     fill_corners_xy_3d_r8
      subroutine fill_corners_xy_3d_r8(x, y, npx, npy, npz, DGRID, AGRID, CGRID, VECTOR)
         real(kind=8), DIMENSION(isd:,jsd:,:), intent(INOUT):: x !<(isd:ied  ,jsd:jed+1)
         real(kind=8), DIMENSION(isd:,jsd:,:), intent(INOUT):: y !<(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy,npz
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j,k

         real(kind=8) :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            do k=1,npz
               call fill_corners_dgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(CGRID)) then
            do k=1,npz
               call fill_corners_cgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(AGRID)) then
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         else
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         endif

      end subroutine fill_corners_xy_3d_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !!
!     fill_corners_xy_3d_r4
      subroutine fill_corners_xy_3d_r4(x, y, npx, npy, npz, DGRID, AGRID, CGRID, VECTOR)
         real(kind=4), DIMENSION(isd:,jsd:,:), intent(INOUT):: x !<(isd:ied  ,jsd:jed+1)
         real(kind=4), DIMENSION(isd:,jsd:,:), intent(INOUT):: y !<(isd:ied+1,jsd:jed  )
         integer, intent(IN):: npx,npy,npz
         logical, OPTIONAL, intent(IN) :: DGRID, AGRID, CGRID, VECTOR
         integer :: i,j,k

         real(kind=4) :: mySign

         mySign = 1.0
         if (present(VECTOR)) then
            if (VECTOR) mySign = -1.0
         endif

         if (present(DGRID)) then
            do k=1,npz
               call fill_corners_dgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(CGRID)) then
            do k=1,npz
               call fill_corners_cgrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         elseif (present(AGRID)) then
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         else
            do k=1,npz
               call fill_corners_agrid(x(:,:,k), y(:,:,k), npx, npy, mySign)
            enddo
         endif

      end subroutine fill_corners_xy_3d_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_dgrid_r8(x, y, npx, npy, mySign)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=8), intent(IN) :: mySign 
         integer :: i,j

               do j=1,ng
                  do i=1,ng
                   !   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !SW Corner 
                   !   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mySign*y(j+1  ,npy-1+i)  !NW Corner
                   !   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mySign*y(npx-j,1-i    )  !SE Corner
                   !   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !NE Corner
                      if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) = mySign*y(1-j  ,i    )  !SW Corner 
                      if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) =        y(1-j  ,npy-i)  !NW Corner
                      if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) =        y(npx+j,i    )  !SE Corner
                      if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) = mySign*y(npx+j,npy-i)  !NE Corner
                  enddo
               enddo
               do j=1,ng
                  do i=1,ng
                   !  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !SW Corner 
                   !  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mySign*x(1-j    ,npy-i)  !NW Corner
                   !  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mySign*x(npx-1+j,i+1  )  !SE Corner
                   !  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !NE Corner
                     if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) = mySign*x(j      ,1-i  )  !SW Corner 
                     if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) =        x(j      ,npy+i)  !NW Corner
                     if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) =        x(npx-j  ,1-i  )  !SE Corner
                     if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) = mySign*x(npx-j  ,npy+i)  !NE Corner
                  enddo
               enddo

      end subroutine fill_corners_dgrid_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_dgrid_r4(x, y, npx, npy, mySign)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=4), intent(IN) :: mySign 
         integer :: i,j

               do j=1,ng
                  do i=1,ng
                   !   if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) =        y(j+1  ,1-i    )  !SW Corner 
                   !   if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) = mySign*y(j+1  ,npy-1+i)  !NW Corner
                   !   if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) = mySign*y(npx-j,1-i    )  !SE Corner
                   !   if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) =        y(npx-j,npy-1+i)  !NE Corner
                      if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j  ) = mySign*y(1-j  ,i    )  !SW Corner 
                      if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy+j) =        y(1-j  ,npy-i)  !NW Corner
                      if ((ie+1==npx) .and. (js  ==  1)) x(npx-1+i,1-j  ) =        y(npx+j,i    )  !SE Corner
                      if ((ie+1==npx) .and. (je+1==npy)) x(npx-1+i,npy+j) = mySign*y(npx+j,npy-i)  !NE Corner
                  enddo
               enddo
               do j=1,ng
                  do i=1,ng
                   !  if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) =        x(1-j    ,i+1  )  !SW Corner 
                   !  if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) = mySign*x(1-j    ,npy-i)  !NW Corner
                   !  if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) = mySign*x(npx-1+j,i+1  )  !SE Corner
                   !  if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) =        x(npx-1+j,npy-i)  !NE Corner
                     if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j    ) = mySign*x(j      ,1-i  )  !SW Corner 
                     if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy-1+j) =        x(j      ,npy+i)  !NW Corner
                     if ((ie+1==npx) .and. (js  ==  1)) y(npx+i  ,1-j    ) =        x(npx-j  ,1-i  )  !SE Corner
                     if ((ie+1==npx) .and. (je+1==npy)) y(npx+i  ,npy-1+j) = mySign*x(npx-j  ,npy+i)  !NE Corner
                  enddo
               enddo

      end subroutine fill_corners_dgrid_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_cgrid_r4(x, y, npx, npy, mySign)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=4), intent(IN) :: mySign
         integer :: i,j

                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j    ) =        y(j      ,1-i  )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy-1+j) = mySign*y(j      ,npy+i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) x(npx+i  ,1-j    ) = mySign*y(npx-j  ,1-i  )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) x(npx+i  ,npy-1+j) =        y(npx-j  ,npy+i)  !NE Corner
                     enddo
                  enddo
                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j  ) =        x(1-j  ,i    )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy+j) = mySign*x(1-j  ,npy-i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) y(npx-1+i,1-j  ) = mySign*x(npx+j,i    )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) y(npx-1+i,npy+j) =        x(npx+j,npy-i)  !NE Corner
                     enddo
                  enddo
      
      end subroutine fill_corners_cgrid_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_cgrid_r8(x, y, npx, npy, mySign)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=8), intent(IN) :: mySign
         integer :: i,j

                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) x(1-i    ,1-j    ) =        y(j      ,1-i  )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) x(1-i    ,npy-1+j) = mySign*y(j      ,npy+i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) x(npx+i  ,1-j    ) = mySign*y(npx-j  ,1-i  )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) x(npx+i  ,npy-1+j) =        y(npx-j  ,npy+i)  !NE Corner
                     enddo
                  enddo
                  do j=1,ng
                     do i=1,ng
                        if ((is  ==  1) .and. (js  ==  1)) y(1-i    ,1-j  ) =        x(1-j  ,i    )  !SW Corner 
                        if ((is  ==  1) .and. (je+1==npy)) y(1-i    ,npy+j) = mySign*x(1-j  ,npy-i)  !NW Corner
                        if ((ie+1==npx) .and. (js  ==  1)) y(npx-1+i,1-j  ) = mySign*x(npx+j,i    )  !SE Corner
                        if ((ie+1==npx) .and. (je+1==npy)) y(npx-1+i,npy+j) =        x(npx+j,npy-i)  !NE Corner
                     enddo
                  enddo
      
      end subroutine fill_corners_cgrid_r8
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_agrid_r4(x, y, npx, npy, mySign)
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=4), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=4), intent(IN) :: mySign
         integer :: i,j

                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) x(1-i    ,1-j    ) = mySign*y(1-j    ,i        )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) x(1-i    ,npy-1+j) =        y(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) x(npx-1+i,1-j    ) =        y(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) x(npx-1+i,npy-1+j) = mySign*y(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) y(1-j    ,1-i    ) = mySign*x(i        ,1-j    )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) y(1-j    ,npy-1+i) =        x(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) y(npx-1+j,1-i    ) =        x(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) y(npx-1+j,npy-1+i) = mySign*x(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo

      end subroutine fill_corners_agrid_r4
!
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ !
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv !
!
      subroutine fill_corners_agrid_r8(x, y, npx, npy, mySign)
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: x
         real(kind=8), DIMENSION(isd:,jsd:), intent(INOUT):: y
         integer, intent(IN):: npx,npy
         real(kind=8), intent(IN) :: mySign
         integer :: i,j

                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) x(1-i    ,1-j    ) = mySign*y(1-j    ,i        )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) x(1-i    ,npy-1+j) =        y(1-j    ,npy-1-i+1)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) x(npx-1+i,1-j    ) =        y(npx-1+j,i        )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) x(npx-1+i,npy-1+j) = mySign*y(npx-1+j,npy-1-i+1)  !NE Corner
                    enddo
                 enddo
                 do j=1,ng
                    do i=1,ng
                       if ((is==    1) .and. (js==    1)) y(1-j    ,1-i    ) = mySign*x(i        ,1-j    )  !SW Corner
                       if ((is==    1) .and. (je==npy-1)) y(1-j    ,npy-1+i) =        x(i        ,npy-1+j)  !NW Corner
                       if ((ie==npx-1) .and. (js==    1)) y(npx-1+j,1-i    ) =        x(npx-1-i+1,1-j    )  !SE Corner
                       if ((ie==npx-1) .and. (je==npy-1)) y(npx-1+j,npy-1+i) = mySign*x(npx-1-i+1,npy-1+j)  !NE Corner
                    enddo
                 enddo

      end subroutine fill_corners_agrid_r8
!



      end module fv_mp_stub_mod
!-------------------------------------------------------------------------------


