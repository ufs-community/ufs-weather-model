
!***********************************************************************
!*                   GNU Lesser General Public License                 
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it 
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or 
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be 
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty 
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.  
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'atmosphere' provides the interface for the
!! Cubed-Sphere FV dynamical core

module atmosphere_stub_mod

! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>block_control_mod</td>
!     <td>block_control_type</td>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>cp_air, rdgas, grav, rvgas, kappa, pstd_mks</td>
!   </tr>
!   <tr>
!     <td>field_manager_mod</td>
!     <td>MODEL_ATMOS</td>
!   </tr>
!   <tr>
!     <td>fms_mod</td>
!     <td>file_exist, open_namelist_file,close_file, error_mesg, FATAL,
!         check_nml_error, stdlog,write_version_number,set_domain,
!         mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_SUBCOMPONENT, 
!         clock_flag_default, nullify_domain</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_atmos_type, R_GRID</td>
!   </tr>
!   <tr>
!     <td>fv_control_mod</td>
!   </tr>
!   <tr>
!     <td>fv_dynamics_mod</td>
!     <td>fv_dynamics</td>
!   </tr>
!   <tr>
!     <td>fv_eta_mod</td>
!   </tr>
!   <tr>
!     <td>fv_fill_mod</td>
!     <td>fill_gfs</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>switch_current_Atm</td>
!   </tr>
!   <tr>
!     <td>fv_nesting_mod</td>
!     <td>twoway_nesting</td>
!   </tr>
!   <tr>
!   <tr>
!     <td>fv_restart_mod</td>
!     <td>fv_restart, fv_write_restart</td>
!   </tr>
!   <tr>
!     <td>fv_sg_mod</td>
!     <td>fv_subgrid_z</td>
!   </tr>
!   <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!   <tr>
!     <td>fv_update_phys_mod</td>
!     <td>fv_update_phys</td>
!   </tr>
!   <tr>
!     <td>mpp_mod</td>
!     <td>mpp_error, stdout, FATAL, NOTE, input_nml_file, mpp_root_pe,
!                    mpp_npes, mpp_pe, mpp_chksum,mpp_get_current_pelist,
!                    mpp_set_current_pelist</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod</td>
!     <td>mpp_get_data_domain, mpp_get_compute_domain, domain2d, mpp_update_domains</td>
!   </tr>
!   <tr>
!     <td>mpp_parameter_mod</td>
!     <td>EUPDATE, WUPDATE, SUPDATE, NUPDATE</td>
!   </tr>
!   <tr>
!     <td>tracer_manager_mod</td>
!     <td>get_tracer_index, get_number_tracers, NO_TRACER</td>
!   </tr>
!   <tr>
!     <td>xgrid_mod</td>
!     <td>grid_box_type</td>
!   </tr>
! </table>

#include <fms_platform.h>

!-----------------
! FMS modules:
!-----------------
use block_control_mod,      only: block_control_type
use constants_mod,          only: cp_air, rdgas, grav, rvgas, kappa, pstd_mks
use fms_mod,                only: file_exist, open_namelist_file,    &
                                  close_file, error_mesg, FATAL,     &
                                  check_nml_error, stdlog,           &
                                  write_version_number,              &
                                  set_domain,                        &
                                  mpp_clock_id, mpp_clock_begin,     &
                                  mpp_clock_end, CLOCK_SUBCOMPONENT, &
                                  clock_flag_default, nullify_domain
use mpp_mod,                only: mpp_error, stdout, FATAL, WARNING, NOTE, &
                                  input_nml_file, mpp_root_pe,       &
                                  mpp_npes, mpp_pe, mpp_chksum,      &
                                  mpp_get_current_pelist,            &
                                  mpp_set_current_pelist, mpp_sync
use mpp_parameter_mod,      only: EUPDATE, WUPDATE, SUPDATE, NUPDATE
use mpp_domains_mod,        only: domain2d, mpp_update_domains
use xgrid_mod,              only: grid_box_type
use kinddef

!-----------------
! FV core modules:
!-----------------
use fv_arrays_stub_mod,      only: fv_atmos_type, R_GRID, fv_grid_bounds_type
use fv_control_stub_mod,only: fv_control_init, ngrids
use mpp_domains_mod,    only:  mpp_get_data_domain, mpp_get_compute_domain

implicit none
private

!--- driver routines
public :: atmosphere_init_stub

!--- utility routines
public :: atmosphere_resolution,                                 &
          atmosphere_control_data, atmosphere_scalar_field_halo

!-----------------------------------------------------------------------
! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>
character(len=20)   :: mod_name = 'fvGFS/atmosphere_mod'

!---- private data ----
  public Atm, mygrid

  !These are convenience variables for local use only, and are set to values in Atm%
  real    :: dt_atmos
  real    :: zvir
  integer :: npx, npy, npz
  integer :: isc, iec, jsc, jec
  integer :: isd, ied, jsd, jed
  integer :: sec, seconds, days
  integer, dimension(:), allocatable :: id_tracerdt_dyn
  integer :: sphum, liq_wat, rainwat, ice_wat, snowwat, graupel, cld_amt  ! condensate species tracer indices

  integer :: mygrid = 1
  integer, allocatable :: pelist(:)
  logical, allocatable :: grids_on_this_pe(:)
  type(fv_atmos_type), allocatable, target :: Atm(:)

  integer :: id_udt_dyn, id_vdt_dyn

  real, parameter:: w0_big = 60.  ! to prevent negative w-tracer diffusion

!---dynamics tendencies for use in fv_subgrid_z and during fv_update_phys
  real, allocatable, dimension(:,:,:)   :: u_dt, v_dt, t_dt
  real, allocatable                     :: pref(:,:), dum1d(:)


contains


!>@brief The subroutine 'atmosphere_init' is an API to initialize the FV3 dynamical core,
!! including the grid structures, memory, initial state (self-initialization or restart), 
 subroutine atmosphere_init_stub (Grid_box)

#ifdef OPENMP
   use omp_lib
#endif

   type(grid_box_type), intent(inout) :: Grid_box
!--- local variables ---
   integer :: i, n
!  integer :: itrac
   logical :: do_atmos_nudge
   character(len=32) :: tracer_name, tracer_units
   real    :: ps1, ps2
   integer :: nthreads, ierr
   integer :: nlunit = 9999
   character (len = 64) :: fn_nml = 'input.nml'

   allocate(pelist(mpp_npes()))
   call mpp_get_current_pelist(pelist)


   call fv_control_init( Atm, dt_atmos, mygrid, grids_on_this_pe)  ! allocates Atm components; sets mygrid


!----- write version and namelist to log file -----
   call write_version_number ( 'fvGFS/ATMOSPHERE_MOD', version )

!-----------------------------------

   npx   = Atm(mygrid)%npx
   npy   = Atm(mygrid)%npy
   npz   = Atm(mygrid)%npz

   isc = Atm(mygrid)%bd%isc
   iec = Atm(mygrid)%bd%iec
   jsc = Atm(mygrid)%bd%jsc
   jec = Atm(mygrid)%bd%jec

   isd = isc - Atm(mygrid)%bd%ng
   ied = iec + Atm(mygrid)%bd%ng
   jsd = jsc - Atm(mygrid)%bd%ng
   jed = jec + Atm(mygrid)%bd%ng

   ! Allocate grid variables to be used to calculate gradient in 2nd order flux exchange
   ! This data is only needed for the COARSEST grid.
   !call switch_current_Atm(Atm(mygrid))
   call set_domain(Atm(mygrid)%domain)

   allocate(Grid_box%dx    (   isc:iec  , jsc:jec+1))
   allocate(Grid_box%dy    (   isc:iec+1, jsc:jec  ))
   Grid_box%dx    (   isc:iec  , jsc:jec+1) = Atm(mygrid)%gridstruct%dx    (   isc:iec,   jsc:jec+1)
   Grid_box%dy    (   isc:iec+1, jsc:jec  ) = Atm(mygrid)%gridstruct%dy    (   isc:iec+1, jsc:jec  )

#ifdef OPENMP
   nthreads = omp_get_max_threads()
#else
   nthreads = 1
#endif
!  --- initiate the start for a restarted regional forecast


#ifdef DEBUG
   call nullify_domain()
#endif

   call set_domain(Atm(mygrid)%domain)
      
 end subroutine atmosphere_init_stub


!>@brief The subroutine 'atmospehre_resolution' is an API to return the local 
!! extents of the current MPI-rank or the global extents of the current 
!! cubed-sphere tile.
 subroutine atmosphere_resolution (i_size, j_size, global)
   integer, intent(out)          :: i_size, j_size
   logical, intent(in), optional :: global
   logical :: local

   local = .true.
   if( PRESENT(global) ) local = .NOT.global

   if( local ) then
       i_size = iec - isc + 1
       j_size = jec - jsc + 1
   else
       i_size = npx - 1
       j_size = npy - 1
   end if

 end subroutine atmosphere_resolution

 subroutine atmosphere_control_data (i1, i2, j1, j2, kt, p_hydro, hydro, tile_num)
   integer, intent(out)           :: i1, i2, j1, j2, kt
   logical, intent(out), optional :: p_hydro, hydro
   integer, intent(out), optional :: tile_num
   i1 = Atm(mygrid)%bd%isc
   i2 = Atm(mygrid)%bd%iec
   j1 = Atm(mygrid)%bd%jsc
   j2 = Atm(mygrid)%bd%jec
   kt = Atm(mygrid)%npz

   if (present(tile_num)) tile_num = Atm(mygrid)%tile_of_mosaic

 end subroutine atmosphere_control_data



 subroutine set_atmosphere_pelist ()
   call mpp_set_current_pelist(Atm(mygrid)%pelist, no_sync=.TRUE.)
 end subroutine set_atmosphere_pelist


!>@brief The subroutine 'atmosphere_scalar_field_halo' is an API to return halo information 
!! of the current MPI_rank for an input scalar field.
!>@detail Up to three point haloes can be returned by this API which includes special handling for
!! the cubed-sphere tile corners. Output will be in (i,j,k) while input can be in (i,j,k) or 
!! horizontally-packed form (ix,k).
 subroutine atmosphere_scalar_field_halo (data, halo, isize, jsize, ksize, data_p)
   !--------------------------------------------------------------------
   ! data   - output array to return the field with halo (i,j,k)
   !          optionally input for field already in (i,j,k) form
   !          sized to include the halo of the field (+ 2*halo)
   ! halo   - size of the halo (must be less than 3)
   ! ied    - horizontal resolution in i-dir with haloes
   ! jed    - horizontal resolution in j-dir with haloes
   ! ksize  - vertical resolution
   ! data_p - optional input field in packed format (ix,k)  
   !--------------------------------------------------------------------
   !--- interface variables ---
   real(kind=kind_phys), dimension(1:isize,1:jsize,ksize), intent(inout) :: data !< output array to return the field with halo (i,j,k)
                                                                                 !< optionally input for field already in (i,j,k) form
                                                                                 !< sized to include the halo of the field (+ 2*halo)
   integer, intent(in) :: halo  !< size of the halo (must be less than 3)
   integer, intent(in) :: isize !< horizontal resolution in i-dir with haloes
   integer, intent(in) :: jsize !< horizontal resolution in j-dir with haloes
   integer, intent(in) :: ksize !< vertical resolution
   real(kind=kind_phys), dimension(:,:), optional, intent(in) :: data_p !< optional input field in packed format (ix,k)
   !--- local variables ---
   integer :: i, j, k
   integer :: ic, jc
   character(len=44) :: modname = 'atmosphere_mod::atmosphere_scalar_field_halo'
   integer :: mpp_flags

   !--- perform error checking
   if (halo .gt. 3) call mpp_error(FATAL, modname//' - halo.gt.3 requires extending the MPP domain to support')
   ic = isize - 2 * halo
   jc = jsize - 2 * halo

   !--- if packed data is present, unpack it into the two-dimensional data array
   if (present(data_p)) then
     if (ic*jc .ne. size(data_p,1)) call mpp_error(FATAL, modname//' - incorrect sizes for incoming &
                                                  &variables data and data_p')
     data = 0.
!$OMP parallel do default (none) &
!$OMP              shared (data, data_p, halo, ic, jc, ksize) &
!$OMP             private (i, j, k)
     do k = 1, ksize
       do j = 1, jc
         do i = 1, ic
           data(i+halo, j+halo, k) = data_p(i + (j-1)*ic, k)
         enddo
       enddo
     enddo
   endif

   mpp_flags = EUPDATE + WUPDATE + SUPDATE + NUPDATE
   if (halo == 1) then
     call mpp_update_domains(data, Atm(mygrid)%domain_for_coupler, flags=mpp_flags, complete=.true.)
   elseif (halo == 3) then
     call mpp_update_domains(data, Atm(mygrid)%domain, flags=mpp_flags, complete=.true.)
   else
     call mpp_error(FATAL, modname//' - unsupported halo size')
   endif

   !--- fill the halo points when at a corner of the cubed-sphere tile 
   !--- interior domain corners are handled correctly
   if ( (isc==1) .or. (jsc==1) .or. (iec==npx-1) .or. (jec==npy-1) ) then
     do k = 1, ksize
       do j=1,halo
         do i=1,halo
           if ((isc==    1) .and. (jsc==    1)) data(halo+1-j ,halo+1-i ,k) = data(halo+i     ,halo+1-j ,k)  !SW Corner
           if ((isc==    1) .and. (jec==npy-1)) data(halo+1-j ,halo+jc+i,k) = data(halo+i     ,halo+jc+j,k)  !NW Corner
           if ((iec==npx-1) .and. (jsc==    1)) data(halo+ic+j,halo+1-i ,k) = data(halo+ic-i+1,halo+1-j ,k)  !SE Corner
           if ((iec==npx-1) .and. (jec==npy-1)) data(halo+ic+j,halo+jc+i,k) = data(halo+ic-i+1,halo+jc+j,k)  !NE Corner
         enddo
       enddo
     enddo
   endif

   return
 end subroutine atmosphere_scalar_field_halo


end module atmosphere_stub_mod
