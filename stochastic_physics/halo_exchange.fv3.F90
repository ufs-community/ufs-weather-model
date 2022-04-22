module halo_exchange

!> This halo exchange routine is for host model fv3atm and requires FMS

use mpp_mod,                only: mpp_error, FATAL
use mpp_parameter_mod,      only: EUPDATE, WUPDATE, SUPDATE, NUPDATE
use mpp_domains_mod,        only: domain2d, mpp_update_domains
use mpp_domains_mod,        only: mpp_update_domains

implicit none
private

!--- utility routines
public ::  atmosphere_scalar_field_halo

contains

!>@brief The subroutine 'atmosphere_scalar_field_halo' is an API to return halo information 
!! of the current MPI_rank for an input scalar field.
!>@detail Up to three point haloes can be returned by this API which includes special handling for
!! the cubed-sphere tile corners. Output will be in (i,j,k) while input can be in (i,j,k) or 
!! horizontally-packed form (ix,k).
 
 subroutine atmosphere_scalar_field_halo (data, halo, isize, jsize, ksize, &
                                          isc, iec, jsc, jec, npx, npy, domain_for_coupler,data_p)
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
   real*8, dimension(1:isize,1:jsize,ksize), intent(inout) :: data !< output array to return the field with halo (i,j,k)
                                                                   !< optionally input for field already in (i,j,k) form
                                                                   !< sized to include the halo of the field (+ 2*halo)
   integer, intent(in) :: halo  !< size of the halo (must be less than 3)
   integer, intent(in) :: isize !< horizontal resolution in i-dir with haloes
   integer, intent(in) :: jsize !< horizontal resolution in j-dir with haloes
   integer, intent(in) :: ksize !< vertical resolution
   real*8, dimension(:,:), optional, intent(in) :: data_p !< optional input field in packed format (ix,k)
   integer, intent(in) :: isc, iec, jsc, jec, npx, npy
   type(domain2d), intent(inout) :: domain_for_coupler
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
     call mpp_update_domains(data, domain_for_coupler, flags=mpp_flags, complete=.true.)
   ! Not needed for cellular automata code
   !elseif (halo == 3) then
   !  call mpp_update_domains(data, Atm(mytile)%domain, flags=mpp_flags, complete=.true.)
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

end module halo_exchange
