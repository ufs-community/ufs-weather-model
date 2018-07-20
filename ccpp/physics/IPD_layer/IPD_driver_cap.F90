!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Semi-auto-generated cap module for the IPD_driver scheme
!!
!
module IPD_driver_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_api,                                        &
                      only: ccpp_t,                                    &
                            ccpp_field_get,                            &
                            ccpp_error
    use            :: IPD_typedefs,                                    &
                      only: IPD_init_type,                             &
                            IPD_control_type,                          &
                            IPD_data_type,                             &
                            IPD_restart_type,                          &
                            IPD_diag_type,                             &
                            IPD_interstitial_type
    use            :: IPD_driver,                                      &
                      only: IPD_initialize,                            &
                            IPD_setup_step,                            &
                            IPD_finalize
    use            :: machine,                                         &
                      only: kind_phys
    use            :: namelist_soilveg,                                &
                      only: salp_data, snupx, max_vegtyp
    implicit none

    private

    public :: ipd_initialize_cap,     &
              ipd_setup_step_cap,     &
              ipd_finalize_cap

    contains

    function ipd_initialize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        integer, allocatable                  :: dims(:)
        type(ccpp_t),                pointer  :: cdata
        type(IPD_control_type),      pointer  :: IPD_Control
        type(IPD_data_type),         pointer  :: IPD_Data(:)
        type(IPD_diag_type),         pointer  :: IPD_Diag(:)
        type(IPD_restart_type),      pointer  :: IPD_Restart
        type(IPD_interstitial_type), pointer  :: IPD_Interstitial(:)
        type(IPD_init_type),         pointer  :: Init_parm
        type(c_ptr)                           :: tmp
        real(kind=kind_phys),        pointer  :: l_snupx(:)
        real(kind=kind_phys),        pointer  :: l_salp_data

        ierr = 0

        call c_f_pointer(ptr, cdata)

        call ccpp_field_get(cdata, 'IPD_Control', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Control')
        end if
        call c_f_pointer(tmp, IPD_Control)

        call ccpp_field_get(cdata, 'IPD_Data', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Data')
        end if
        call c_f_pointer(tmp, IPD_Data, dims)
        deallocate(dims)

        call ccpp_field_get(cdata, 'IPD_Diag', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Diag')
        end if
        call c_f_pointer(tmp, IPD_Diag, dims)
        deallocate(dims)

        call ccpp_field_get(cdata, 'IPD_Restart', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Restart')
        end if
        call c_f_pointer(tmp, IPD_Restart)

        call ccpp_field_get(cdata, 'IPD_Interstitial', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Interstitial')
        end if
        call c_f_pointer(tmp, IPD_Interstitial, dims)
        deallocate(dims)

        call ccpp_field_get(cdata, 'Init_parm', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Init_parm')
        end if
        call c_f_pointer(tmp, Init_parm)

        call ccpp_field_get(cdata, 'salp_data', l_salp_data, ierr)
        call ccpp_field_get(cdata, 'snupx', l_snupx, ierr)

        call IPD_initialize(IPD_Control=IPD_Control,           &
                            IPD_Data=IPD_Data,                 &
                            IPD_Diag=IPD_Diag,                 &
                            IPD_Restart=IPD_Restart,           &
                            IPD_Interstitial=IPD_Interstitial, &
                            IPD_init_parm=Init_parm)

        l_snupx = snupx
        l_salp_data = salp_data

    end function ipd_initialize_cap

    function ipd_setup_step_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        integer, allocatable             :: dims(:)
        type(ccpp_t),           pointer  :: cdata
        type(IPD_control_type), pointer  :: IPD_Control
        type(IPD_data_type),    pointer  :: IPD_Data(:)
        type(IPD_diag_type),    pointer  :: IPD_Diag(:)
        type(IPD_restart_type), pointer  :: IPD_Restart
        type(c_ptr)                      :: tmp

        ierr = 0

        call c_f_pointer(ptr, cdata)

        call ccpp_field_get(cdata, 'IPD_Control', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Control')
        end if
        call c_f_pointer(tmp, IPD_Control)

        call ccpp_field_get(cdata, 'IPD_Data', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Data')
        end if
        call c_f_pointer(tmp, IPD_Data, dims)
        deallocate(dims)

        call ccpp_field_get(cdata, 'IPD_Diag', tmp, ierr, dims=dims)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Diag')
        end if
        call c_f_pointer(tmp, IPD_Diag, dims)
        deallocate(dims)

        call ccpp_field_get(cdata, 'IPD_Restart', tmp, ierr)
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve IPD_Restart')
        end if
        call c_f_pointer(tmp, IPD_Restart)

        call IPD_setup_step(IPD_Control=IPD_Control, &
                            IPD_Data=IPD_Data,       &
                            IPD_Diag=IPD_Diag,       &
                            IPD_Restart=IPD_Restart)

    end function IPD_setup_step_cap

    function ipd_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        ierr = 0

        call IPD_finalize()

    end function ipd_finalize_cap

end module IPD_driver_cap
