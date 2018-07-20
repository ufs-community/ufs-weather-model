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
!! @brief CCPP API module.
!!
!! @details This module exposes calls to CCPP routines to the users.
!
module ccpp_api

    use ccpp_types,         only: CCPP_STR_LEN,                        &
                                  ccpp_t
    use ccpp_errors,        only: ccpp_error,                          &
                                  ccpp_debug
    use ccpp,               only: ccpp_init,                           &
                                  ccpp_finalize,                       &
                                  ccpp_initialized
    use ccpp_fcall,         only: ccpp_physics_init,                   &
                                  ccpp_physics_run,                    &
                                  ccpp_physics_finalize
    use ccpp_fields,        only: ccpp_field_add,                      &
                                  ccpp_field_get
    use ccpp_memory,        only: ccpp_memory_usage

    implicit none

    public :: CCPP_STR_LEN,                                            &
              ccpp_t,                                                  &
              ccpp_error,                                              &
              ccpp_debug,                                              &
              ccpp_init,                                               &
              ccpp_finalize,                                           &
              ccpp_physics_init,                                       &
              ccpp_physics_run,                                        &
              ccpp_physics_finalize,                                   &
              ccpp_field_add,                                          &
              ccpp_initialized,                                        &
              ccpp_memory_usage

end module ccpp_api
