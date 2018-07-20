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
!! @brief Type definitions module.
!!
!! @details The types module provides definitions for
!!          atmospheic driver to call the CCPP.
!
module ccpp_types

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_ptr, c_funptr

    implicit none

    private
    public :: CCPP_STR_LEN,                                            &
              ccpp_t,                                                  &
              ccpp_field_t,                                            &
              ccpp_scheme_t,                                           &
              ccpp_suite_t,                                            &
              ccpp_group_t,                                            &
              ccpp_subcycle_t

    !> @var CCPP_STR_LEN Parameter defined for string lengths.
    integer, parameter                          :: CCPP_STR_LEN = 256

    !>
    !! @breif CCPP field type
    !!
    !! The field type contains all the information/meta-data and data
    !! for fields that need to be passed between the atmosphere driver
    !! and the physics drivers.
    type :: ccpp_field_t
            character(len=CCPP_STR_LEN)                       :: standard_name
            character(len=CCPP_STR_LEN)                       :: long_name
            character(len=CCPP_STR_LEN)                       :: units
            integer                                           :: rank
            integer, allocatable, dimension(:)                :: dims
            type(c_ptr)                                       :: ptr
    end type ccpp_field_t

    !>
    !! @breif CCPP scheme type
    !!
    !! The scheme type contains all the scheme information.
    !
    type :: ccpp_scheme_t
            character(:), allocatable                         :: name
            character(:), allocatable                         :: library
            character(:), allocatable                         :: version
            type(c_ptr)                                       :: scheme_hdl
            type(c_ptr)                                       :: library_hdl
    end type ccpp_scheme_t

    !>
    !! @breif CCPP subcycle type
    !!
    !! The subcycle type contains all the scheme names and the number of
    !! times the subcycle will loop. It is a direct mapping to the group
    !! suite subcycle XML.
    !
    type :: ccpp_subcycle_t
            integer                                           :: loop
            integer                                           :: schemes_max
            type(ccpp_scheme_t), allocatable, dimension(:)    :: schemes
    end type ccpp_subcycle_t

    !>
    !! @breif CCPP group type
    !!
    !! The group type contains all the subcycles and the name of
    !! the group call. It is a direct mapping to the group element in XML.
    !
    type :: ccpp_group_t
            character(:), allocatable                           :: name
            integer                                             :: subcycles_max
            type(ccpp_subcycle_t), allocatable, dimension(:)    :: subcycles
    end type ccpp_group_t

    !>
    !! @breif CCPP suite type
    !!
    !! The suite type contains all the group parts names and number of
    !! times the subcycle will loop. It is a direct mapping to the
    !! suite element in XML.
    !
    type :: ccpp_suite_t
            character(:), allocatable                           :: name
            character(:), allocatable                           :: library
            character(:), allocatable                           :: version
            type(ccpp_scheme_t)                                 :: init
            type(ccpp_scheme_t)                                 :: finalize
            integer                                             :: groups_max
            type(ccpp_group_t), allocatable, dimension(:)       :: groups
            logical                                             :: iscopy
    end type ccpp_suite_t

    !>
    !! @breif CCPP physics type.
    !!
    !! Generic type that contains all components to run the CCPP.
    !!
    !! - Array of fields to all the data needing to go
    !!   the physics drivers.
    !! - The suite definitions in a ccpp_suite_t type.
    !
    type :: ccpp_t
            type(c_ptr)                                         :: fields_idx
            type(ccpp_field_t), allocatable, dimension(:)       :: fields
            type(ccpp_suite_t)                                  :: suite
            logical                                             :: initialized = .false.
    end type ccpp_t

end module ccpp_types
