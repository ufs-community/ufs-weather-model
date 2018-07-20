!> \file GFS_zhao_carr_pre.f90
!! This file contains the subroutines that calculates physics/diagnotics variables 
!! before calling microphysics scheme:

      module GFS_zhao_carr_pre
      contains 

!> \defgroup GFS_zhao_carr_pre GFS Zhao-Carr Scheme pre
!! @{
!! \section arg_table_GFS_zhao_carr_pre_init Argument Table
!!
      subroutine GFS_zhao_carr_pre_init
      end subroutine GFS_zhao_carr_pre_init


!> \section arg_table_GFS_zhao_carr_pre_run Argument Table
!! | local_name     | standard_name                                              | long_name                                               | units       | rank |  type     |   kind    | intent | optional |
!! |----------------|------------------------------------------------------------|---------------------------------------------------------|-------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                     | horizontal loop extent                                  | count       |    0 | integer   |           | in     | F        |
!! | ix             | horizontal_dimension                                       | horizontal dimension                                    | count       |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                         | vertical layer dimension                                | count       |    0 | integer   |           | in     | F        |
!! | cwm            | cloud_condensed_water_specific_humidity_updated_by_physics | cloud condensed water specific humidity                 | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | clw1           | cloud_ice_specific_humidity                                | cloud ice specific humidity                             | kg kg-1     |    2 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                              | error message for error handling in CCPP                | none        |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                 | error flag for error handling in CCPP                   | flag        |    0 | integer   |           | out    | F        |
!!
      subroutine GFS_zhao_carr_pre_run (im, ix, levs, cwm, clw1, errmsg, errflg)
!
      use machine,               only: kind_phys

      implicit none

!     declare variables.

      integer,intent(in)   :: im, ix, levs
      integer              :: i,k
      real(kind=kind_phys),dimension(ix,levs), intent(in)  ::  cwm 
      real(kind=kind_phys),dimension(ix,levs), intent(out) ::  clw1

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i = 1, im
        do k = 1, levs
             clw1(i,k) = cwm(i,k) ! Stateout%gq0(:,:,Model%ntcw)
        enddo
      enddo

      end subroutine GFS_zhao_carr_pre_run

!> \section arg_table_GFS_zhao_carr_pre_finalize Argument Table
!!
      subroutine GFS_zhao_carr_pre_finalize
      end subroutine GFS_zhao_carr_pre_finalize
!! @}
      end module GFS_zhao_carr_pre

