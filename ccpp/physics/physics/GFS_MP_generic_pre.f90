!> \file GFS_MP_generic_pre.f90
!! This file contains the subroutines that calculates physics/diagnotics variables 
!! before calling microphysics scheme:

      module GFS_MP_generic_pre
      contains 

!> \defgroup GFS_MP_generic_pre GFS MP generic pre
!! @{
!! \section arg_table_GFS_MP_generic_pre_init Argument Table
!!
      subroutine GFS_MP_generic_pre_init
      end subroutine GFS_MP_generic_pre_init


!> \section arg_table_GFS_MP_generic_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                                | units       | rank |  type     |   kind    | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------------------------------------|-------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                 | horizontal loop extent                                                   | count       |    0 | integer   |           | in     | F        |
!! | ix             | horizontal_dimension                                   | horizontal dimension                                                     | count       |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                     | vertical layer dimension                                                 | count       |    0 | integer   |           | in     | F        |
!! | clw1           | cloud_ice_specific_humidity                            | cloud ice specific humidity                                              | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | clw2           | cloud_liquid_water_specific_humidity                   | cloud water specific humidity                                            | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | ldiag3d        | flag_diagnostics_3D                                    | logical flag for 3D diagnostics                                          | flag        |    0 | logical   |           | in     | F        |
!! | ntcw           | index_for_liquid_cloud_condensate                      | cloud condensate index in tracer array(3)                                | index       |    0 | integer   |           | in     | F        |
!! | ncld           | number_of_hydrometeors                                 | number of hydrometeors(1 for Z-C)                                        | count       |    0 | integer   |           | in     | F        |
!! | num_p3d        | array_dimension_of_microphysics                        | number of 3D arrays needed for microphysics                              | count       |    0 | integer   |           | in     | F        |
!! | t              | air_temperature_updated_by_physics                     | layer mean air temperature                                               | K           |    2 | real      | kind_phys | in     | F        |
!! | q              | water_vapor_specific_humidity_updated_by_physics       | water vapor specific humidity                                            | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | save_t         | air_temperature_save                                   | air temperature before entering a physics scheme                         | K           |    2 | real      | kind_phys | out    | F        |
!! | save_qv        | water_vapor_specific_humidity_save                     | water vapor specific humidity before entering a physics scheme           | kg kg-1     |    2 | real      | kind_phys | out    | F        |
!! | save_qcw       | cloud_condensed_water_specific_humidity_save           | cloud condensed water specific humidity before entering a physics scheme | kg kg-1     |    2 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                                 | none        |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                                    | flag        |    0 | integer   |           | out    | F        |
!!
      subroutine GFS_MP_generic_pre_run(im, ix, levs, clw1, clw2, &
                  ldiag3d, ntcw, ncld, num_p3d, t, q,             &
                  save_t, save_qv, save_qcw, errmsg, errflg)
!
      use machine,               only: kind_phys
      use physcons,              only:  con_g

      implicit none
!    
!     declare variables.
!
      integer, intent(in) :: im, ix, levs, ntcw, ncld, num_p3d
      integer             :: n,i,k
      logical, intent(in) :: ldiag3d
      real(kind=kind_phys),dimension(ix,levs), intent(in) :: t,q,  &
                                                     clw1,clw2
      real(kind=kind_phys),dimension(ix,levs), intent(out) ::   &
                                          save_t, save_qv, save_qcw
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (ldiag3d) then
        do i = 1, ix
          do k = 1, levs
            !CCPP dtdt(i,k)   = t(i,k)
            !CCPP dqdt1(i,k) =  q(i,k)
            save_t(i,k) = t(i,k)
            save_qv(i,k) = q(i,k)
          enddo
        end do
        !in FV3GFS v0 OP: ntcw=3, ncld=1, num_p3d=4, ntrac=3
        do n=ntcw,ntcw+ncld-1
          if (n == ntcw .and. num_p3d == 4 ) then
              do i = 1, ix
                do k = 1, levs
                   !CCPP dqdt3(i,k) = clw1(i,k)+clw2(i,k)   !
                   save_qcw(i,k) = clw1(i,k)+clw2(i,k)
                enddo
               enddo
          endif
        enddo
      else
          save_t = 0.0
          save_qv = 0.0
          save_qcw = 0.0
      endif

      end subroutine GFS_MP_generic_pre_run

!> \section arg_table_GFS_MP_generic_pre_finalize Argument Table
!!
      subroutine GFS_MP_generic_pre_finalize
      end subroutine GFS_MP_generic_pre_finalize
!! @}
      end module GFS_MP_generic_pre

