!>\file rrtmg_lw_pre.f90
!! This file contains a call to module_radiation_surface::setemis() to
!! setup surface emissivity for LW radiation.
      module rrtmg_lw_pre
      contains

!>\defgroup rrtmg_lw_pre GFS RRTMG scheme pre
!! @{
!> \section arg_table_rrtmg_lw_pre_init Argument Table
!!
      subroutine rrtmg_lw_pre_init ()
      end subroutine rrtmg_lw_pre_init 

!> \section arg_table_rrtmg_lw_pre_run Argument Table
!! | local_name     | standard_name                             | long_name                                                          | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|-------------------------------------------|--------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                      | Fortran DDT containing FV3-GFS model control parameters            | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                         | Fortran DDT containing FV3-GFS grid and interpolation related data | DDT      |    0 | GFS_grid_type         |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                      | Fortran DDT containing FV3-GFS surface fields                      | DDT      |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                      | Fortran DDT containing FV3-GFS radiation tendencies                | DDT      |    0 | GFS_radtend_type      |           | inout  | F        |
!! | im             | horizontal_loop_extent                    | horizontal loop extent                                             | count    |    0 | integer               |           | in     | F        |
!! | tsfg           | surface_ground_temperature_for_radiation  | surface ground temperature for radiation                           | K        |    1 | real                  | kind_phys | in     | F        |
!! | tsfa           | surface_air_temperature_for_radiation     | lowest model layer air temperature for radiation                   | K        |    1 | real                  | kind_phys | in     | F        |
!! | errmsg         | error_message                             | error message for error handling in CCPP                           | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                | error flag for error handling in CCPP                              | flag     |    0 | integer               |           | out    | F        |
!!
      subroutine rrtmg_lw_pre_run (Model, Grid, Sfcprop, Radtend, im, tsfg, tsfa, errmsg, errflg)
    
      use machine,                   only: kind_phys

      use GFS_typedefs,              only: GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_sfcprop_type         
      use module_radiation_surface,  only: setemis

      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
      type(GFS_grid_type),            intent(in)    :: Grid
      integer, intent(in)                           :: im
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) ::  tsfa, tsfg
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lslwr) then
!>  - Call module_radiation_surface::setemis(),to setup surface
!! emissivity for LW radiation.
        call setemis (Grid%xlon, Grid%xlat, Sfcprop%slmsk,        &        !  ---  inputs
                     Sfcprop%snowd, Sfcprop%sncovr, Sfcprop%zorl, &
                     tsfg, tsfa, Sfcprop%hprim, IM,               &
                      Radtend%semis)                              !  ---  outputs
      endif

       end subroutine rrtmg_lw_pre_run

!> \section arg_table_rrtmg_lw_pre_finalize Argument Table
!!
       subroutine rrtmg_lw_pre_finalize ()
       end subroutine rrtmg_lw_pre_finalize
!! @}
       end module rrtmg_lw_pre
