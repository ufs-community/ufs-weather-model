!>\file rrtmg_sw_post
!! This file contains
      module rrtmg_sw_post
      contains

!>\defgroup rrtmg_sw_post GFS RRTMG scheme post
!! @{
!> \section arg_table_rrtmg_sw_post_init Argument Table
!!
      subroutine rrtmg_sw_post_init ()
      end subroutine rrtmg_sw_post_init
! PGI compiler does not accept lines longer than 264 characters, remove during pre-processing
#ifndef __PGI
!> \section arg_table_rrtmg_sw_post_run Argument Table
!! | local_name     | standard_name                                                                                  | long_name                                                                    | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                                                           | Fortran DDT containing FV3-GFS model control parameters                      | DDT      |    0 | GFS_control_type      |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                                                              | Fortran DDT containing FV3-GFS grid and interpolation related data           | DDT      |    0 | GFS_grid_type         |           | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                                                              | Fortran DDT containing FV3-GFS diagnotics data                               | DDT      |    0 | GFS_diag_type         |           | inout  | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                                                           | Fortran DDT containing FV3-GFS fields targetted for diagnostic output        | DDT      |    0 | GFS_radtend_type      |           | inout  | F        |
!! | Coupling       | FV3-GFS_Coupling_type                                                                          | Fortran DDT containing FV3-GFS fields to/from coupling with other components | DDT      |    0 | GFS_coupling_type     |           | inout  | F        |
!! | ltp            | extra_top_layer                                                                                | extra top layers                                                             | none     |    0 | integer               |           | in     | F        |
!! | nday           | daytime_points_dimension                                                                       | daytime points dimension                                                     | count    |    0 | integer               |           | in     | F        |
!! | lm             | vertical_layer_dimension_for_radiation                                                         | number of vertical layers for radiation calculation                          | count    |    0 | integer               |           | in     | F        |
!! | kd             | vertical_index_difference_between_inout_and_local                                              | vertical index difference between in/out and local                           | index    |    0 | integer               |           | in     | F        |
!! | htswc          | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | total sky heating rate due to shortwave radiation                            | K s-1    |    2 | real                  | kind_phys | in     | F        |
!! | htsw0          | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | clear sky heating rates due to shortwave radiation                           | K s-1    |    2 | real                  | kind_phys | in     | F        |
!! | sfcalb1        | surface_albedo_due_to_near_IR_direct                                                           | surface albedo due to near IR direct beam                                    | frac     |    1 | real                  | kind_phys | in     | F        |
!! | sfcalb2        | surface_albedo_due_to_near_IR_diffused                                                         | surface albedo due to near IR diffused beam                                  | frac     |    1 | real                  | kind_phys | in     | F        |
!! | sfcalb3        | surface_albedo_due_to_UV_and_VIS_direct                                                        | surface albedo due to UV+VIS direct beam                                     | frac     |    1 | real                  | kind_phys | in     | F        |
!! | sfcalb4        | surface_albedo_due_to_UV_and_VIS_diffused                                                      | surface albedo due to UV+VIS diffused beam                                   | frac     |    1 | real                  | kind_phys | in     | F        |
!! | scmpsw         | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes     | W m-2    |    1 | cmpfsw_type           |           | inout  | F        |
!! | errmsg         | error_message                                                                                  | error message for error handling in CCPP                                     | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                                                                     | error flag for error handling in CCPP                                        | flag     |    0 | integer               |           | out    | F        |
!!
#endif
      subroutine rrtmg_sw_post_run (Model, Grid, Diag, Radtend, Coupling, &
                 ltp, nday, lm, kd, htswc, htsw0,                         &
                 sfcalb1, sfcalb2, sfcalb3, sfcalb4, scmpsw, errmsg, errflg)

      use machine,                   only: kind_phys
      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type,   &
                                           cmpfsw_type
      use GFS_typedefs,              only: GFS_coupling_type,          &
                                           GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_diag_type

      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_grid_type),            intent(in)    :: Grid
      type(GFS_diag_type),            intent(inout) :: Diag
      integer, intent(in)                           ::  lm, kd, nday, ltp
      type(cmpfsw_type), dimension(size(Grid%xlon,1)), intent(inout) :: scmpsw
      real(kind=kind_phys), dimension(Size(Grid%xlon,1), Model%levr+LTP), intent(in) ::  htswc, htsw0
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: sfcalb1, sfcalb2, sfcalb3, sfcalb4
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer ::  k1, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lsswr) then
        if (nday > 0) then
          do k = 1, LM
            k1 = k + kd
            Radtend%htrsw(:,k) = htswc(:,k1)
          enddo
          ! --- repopulate the points above levr
          if (Model%levr < Model%levs) then
            do k = LM,Model%levs
              Radtend%htrsw (:,k) = Radtend%htrsw (:,LM)
            enddo
          endif
 
          if (Model%swhtr) then
            do k = 1, lm
               k1 = k + kd
               Radtend%swhc(:,k) = htsw0(:,k1)
             enddo
             ! --- repopulate the points above levr
             if (Model%levr < Model%levs) then
               do k = LM,Model%levs
                 Radtend%swhc(:,k) = Radtend%swhc(:,LM)
               enddo
             endif
          endif
 
!  --- surface down and up spectral component fluxes
!>  - Save two spectral bands' surface downward and upward fluxes for
!!    output.

          Coupling%nirbmdi(:) = scmpsw(:)%nirbm
          Coupling%nirdfdi(:) = scmpsw(:)%nirdf
          Coupling%visbmdi(:) = scmpsw(:)%visbm
          Coupling%visdfdi(:) = scmpsw(:)%visdf
 
          Coupling%nirbmui(:) = scmpsw(:)%nirbm * sfcalb1(:)
          Coupling%nirdfui(:) = scmpsw(:)%nirdf * sfcalb2(:)
          Coupling%visbmui(:) = scmpsw(:)%visbm * sfcalb3(:)
          Coupling%visdfui(:) = scmpsw(:)%visdf * sfcalb4(:)
 
        else                   ! if_nday_block
 
          Radtend%htrsw(:,:) = 0.0
 
          Radtend%sfcfsw = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          Diag%topfsw    = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw         = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
 
          Coupling%nirbmdi(:) = 0.0
          Coupling%nirdfdi(:) = 0.0
          Coupling%visbmdi(:) = 0.0
          Coupling%visdfdi(:) = 0.0
 
          Coupling%nirbmui(:) = 0.0
          Coupling%nirdfui(:) = 0.0
          Coupling%visbmui(:) = 0.0
          Coupling%visdfui(:) = 0.0
 
          if (Model%swhtr) then
            Radtend%swhc(:,:) = 0
          endif
 
        endif                  ! end_if_nday
 
! --- radiation fluxes for other physics processes
        Coupling%sfcnsw(:) = Radtend%sfcfsw(:)%dnfxc - Radtend%sfcfsw(:)%upfxc
        Coupling%sfcdsw(:) = Radtend%sfcfsw(:)%dnfxc
 
      endif                                ! end_if_lsswr

      end subroutine rrtmg_sw_post_run
 
!> \section arg_table_rrtmg_sw_post_finalize Argument Table
!!
      subroutine rrtmg_sw_post_finalize ()
      end subroutine rrtmg_sw_post_finalize
!! @}
      end module rrtmg_sw_post
