!>\file rrtmg_sw_pre.f90
!! This file contains a subroutine to module_radiation_surface::setalb() to
!! setup surface albedo for SW radiation.
      module rrtmg_sw_pre
      contains

!>\defgroup rrtmg_sw_pre GFS RRTMG scheme Pre
!! @{
!> \section arg_table_rrtmg_sw_pre_init Argument Table
!!
      subroutine rrtmg_sw_pre_init ()
      end subroutine rrtmg_sw_pre_init

!> \section arg_table_rrtmg_sw_pre_run Argument Table
!! | local_name     | standard_name                             | long_name                                                          | units    | rank |  type            |   kind    | intent | optional |
!! |----------------|-------------------------------------------|--------------------------------------------------------------------|----------|------|------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                      | Fortran DDT containing FV3-GFS model control parameters            | DDT      |    0 | GFS_control_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                         | Fortran DDT containing FV3-GFS grid and interpolation related data | DDT      |    0 | GFS_grid_type    |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                      | Fortran DDT containing FV3-GFS surface fields                      | DDT      |    0 | GFS_sfcprop_type |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                      | Fortran DDT containing FV3-GFS radiation tendencies                | DDT      |    0 | GFS_radtend_type |           | inout  | F        |
!! | im             | horizontal_loop_extent                    | horizontal loop extent                                             | count    |    0 | integer          |           | in     | F        |
!! | nday           | daytime_points_dimension                  | daytime points dimension                                           | count    |    0 | integer          |           | out    | F        |
!! | idxday         | daytime_points                            | daytime points                                                     | index    |    1 | integer          |           | out    | F        |
!! | tsfg           | surface_ground_temperature_for_radiation  | surface ground temperature for radiation                           | K        |    1 | real             | kind_phys | in     | F        |
!! | tsfa           | surface_air_temperature_for_radiation     | lowest model layer air temperature for radiation                   | K        |    1 | real             | kind_phys | in     | F        |
!! | sfcalb1        | surface_albedo_due_to_near_IR_direct      | surface albedo due to near IR direct beam                          | frac     |    1 | real             | kind_phys | out    | F        |
!! | sfcalb2        | surface_albedo_due_to_near_IR_diffused    | surface albedo due to near IR diffused beam                        | frac     |    1 | real             | kind_phys | out    | F        |
!! | sfcalb3        | surface_albedo_due_to_UV_and_VIS_direct   | surface albedo due to UV+VIS direct beam                           | frac     |    1 | real             | kind_phys | out    | F        |
!! | sfcalb4        | surface_albedo_due_to_UV_and_VIS_diffused | surface albedo due to UV+VIS diffused beam                         | frac     |    1 | real             | kind_phys | out    | F        |
!! | errmsg         | error_message                             | error message for error handling in CCPP                           | none     |    0 | character        | len=*     | out    | F        |
!! | errflg         | error_flag                                | error flag for error handling in CCPP                              | flag     |    0 | integer          |           | out    | F        |
!!
      subroutine rrtmg_sw_pre_run (Model, Grid, Sfcprop, Radtend, im, &
        nday, idxday, tsfg, tsfa, sfcalb1, sfcalb2, sfcalb3, sfcalb4, &
        errmsg, errflg)

      use machine,                   only: kind_phys

      use GFS_typedefs,              only: GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_sfcprop_type
      use module_radiation_surface,  only: NF_ALBD, setalb

      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
      type(GFS_grid_type),            intent(in)    :: Grid
      integer,                        intent(in)    :: im
      integer,                        intent(out)   :: nday
      integer, dimension(size(Grid%xlon,1)), intent(out) :: idxday
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) ::  tsfa, tsfg
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(out) :: sfcalb1, sfcalb2, sfcalb3, sfcalb4
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      ! Local variables
      integer :: i
      real(kind=kind_phys), dimension(size(Grid%xlon,1),NF_ALBD) :: sfcalb

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- ...  start radiation calculations
!           remember to set heating rate unit to k/sec!
!> -# Start SW radiation calculations
      if (Model%lsswr) then

!>  - Check for daytime points for SW radiation.
        nday = 0
        idxday = 0
        do i = 1, IM
          if (Radtend%coszen(i) >= 0.0001) then
            nday = nday + 1
            idxday(nday) = i
          endif
        enddo

!>  - Call module_radiation_surface::setalb() to setup surface albedo.
!!  for SW radiation.

        call setalb (Sfcprop%slmsk, Sfcprop%snowd, Sfcprop%sncovr,&       !  ---  inputs:
                     Sfcprop%snoalb, Sfcprop%zorl, Radtend%coszen,&
                     tsfg, tsfa, Sfcprop%hprim, Sfcprop%alvsf,    &
                     Sfcprop%alnsf, Sfcprop%alvwf, Sfcprop%alnwf, &
                     Sfcprop%facsf, Sfcprop%facwf, Sfcprop%fice,  &
                     Sfcprop%tisfc, IM,                           &
                     sfcalb)                                              !  ---  outputs

!> -# Approximate mean surface albedo from vis- and nir-  diffuse values.
        Radtend%sfalb(:) = max(0.01, 0.5 * (sfcalb(:,2) + sfcalb(:,4)))
      else
        nday = 0
        idxday = 0
        sfcalb = 0.0
      endif

      do i = 1, im
        sfcalb1(i) = sfcalb(i,1)
        sfcalb2(i) = sfcalb(i,2)
        sfcalb3(i) = sfcalb(i,3)
        sfcalb4(i) = sfcalb(i,4)
      enddo

      end subroutine rrtmg_sw_pre_run

!> \section arg_table_rrtmg_sw_pre_finalize Argument Table
!!
      subroutine rrtmg_sw_pre_finalize ()
      end subroutine rrtmg_sw_pre_finalize

!! @}
      end module rrtmg_sw_pre
