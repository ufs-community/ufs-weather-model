!> \file GFS_suite_interstitial.f90
!!  Contains code related to more than one scheme in the GFS physics suite.

    module GFS_suite_interstitial_rad_reset

    contains

    subroutine GFS_suite_interstitial_rad_reset_init ()
    end subroutine GFS_suite_interstitial_rad_reset_init

    subroutine GFS_suite_interstitial_rad_reset_finalize()
    end subroutine GFS_suite_interstitial_rad_reset_finalize

!> \section arg_table_GFS_suite_interstitial_rad_reset_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Interstitial   | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_rad_reset_run (Interstitial, errmsg, errflg)

      use GFS_typedefs, only: GFS_interstitial_type

      implicit none

      ! interface variables
      type(GFS_interstitial_type), intent(inout) :: Interstitial
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call Interstitial%rad_reset()

    end subroutine GFS_suite_interstitial_rad_reset_run

    end module GFS_suite_interstitial_rad_reset


    module GFS_suite_interstitial_phys_reset

    contains

    subroutine GFS_suite_interstitial_phys_reset_init ()
    end subroutine GFS_suite_interstitial_phys_reset_init

    subroutine GFS_suite_interstitial_phys_reset_finalize()
    end subroutine GFS_suite_interstitial_phys_reset_finalize

!> \section arg_table_GFS_suite_interstitial_phys_reset_run Argument Table
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Interstitial   | FV3-GFS_Interstitial_type                              | derived type GFS_interstitial_type in FV3               | DDT           |    0 | GFS_interstitial_type |           | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                   | flag          |    0 | integer               |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_phys_reset_run (Interstitial, errmsg, errflg)

      use GFS_typedefs, only: GFS_interstitial_type

      implicit none

      ! interface variables
      type(GFS_interstitial_type), intent(inout) :: Interstitial
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      errmsg = ''
      errflg = 0

      call Interstitial%phys_reset()

    end subroutine GFS_suite_interstitial_phys_reset_run

    end module GFS_suite_interstitial_phys_reset


    module GFS_suite_interstitial_1

    contains

    subroutine GFS_suite_interstitial_1_init ()
    end subroutine GFS_suite_interstitial_1_init

    subroutine GFS_suite_interstitial_1_finalize()
    end subroutine GFS_suite_interstitial_1_finalize

!> \section arg_table_GFS_suite_interstitial_1_run Argument Table
!! | local_name     | standard_name                                                            | long_name                                                               | units         | rank | type             |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------------------------|-------------------------------------------------------------------------|---------------|------|------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                                     | Fortran DDT containing FV3-GFS model control parameters                 | DDT           |    0 | GFS_control_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                                        | Fortran DDT containing FV3-GFS grid and interpolation related data      | DDT           |    0 | GFS_grid_type    |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                                     | Fortran DDT containing FV3-GFS surface fields                           | DDT           |    0 | GFS_sfcprop_type |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                                     | Fortran DDT containing FV3-GFS prognostic state data in from dycore     | DDT           |    0 | GFS_statein_type |           | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                                        | Fortran DDT containing FV3-GFS fields targeted for diagnostic output    | DDT           |    0 | GFS_diag_type    |           | inout  | F        |
!! | rhbbot         | critical_relative_humidity_at_surface                                    | critical relative humidity at the surface                               | frac          |    0 | real             | kind_phys | out    | F        |
!! | rhpbl          | critical_relative_humidity_at_PBL_top                                    | critical relative humidity at the PBL top                               | frac          |    0 | real             | kind_phys | out    | F        |
!! | rhbtop         | critical_relative_humidity_at_top_of_atmosphere                          | critical relative humidity at the top of atmosphere                     | frac          |    0 | real             | kind_phys | out    | F        |
!! | frain          | dynamics_to_physics_timestep_ratio                                       | ratio of dynamics timestep to physics timestep                          | none          |    0 | real             | kind_phys | out    | F        |
!! | islmsk         | sea_land_ice_mask                                                        | landmask: sea/land/ice=0/1/2                                            | flag          |    1 | integer          |           | out    | F        |
!! | work1          | grid_size_related_coefficient_used_in_scale-sensitive_schemes            | grid size related coefficient used in scale-sensitive schemes           | none          |    1 | real             | kind_phys | out    | F        |
!! | work2          | grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement | complement to work1                                                     | none          |    1 | real             | kind_phys | out    | F        |
!! | dudt           | tendency_of_x_wind_due_to_model_physics                                  | updated tendency of the x wind                                          | m s-2         |    2 | real             | kind_phys | out    | F        |
!! | dvdt           | tendency_of_y_wind_due_to_model_physics                                  | updated tendency of the y wind                                          | m s-2         |    2 | real             | kind_phys | out    | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics                         | updated tendency of the temperature                                     | K s-1         |    2 | real             | kind_phys | out    | F        |
!! | dtdtc          | tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky  | clear sky radiative (shortwave + longwave) heating rate at current time | K s-1         |    2 | real             | kind_phys | out    | F        |
!! | dqdt           | tendency_of_tracers_due_to_model_physics                                 | updated tendency of the tracers                                         | kg kg-1 s-1   |    3 | real             | kind_phys | out    | F        |
!! | errmsg         | error_message                                                            | error message for error handling in CCPP                                | none          |    0 | character        | len=*     | out    | F        |
!! | errflg         | error_flag                                                               | error flag for error handling in CCPP                                   | flag          |    0 | integer          |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_1_run (Model, Grid, Sfcprop, Statein, Diag, rhbbot, rhpbl, rhbtop, frain, islmsk, &
                                             work1, work2, dudt, dvdt, dtdt, dtdtc, dqdt, errmsg, errflg)

      use machine,               only: kind_phys
      use physcons,              only: dxmin, dxinv
      use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_sfcprop_type, GFS_statein_type, GFS_diag_type

      implicit none

      ! interface variables
      type(GFS_control_type),           intent(in) :: Model
      type(GFS_grid_type),              intent(in) :: Grid
      type(GFS_sfcprop_type),           intent(in) :: Sfcprop
      type(GFS_statein_type),           intent(in) :: Statein
      type(GFS_diag_type),              intent(inout) :: Diag

      real(kind=kind_phys), intent(out) :: rhbbot, rhpbl, rhbtop, frain
      integer, dimension(size(Grid%xlon,1)), intent(out) :: islmsk
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(out)  :: work1, work2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(out) :: dudt, dvdt, dtdt, dtdtc
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,Model%ntrac), intent(out) ::  dqdt
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      rhbbot = Model%crtrh(1)
      rhpbl  = Model%crtrh(2)
      rhbtop = Model%crtrh(3)

      frain = Model%dtf / Model%dtp

      do i = 1, size(Grid%xlon,1)
        islmsk(i)   = nint(Sfcprop%slmsk(i))
        work1(i)   = (log(Grid%area(i)) - dxmin) * dxinv
        work1(i)   = max(0.0, min(1.0,work1(i)))
        work2(i)   = 1.0 - work1(i)
        Diag%psurf(i)   = Statein%pgr(i)
      end do

      dudt(:,:)  = 0.
      dvdt(:,:)  = 0.
      dtdt(:,:)  = 0.
      dtdtc(:,:) = 0.
      dqdt(:,:,:) = 0.

    end subroutine GFS_suite_interstitial_1_run

  end module GFS_suite_interstitial_1


  module GFS_suite_interstitial_2

  contains

    subroutine GFS_suite_interstitial_2_init ()
    end subroutine GFS_suite_interstitial_2_init

    subroutine GFS_suite_interstitial_2_finalize()
    end subroutine GFS_suite_interstitial_2_finalize

!> \section arg_table_GFS_suite_interstitial_2_run Argument Table
!! | local_name     | standard_name                                                | long_name                                                             | units         | rank | type             |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                         | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_control_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                            | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_grid_type    |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                         | Fortran DDT containing FV3-GFS prognostic state data in from dycore   | DDT           |    0 | GFS_statein_type |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                         | Fortran DDT containing FV3-GFS radiation tendencies needed in physics | DDT           |    0 | GFS_radtend_type |           | in     | F        |
!! | xcosz          | instantaneous_cosine_of_zenith_angle                         | cosine of zenith angle at current time                                | none          |    1 | real             | kind_phys | in     | F        |
!! | adjsfcdsw      | surface_downwelling_shortwave_flux                           | surface downwelling shortwave flux at current time                    | W m-2         |    1 | real             | kind_phys | in     | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux                            | surface downwelling longwave flux at current time                     | W m-2         |    1 | real             | kind_phys | in     | F        |
!! | adjsfculw      | surface_upwelling_longwave_flux                              | surface upwelling longwave flux at current time                       | W m-2         |    1 | real             | kind_phys | in     | F        |
!! | xmu            | zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes | zenith angle temporal adjustment factor for shortwave fluxes          | none          |    1 | real             | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                            | Fortran DDT containing FV3-GFS fields targeted for diagnostic output  | DDT           |    0 | GFS_diag_type    |           | inout  | F        |
!! | kcnv           | flag_deep_convection                                         | flag indicating whether convection occurs in column (0 or 1)          | flag          |    1 | integer          |           | out    | F        |
!! | hflx           | kinematic_surface_upward_sensible_heat_flux                  | kinematic surface upward sensible heat flux                           | K m s-1       |    1 | real             | kind_phys | out    | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux                    | kinematic surface upward latent heat flux                             | kg kg-1 m s-1 |    1 | real             | kind_phys | out    | F        |
!! | errmsg         | error_message                                                | error message for error handling in CCPP                              | none          |    0 | character        | len=*     | out    | F        |
!! | errflg         | error_flag                                                   | error flag for error handling in CCPP                                 | flag          |    0 | integer          |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_2_run (Model, Grid, Statein, Radtend, xcosz, adjsfcdsw, adjsfcdlw, adjsfculw, xmu, &
                                             Diag, kcnv, hflx, evap, errmsg, errflg)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_statein_type, GFS_radtend_type, GFS_diag_type

      implicit none

      ! interface variables
      type(GFS_control_type),           intent(in)    :: Model
      type(GFS_grid_type),              intent(in)    :: Grid
      type(GFS_statein_type),           intent(in)    :: Statein
      type(GFS_radtend_type),           intent(in)    :: Radtend
      type(GFS_diag_type),              intent(inout) :: Diag

      integer, dimension(size(Grid%xlon,1)), intent(out) :: kcnv
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: xcosz, adjsfcdsw, adjsfcdlw, adjsfculw, xmu
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(out) :: hflx, evap
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! local variables
      real(kind=kind_phys), parameter :: czmin   = 0.0001      ! cos(89.994)
      integer :: i, k
      real(kind=kind_phys) :: tem1

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lssav) then      !  --- ...  accumulate/save output variables

!  --- ...  sunshine duration time is defined as the length of time (in mdl output
!           interval) that solar radiation falling on a plane perpendicular to the
!           direction of the sun >= 120 w/m2

        do i = 1, size(Grid%xlon,1)
          if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
            tem1 = adjsfcdsw(i) / xcosz(i)
            if ( tem1 >= 120.0 ) then
              Diag%suntim(i) = Diag%suntim(i) + Model%dtf
            endif
          endif
        enddo

!  --- ...  sfc lw fluxes used by atmospheric model are saved for output

        Diag%dlwsfc(:) = Diag%dlwsfc(:) +   adjsfcdlw(:)*Model%dtf
        Diag%ulwsfc(:) = Diag%ulwsfc(:) +   adjsfculw(:)*Model%dtf
        Diag%psmean(:) = Diag%psmean(:) + Statein%pgr(:)*Model%dtf        ! mean surface pressure

        if (Model%ldiag3d) then
          if (Model%lsidea) then
            Diag%dt3dt(:,:,1) = Diag%dt3dt(:,:,1) + Radtend%lwhd(:,:,1)*Model%dtf
            Diag%dt3dt(:,:,2) = Diag%dt3dt(:,:,2) + Radtend%lwhd(:,:,2)*Model%dtf
            Diag%dt3dt(:,:,3) = Diag%dt3dt(:,:,3) + Radtend%lwhd(:,:,3)*Model%dtf
            Diag%dt3dt(:,:,4) = Diag%dt3dt(:,:,4) + Radtend%lwhd(:,:,4)*Model%dtf
            Diag%dt3dt(:,:,5) = Diag%dt3dt(:,:,5) + Radtend%lwhd(:,:,5)*Model%dtf
            Diag%dt3dt(:,:,6) = Diag%dt3dt(:,:,6) + Radtend%lwhd(:,:,6)*Model%dtf
          else
            do k = 1, Model%levs
              Diag%dt3dt(:,k,1) = Diag%dt3dt(:,k,1) + Radtend%htrlw(:,k)*Model%dtf
              Diag%dt3dt(:,k,2) = Diag%dt3dt(:,k,2) + Radtend%htrsw(:,k)*Model%dtf*xmu(:)
            enddo
          endif
        endif
      endif    ! end if_lssav_block

      kcnv(:)   = 0

      hflx(:)       = 0.0
      evap(:)       = 0.0

      Diag%t1(:)      = Statein%tgrs(:,1)
      Diag%q1(:)      = Statein%qgrs(:,1,1)
      Diag%u1(:)      = Statein%ugrs(:,1)
      Diag%v1(:)      = Statein%vgrs(:,1)

    end subroutine GFS_suite_interstitial_2_run

  end module GFS_suite_interstitial_2


  module GFS_suite_update_stateout

  contains

    subroutine GFS_suite_update_stateout_init ()
    end subroutine GFS_suite_update_stateout_init

    subroutine GFS_suite_update_stateout_finalize()
    end subroutine GFS_suite_update_stateout_finalize

!> \section arg_table_GFS_suite_update_stateout_run Argument Table
!! | local_name     | standard_name                                                | long_name                                                             | units         | rank | type              |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | Statein        | FV3-GFS_Statein_type                                         | Fortran DDT containing FV3-GFS prognostic state data in from dycore   | DDT           |    0 | GFS_statein_type  |           | in     | F        |
!! | Model          | FV3-GFS_Control_type                                         | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_control_type  |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                            | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_grid_type     |           | in     | F        |
!! | dudt           | tendency_of_x_wind_due_to_model_physics                      | updated tendency of the x wind                                        | m s-2         |    2 | real              | kind_phys | in     | F        |
!! | dvdt           | tendency_of_y_wind_due_to_model_physics                      | updated tendency of the y wind                                        | m s-2         |    2 | real              | kind_phys | in     | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics             | updated tendency of the temperature                                   | K s-1         |    2 | real              | kind_phys | in     | F        |
!! | dqdt           | tendency_of_tracers_due_to_model_physics                     | updated tendency of the tracers                                       | kg kg-1 s-1   |    3 | real              | kind_phys | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                        | Fortran DDT containing FV3-GFS prognostic state to return to dycore   | DDT           |    0 | GFS_stateout_type |           | inout  | F        |
!! | errmsg         | error_message                                                | error message for error handling in CCPP                              | none          |    0 | character         | len=*     | out    | F        |
!! | errflg         | error_flag                                                   | error flag for error handling in CCPP                                 | flag          |    0 | integer           |           | out    | F        |
!!
    subroutine GFS_suite_update_stateout_run (Statein, Model, Grid, dudt, dvdt, dtdt, dqdt, Stateout, errmsg, errflg)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_control_type, GFS_statein_type, GFS_grid_type, GFS_stateout_type

      implicit none

      ! interface variables
      type(GFS_control_type),           intent(in)    :: Model
      type(GFS_statein_type),           intent(in)    :: Statein
      type(GFS_grid_type),              intent(in)    :: Grid
      type(GFS_stateout_type),          intent(inout) :: Stateout

      real(kind=kind_phys), dimension(size(Grid%xlon,1), Model%levs), intent(in) :: dudt, dvdt, dtdt
      real(kind=kind_phys), dimension(size(Grid%xlon,1), Model%levs, Model%ntrac), intent(in) :: dqdt

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      Stateout%gt0(:,:)   = Statein%tgrs(:,:) + dtdt(:,:) * Model%dtp
      Stateout%gu0(:,:)   = Statein%ugrs(:,:) + dudt(:,:) * Model%dtp
      Stateout%gv0(:,:)   = Statein%vgrs(:,:) + dvdt(:,:) * Model%dtp
      Stateout%gq0(:,:,:) = Statein%qgrs(:,:,:) + dqdt(:,:,:) * Model%dtp

    end subroutine GFS_suite_update_stateout_run

  end module GFS_suite_update_stateout


  module GFS_suite_interstitial_3

  contains

    subroutine GFS_suite_interstitial_3_init ()
    end subroutine GFS_suite_interstitial_3_init

    subroutine GFS_suite_interstitial_3_finalize()
    end subroutine GFS_suite_interstitial_3_finalize

!> \section arg_table_GFS_suite_interstitial_3_run Argument Table
!! | local_name     | standard_name                                                            | long_name                                                             | units         | rank | type             |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                                     | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_control_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                                        | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_grid_type    |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                                     | Fortran DDT containing FV3-GFS prognostic state data in from dycore   | DDT           |    0 | GFS_statein_type |           | in     | F        |
!! | rhbbot         | critical_relative_humidity_at_surface                                    | critical relative humidity at the surface                             | frac          |    0 | real             | kind_phys | in     | F        |
!! | rhbtop         | critical_relative_humidity_at_top_of_atmosphere                          | critical relative humidity at the top of atmosphere                   | frac          |    0 | real             | kind_phys | in     | F        |
!! | work1          | grid_size_related_coefficient_used_in_scale-sensitive_schemes            | grid size related coefficient used in scale-sensitive schemes         | none          |    1 | real             | kind_phys | in     | F        |
!! | work2          | grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement | complement to work1                                                   | none          |    1 | real             | kind_phys | in     | F        |
!! | clw            | convective_transportable_tracers                                         | array to contain cloud water and other convective trans. tracers      | kg kg-1       |    3 | real             | kind_phys | inout  | F        |
!! | cnvc           | convective_cloud_cover                                                   | convective cloud cover                                                | frac          |    2 | real             | kind_phys | inout  | F        |
!! | cnvw           | convective_cloud_water_specific_humidity                                 | convective cloud water specific humidity                              | kg kg-1       |    2 | real             | kind_phys | inout  | F        |
!! | ktop           | vertical_index_at_cloud_top                                              | vertical index at cloud top                                           | index         |    1 | integer          |           | inout  | F        |
!! | kbot           | vertical_index_at_cloud_base                                             | vertical index at cloud base                                          | index         |    1 | integer          |           | inout  | F        |
!! | rhc            | critical_relative_humidity                                               | critical relative humidity                                            | frac          |    2 | real             | kind_phys | out    | F        |
!! | errmsg         | error_message                                                            | error message for error handling in CCPP                              | none          |    0 | character        | len=*     | out    | F        |
!! | errflg         | error_flag                                                               | error flag for error handling in CCPP                                 | flag          |    0 | integer          |           | out    | F        |
!!
    subroutine GFS_suite_interstitial_3_run (Model, Grid, Statein, rhbbot, rhbtop, work1, work2, clw, cnvc, cnvw, &
                                             ktop, kbot, rhc, errmsg, errflg)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_statein_type
      use physcons,              only: rhc_max

      implicit none

      ! interface variables
      type(GFS_control_type),           intent(in)    :: Model
      type(GFS_grid_type),              intent(in)    :: Grid
      type(GFS_statein_type),           intent(in)    :: Statein

      real(kind=kind_phys), intent(in)                                           :: rhbbot, rhbtop
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)             :: work1, work2
      real(kind=kind_phys), intent(inout)                                        :: clw(:,:,:), cnvc(:,:), cnvw(:,:)
      integer, dimension(size(Grid%xlon,1)), intent(out)                         :: ktop, kbot
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(out) :: rhc

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! local variables
      integer :: i,k
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      clw(:,:,1) = 0.0
      clw(:,:,2) = -999.9
      if ((Model%imfdeepcnv >= 0) .or. (Model%imfshalcnv > 0)) then
        cnvc(:,:)  = 0.0
        cnvw(:,:)  = 0.0
      endif

      ktop(:)  = 1
      kbot(:)  = Model%levs
      rhc(:,:) = 0.0

      if (Model%ntcw > 0) then
        do k=1,Model%levs
          do i=1, size(Grid%xlon,1)
            tem      = rhbbot - (rhbbot-rhbtop) * (1.0-Statein%prslk(i,k))
            tem      = rhc_max * work1(i) + tem * work2(i)
            rhc(i,k) = max(0.0, min(1.0,tem))
          enddo
        enddo
      endif

    end subroutine GFS_suite_interstitial_3_run

  end module GFS_suite_interstitial_3
