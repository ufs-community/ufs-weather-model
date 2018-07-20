!> \file GFS_DCNV_generic.f90
!!  Contains code related to deep convective schemes to be used within the GFS physics suite.

      module GFS_DCNV_generic_pre

      contains

      subroutine GFS_DCNV_generic_pre_init ()
      end subroutine GFS_DCNV_generic_pre_init

      subroutine GFS_DCNV_generic_pre_finalize()
      end subroutine GFS_DCNV_generic_pre_finalize

!> \section arg_table_GFS_DCNV_generic_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                                | units         | rank | type              |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                  | DDT           |    0 | GFS_control_type  |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | Fortran DDT containing FV3-GFS prognostic state to return to dycore      | DDT           |    0 | GFS_stateout_type |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | Fortran DDT containing FV3-GFS grid and interpolation related data       | DDT           |    0 | GFS_grid_type     |           | in     | F        |
!! | save_u         | x_wind_save                                            | x-wind before entering a physics scheme                                  | m s-1         |    2 | real              | kind_phys | inout  | F        |
!! | save_v         | y_wind_save                                            | y-wind before entering a physics scheme                                  | m s-1         |    2 | real              | kind_phys | inout  | F        |
!! | save_t         | air_temperature_save                                   | air temperature before entering a physics scheme                         | K             |    2 | real              | kind_phys | inout  | F        |
!! | save_qv        | water_vapor_specific_humidity_save                     | water vapor specific humidity before entering a physics scheme           | kg kg-1       |    2 | real              | kind_phys | inout  | F        |
!! | save_qcw       | cloud_condensed_water_specific_humidity_save           | cloud condensed water specific humidity before entering a physics scheme | kg kg-1       |    2 | real              | kind_phys | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                                 | none          |    0 | character         | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                                    | flag          |    0 | integer           |           | out    | F        |
!!
    subroutine GFS_DCNV_generic_pre_run (Model, Stateout, Grid, save_u, save_v, save_t, save_qv, save_qcw, errmsg, errflg)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_control_type, GFS_stateout_type, GFS_grid_type

      implicit none

      type(GFS_control_type),           intent(in) :: Model
      type(GFS_stateout_type),          intent(in) :: Stateout
      type(GFS_grid_type),              intent(in) :: Grid
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(inout) :: save_u
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(inout) :: save_v
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(inout) :: save_t
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(inout) :: save_qv
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(inout) :: save_qcw
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%ldiag3d) then
        save_t(:,:) = Stateout%gt0(:,:)
        save_u(:,:) = Stateout%gu0(:,:)
        save_v(:,:) = Stateout%gv0(:,:)
      elseif (Model%cnvgwd) then
        save_t(:,:) = Stateout%gt0(:,:)
      endif   ! end if_ldiag3d/cnvgwd

      if (Model%ldiag3d .or. Model%lgocart) then
        save_qv(:,:) = Stateout%gq0(:,:,1)
        save_qcw(:,:) = Stateout%gq0(:,:,3)
      endif   ! end if_ldiag3d/lgocart

    end subroutine GFS_DCNV_generic_pre_run

    end module GFS_DCNV_generic_pre

    module GFS_DCNV_generic_post

    contains

    subroutine GFS_DCNV_generic_post_init ()
    end subroutine GFS_DCNV_generic_post_init

    subroutine GFS_DCNV_generic_post_finalize ()
    end subroutine GFS_DCNV_generic_post_finalize

!> \section arg_table_GFS_DCNV_generic_post_run Argument Table
!! | local_name     | standard_name                                             | long_name                                                                | units         | rank | type              |    kind   | intent | optional |
!! |----------------|-----------------------------------------------------------|--------------------------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | Grid           | FV3-GFS_Grid_type                                         | Fortran DDT containing FV3-GFS grid and interpolation related data       | DDT           |    0 | GFS_grid_type     |           | in     | F        |
!! | Model          | FV3-GFS_Control_type                                      | Fortran DDT containing FV3-GFS model control parameters                  | DDT           |    0 | GFS_control_type  |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                     | Fortran DDT containing FV3-GFS prognostic state to return to dycore      | DDT           |    0 | GFS_stateout_type |           | in     | F        |
!! | frain          | dynamics_to_physics_timestep_ratio                        | ratio of dynamics timestep to physics timestep                           | none          |    0 | real              | kind_phys | in     | F        |
!! | rain1          | lwe_thickness_of_deep_convective_precipitation_amount     | deep convective rainfall amount on physics timestep                      | m             |    1 | real              | kind_phys | in     | F        |
!! | cld1d          | cloud_work_function                                       | cloud work function                                                      | m2 s-2        |    1 | real              | kind_phys | in     | F        |
!! | save_u         | x_wind_save                                               | x-wind before entering a physics scheme                                  | m s-1         |    2 | real              | kind_phys | in     | F        |
!! | save_v         | y_wind_save                                               | y-wind before entering a physics scheme                                  | m s-1         |    2 | real              | kind_phys | in     | F        |
!! | save_t         | air_temperature_save                                      | air temperature before entering a physics scheme                         | K             |    2 | real              | kind_phys | in     | F        |
!! | save_qv        | water_vapor_specific_humidity_save                        | water vapor specific humidity before entering a physics scheme           | kg kg-1       |    2 | real              | kind_phys | in     | F        |
!! | ud_mf          | instantaneous_atmosphere_updraft_convective_mass_flux     | (updraft mass flux) * delt                                               | kg m-2        |    2 | real              | kind_phys | in     | F        |
!! | dd_mf          | instantaneous_atmosphere_downdraft_convective_mass_flux   | (downdraft mass flux) * delt                                             | kg m-2        |    2 | real              | kind_phys | in     | F        |
!! | dt_mf          | instantaneous_atmosphere_detrainment_convective_mass_flux | (detrainment mass flux) * delt                                           | kg m-2        |    2 | real              | kind_phys | in     | F        |
!! | cnvw           | convective_cloud_water_specific_humidity                  | convective cloud water specific humidity                                 | kg kg-1       |    2 | real              | kind_phys | in     | F        |
!! | cnvc           | convective_cloud_cover                                    | convective cloud cover                                                   | frac          |    2 | real              | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                         | Fortran DDT containing FV3-GFS fields targeted for diagnostic output     | DDT           |    0 | GFS_diag_type     |           | inout  | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                          | Fortran DDT containing FV3-GFS miscellaneous data                        | DDT           |    0 | GFS_tbd_type      |           | inout  | F        |
!! | errmsg         | error_message                                             | error message for error handling in CCPP                                 | none          |    0 | character         | len=*     | out    | F        |
!! | errflg         | error_flag                                                | error flag for error handling in CCPP                                    | flag          |    0 | integer           |           | out    | F        |
!!
    subroutine GFS_DCNV_generic_post_run (Grid, Model, Stateout, frain, rain1, cld1d, save_u, save_v, save_t, save_qv, &
                                                               ud_mf, dd_mf, dt_mf, cnvw, cnvc, Diag, Tbd, errmsg, errflg)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_grid_type, GFS_control_type, GFS_stateout_type, GFS_diag_type, GFS_tbd_type
      use physcons,              only: con_g

      implicit none

      type(GFS_grid_type),            intent(in) :: Grid
      type(GFS_control_type),         intent(in) :: Model
      type(GFS_stateout_type),        intent(in) :: Stateout
      type(GFS_diag_type),         intent(inout) :: Diag
      type(GFS_tbd_type),          intent(inout) :: Tbd

      real(kind=kind_phys), intent(in) :: frain
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: rain1, cld1d
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(in) :: save_u, save_v, save_t, save_qv
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(in) :: ud_mf, dd_mf, dt_mf
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(in) :: cnvw, cnvc

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, num2, num3

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i = 1, size(Grid%xlon,1)
        Diag%rainc(:) = frain * rain1(:)
      enddo

      if (Model%lssav) then
        Diag%cldwrk (:) = Diag%cldwrk (:) + cld1d(:) * Model%dtf
        Diag%cnvprcp(:) = Diag%cnvprcp(:) + Diag%rainc(:)

        if (Model%ldiag3d) then
          Diag%dt3dt(:,:,4) = Diag%dt3dt(:,:,4) + (Stateout%gt0(:,:)-save_t(:,:)) * frain
          Diag%dq3dt(:,:,2) = Diag%dq3dt(:,:,2) + (Stateout%gq0(:,:,1)-save_qv(:,:)) * frain
          Diag%du3dt(:,:,3) = Diag%du3dt(:,:,3) + (Stateout%gu0(:,:)-save_u(:,:)) * frain
          Diag%dv3dt(:,:,3) = Diag%dv3dt(:,:,3) + (Stateout%gv0(:,:)-save_v(:,:)) * frain

          Diag%upd_mf(:,:)  = Diag%upd_mf(:,:)  + ud_mf(:,:) * (con_g*frain)
          Diag%dwn_mf(:,:)  = Diag%dwn_mf(:,:)  + dd_mf(:,:) * (con_g*frain)
          Diag%det_mf(:,:)  = Diag%det_mf(:,:)  + dt_mf(:,:) * (con_g*frain)
        endif ! if (ldiag3d)

      endif   ! end if_lssav

      if ((Model%npdf3d == 3) .and. (Model%num_p3d == 4)) then
        num2 = Model%num_p3d + 2
        num3 = num2 + 1
        Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
        Tbd%phy_f3d(:,:,num3) = cnvc(:,:)
      elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then
        num2 = Model%num_p3d + 1
        Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
      endif

    end subroutine GFS_DCNV_generic_post_run

    end module GFS_DCNV_generic_post
