!> \file GFS_SCNV_generic.f90
!!  Contains code related to shallow convective schemes to be used within the GFS physics suite.

      module GFS_SCNV_generic_pre

      contains

      subroutine GFS_SCNV_generic_pre_init ()
      end subroutine GFS_SCNV_generic_pre_init

      subroutine GFS_SCNV_generic_pre_finalize()
      end subroutine GFS_SCNV_generic_pre_finalize

!> \section arg_table_GFS_SCNV_generic_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                             | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | Fortran DDT containing FV3-GFS prognostic state to return to dycore   | DDT           |    0 | GFS_stateout_type     |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_grid_type         |           | in     | F        |
!! | save_t         | air_temperature_save                                   | air temperature before entering a physics scheme                      | K             |    2 | real                  | kind_phys | inout  | F        |
!! | save_qv        | water_vapor_specific_humidity_save                     | water vapor specific humidity before entering a physics scheme        | kg kg-1       |    2 | real                  | kind_phys | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                              | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                                 | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_SCNV_generic_pre_run (Model, Stateout, Grid, save_t, save_qv, errmsg, errflg)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_stateout_type, GFS_grid_type

        implicit none

        type(GFS_control_type),           intent(in)  :: Model
        type(GFS_stateout_type),          intent(in)  :: Stateout
        type(GFS_grid_type),              intent(in)  :: Grid
        real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(inout) :: save_t, save_qv
        character(len=*),                 intent(out) :: errmsg
        integer,                          intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (Model%ldiag3d) then
          save_t(:,:)   = Stateout%gt0(:,:)
        endif
        if (Model%ldiag3d .or. Model%lgocart) then
          save_qv(:,:) = Stateout%gq0(:,:,1)
        endif

    end subroutine GFS_SCNV_generic_pre_run

    end module GFS_SCNV_generic_pre

    module GFS_SCNV_generic_post

    contains

    subroutine GFS_SCNV_generic_post_init ()
    end subroutine GFS_SCNV_generic_post_init

    subroutine GFS_SCNV_generic_post_finalize ()
    end subroutine GFS_SCNV_generic_post_finalize

!> \section arg_table_GFS_SCNV_generic_post_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                             | units         | rank | type                  |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_control_type      |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | Fortran DDT containing FV3-GFS prognostic state to return to dycore   | DDT           |    0 | GFS_stateout_type     |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_grid_type         |           | in     | F        |
!! | save_t         | air_temperature_save                                   | air temperature before entering a physics scheme                      | K             |    2 | real                  | kind_phys | in     | F        |
!! | save_qv        | water_vapor_specific_humidity_save                     | water vapor specific humidity before entering a physics scheme        | kg kg-1       |    2 | real                  | kind_phys | in     | F        |
!! | frain          | dynamics_to_physics_timestep_ratio                     | ratio of dynamics timestep to physics timestep                        | none          |    0 | real                  | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                      | Fortran DDT containing FV3-GFS fields targeted for diagnostic output  | DDT           |    0 | GFS_diag_type         |           | inout  | F        |
!! | clw            | convective_transportable_tracers                       | array to contain cloud water and other convective trans. tracers      | kg kg-1       |    3 | real                  | kind_phys | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                              | none          |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                                 | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_SCNV_generic_post_run (Model, Stateout, Grid, save_t, save_qv, frain, Diag, clw, errmsg, errflg)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_grid_type, GFS_control_type, GFS_stateout_type, GFS_diag_type

      implicit none

      type(GFS_grid_type),            intent(in) :: Grid
      type(GFS_control_type),         intent(in) :: Model
      type(GFS_stateout_type),        intent(in) :: Stateout
      type(GFS_diag_type),         intent(inout) :: Diag

      real(kind=kind_phys),           intent(in) :: frain
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(in) :: save_t, save_qv

      real(kind=kind_phys),        intent(inout) :: clw(:,:,:)

      character(len=*),              intent(out) :: errmsg
      integer,                       intent(out) :: errflg

      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lssav) then
        if (Model%ldiag3d) then
          Diag%dt3dt(:,:,5) = Diag%dt3dt(:,:,5) + (Stateout%gt0(:,:)-save_t(:,:)) * frain
          Diag%dq3dt(:,:,3) = Diag%dq3dt(:,:,3) + (Stateout%gq0(:,:,1)-save_qv(:,:)) * frain
        endif
      endif   ! end if_lssav
!
      do k = 1, Model%levs
        do i = 1, size(Grid%xlon,1)
          if (clw(i,k,2) <= -999.0) clw(i,k,2) = 0.0
        enddo
      enddo

      end subroutine GFS_SCNV_generic_post_run

      end module GFS_SCNV_generic_post
