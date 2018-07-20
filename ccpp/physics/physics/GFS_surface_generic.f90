!> \file GFS_surface_generic.f90
!!  Contains code related to all GFS surface schemes.

      module GFS_surface_generic_pre

      contains

      subroutine GFS_surface_generic_pre_init ()
      end subroutine GFS_surface_generic_pre_init

      subroutine GFS_surface_generic_pre_finalize()
      end subroutine GFS_surface_generic_pre_finalize

!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! | local_name     | standard_name                                                                | long_name                                                             | units      | rank | type                  |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|-----------------------------------------------------------------------|------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                                         | Fortran DDT containing FV3-GFS model control parameters               | DDT        |    0 | GFS_control_type      |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                                            | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT        |    0 | GFS_grid_type         |           | in     | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                                         | Fortran DDT containing FV3-GFS surface fields                         | DDT        |    0 | GFS_sfcprop_type      |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                                         | Fortran DDT containing FV3-GFS radiation tendencies needed in physics | DDT        |    0 | GFS_radtend_type      |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                                         | Fortran DDT containing FV3-GFS prognostic state data in from dycore   | DDT        |    0 | GFS_statein_type      |           | in     | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux                                            | surface downwelling longwave flux at current time                     | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                                            | Fortran DDT containing FV3-GFS fields targeted for diagnostic output  | DDT        |    0 | GFS_diag_type         |           | inout  | F        |
!! | sigmaf         | vegetation_area_fraction                                                     | areal fractional cover of green vegetation                            | frac       |    1 | real                  | kind_phys | inout  | F        |
!! | islmsk         | sea_land_ice_mask                                                            | landmask: sea/land/ice=0/1/2                                          | flag       |    1 | integer               |           | inout  | F        |
!! | soiltyp        | cell_soil_type                                                               | soil type at each grid cell                                           | index      |    1 | integer               |           | inout  | F        |
!! | vegtype        | cell_vegetation_type                                                         | vegetation type at each grid cell                                     | index      |    1 | integer               |           | inout  | F        |
!! | slopetyp       | surface_slope_classification                                                 | class of sfc slope                                                    | index      |    1 | integer               |           | inout  | F        |
!! | work3          | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer           | ratio      |    1 | real                  | kind_phys | inout  | F        |
!! | gabsbdlw       | surface_downwelling_longwave_flux_absorbed_by_ground                         | total sky surface downward longwave flux absorbed by the ground       | W m-2      |    1 | real                  | kind_phys | inout  | F        |
!! | tsurf          | surface_skin_temperature_after_iteration                                     | surface skin temperature after iteration                              | K          |    1 | real                  | kind_phys | inout  | F        |
!! | flag_guess     | flag_for_guess_run                                                           | flag for guess run                                                    | flag       |    1 | logical               |           | inout  | F        |
!! | flag_iter      | flag_for_iteration                                                           | flag for iteration                                                    | flag       |    1 | logical               |           | inout  | F        |
!! | ep1d           | surface_upward_potential_latent_heat_flux                                    | surface upward potential latent heat flux                             | W m-2      |    1 | real                  | kind_phys | inout  | F        |
!! | errmsg         | error_message                                                                | error message for error handling in CCPP                              | none       |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                                                   | error flag for error handling in CCPP                                 | flag       |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_surface_generic_pre_run (Model, Grid, Sfcprop, Radtend, Statein, adjsfcdlw, Diag, sigmaf, islmsk, &
                        soiltyp, vegtype, slopetyp, work3, gabsbdlw, tsurf, flag_guess, flag_iter, ep1d, errmsg, errflg)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_sfcprop_type, &
                                         GFS_diag_type, GFS_radtend_type, GFS_statein_type
        use physcons,              only: con_g

        implicit none

        type(GFS_grid_type),              intent(in) :: Grid
        type(GFS_control_type),           intent(in) :: Model
        type(GFS_sfcprop_type),           intent(in) :: Sfcprop
        type(GFS_radtend_type),           intent(in) :: Radtend
        type(GFS_statein_type),           intent(in) :: Statein
        type(GFS_diag_type),              intent(inout) :: Diag

        integer, dimension(size(Grid%xlon,1)), intent(inout) :: islmsk, soiltyp, vegtype, slopetyp

        logical, dimension(size(Grid%xlon,1)), intent(inout) :: flag_guess, flag_iter

        real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)  :: adjsfcdlw
        real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(inout)  :: sigmaf, work3, gabsbdlw, tsurf, ep1d

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        integer :: i
        real(kind=kind_phys), parameter :: onebg   = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        do i = 1, size(Grid%xlon,1)
            sigmaf(i)   = max( Sfcprop%vfrac(i),0.01 )

            if (islmsk(i) == 2) then
              if (Model%isot == 1) then
                soiltyp(i)  = 16
              else
                soiltyp(i)  = 9
              endif
              if (Model%ivegsrc == 1) then
                vegtype(i)  = 15
              elseif(Model%ivegsrc == 2) then
                vegtype(i)  = 13
              endif
              slopetyp(i) = 9
            else
              soiltyp(i)  = int( Sfcprop%stype(i)+0.5 )
              vegtype(i)  = int( Sfcprop%vtype(i)+0.5 )
              slopetyp(i) = int( Sfcprop%slope(i)+0.5 )    !! clu: slope -> slopetyp
            endif

            work3(i)   = Statein%prsik(i,1) / Statein%prslk(i,1)
        end do

        !  ---  convert lw fluxes for land/ocean/sea-ice models
        !  note: for sw: adjsfcdsw and adjsfcnsw are zenith angle adjusted downward/net fluxes.
        !        for lw: adjsfcdlw is (sfc temp adjusted) downward fluxe with no emiss effect.
        !                adjsfculw is (sfc temp adjusted) upward fluxe including emiss effect.
        !        one needs to be aware that that the absorbed downward lw flux (used by land/ocean
        !        models as downward flux) is not the same as adjsfcdlw but a value reduced by
        !        the factor of emissivity.  however, the net effects are the same when seeing
        !        it either above the surface interface or below.
        !
        !   - flux above the interface used by atmosphere model:
        !        down: adjsfcdlw;    up: adjsfculw = sfcemis*sigma*T**4 + (1-sfcemis)*adjsfcdlw
        !        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)
        !   - flux below the interface used by lnd/oc/ice models:
        !        down: sfcemis*adjsfcdlw;  up: sfcemis*sigma*T**4
        !        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)

        !  --- ...  define the downward lw flux absorbed by ground

        gabsbdlw(:) = Radtend%semis(:) * adjsfcdlw(:)

        tsurf(:)      = Sfcprop%tsfc(:)
        flag_guess(:) = .false.
        flag_iter(:)  = .true.

        ep1d(:)       = 0.0

        Diag%zlvl(:)    = Statein%phil(:,1) * onebg

      end subroutine GFS_surface_generic_pre_run

      end module GFS_surface_generic_pre

      module GFS_surface_generic_post

      contains

      subroutine GFS_surface_generic_post_init ()
      end subroutine GFS_surface_generic_post_init

      subroutine GFS_surface_generic_post_finalize()
      end subroutine GFS_surface_generic_post_finalize

!> \section arg_table_GFS_surface_generic_post_run Argument Table
!! | local_name     | standard_name                                       | long_name                                                            | units      | rank | type                  |    kind   | intent | optional |
!! |----------------|-----------------------------------------------------|----------------------------------------------------------------------|------------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                | Fortran DDT containing FV3-GFS model control parameters              | DDT        |    0 | GFS_control_type      |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                   | Fortran DDT containing FV3-GFS grid and interpolation related data   | DDT        |    0 | GFS_grid_type         |           | in     | F        |
!! | ep1d           | surface_upward_potential_latent_heat_flux           | surface upward potential latent heat flux                            | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | gflx           | upward_heat_flux_in_soil                            | upward soil heat flux                                                | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | evbs           | soil_upward_latent_heat_flux                        | soil upward latent heat flux                                         | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | evcw           | canopy_upward_latent_heat_flux                      | canopy upward latent heat flux                                       | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | trans          | transpiration_flux                                  | total plant transpiration rate                                       | kg m-2 s-1 |    1 | real                  | kind_phys | in     | F        |
!! | sbsno          | snow_deposition_sublimation_upward_latent_heat_flux | latent heat flux from snow depo/subl                                 | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | snowc          | surface_snow_area_fraction                          | surface snow area fraction                                           | frac       |    1 | real                  | kind_phys | in     | F        |
!! | snohf          | snow_freezing_rain_upward_latent_heat_flux          | latent heat flux due to snow and frz rain                            | W m-2      |    1 | real                  | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                   | Fortran DDT containing FV3-GFS fields targeted for diagnostic output | DDT        |    0 | GFS_diag_type         |           | inout  | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                | Fortran DDT containing FV3-GFS surface fields                        | DDT        |    0 | GFS_sfcprop_type      |           | inout  | F        |
!! | errmsg         | error_message                                       | error message for error handling in CCPP                             | none       |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                                          | error flag for error handling in CCPP                                | flag       |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_surface_generic_post_run (Model, Grid, ep1d, gflx, evbs, evcw, trans, sbsno, snowc, snohf, Diag, &
                                               Sfcprop, errmsg, errflg)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, GFS_sfcprop_type, GFS_diag_type

        implicit none

        type(GFS_control_type),           intent(in) :: Model
        type(GFS_grid_type),              intent(in) :: Grid
        type(GFS_sfcprop_type),           intent(inout) :: Sfcprop
        type(GFS_diag_type),              intent(inout) :: Diag

        real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in)  :: ep1d, gflx, evbs, evcw, trans, sbsno, snowc, snohf

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        Diag%epi(:)     = ep1d(:)
        Diag%gfluxi(:)  = gflx(:)

        if (Model%lssav) then
          Diag%gflux(:)   = Diag%gflux(:)  + gflx(:)  * Model%dtf
          Diag%evbsa(:)   = Diag%evbsa(:)  + evbs(:)  * Model%dtf
          Diag%evcwa(:)   = Diag%evcwa(:)  + evcw(:)  * Model%dtf
          Diag%transa(:)  = Diag%transa(:) + trans(:) * Model%dtf
          Diag%sbsnoa(:)  = Diag%sbsnoa(:) + sbsno(:) * Model%dtf
          Diag%snowca(:)  = Diag%snowca(:) + snowc(:) * Model%dtf
          Diag%snohfa(:)  = Diag%snohfa(:) + snohf(:) * Model%dtf
          Diag%ep(:)      = Diag%ep(:)     + ep1d(:)  * Model%dtf

          Diag%tmpmax(:)  = max(Diag%tmpmax(:),Sfcprop%t2m(:))
          Diag%tmpmin(:)  = min(Diag%tmpmin(:),Sfcprop%t2m(:))

          Diag%spfhmax(:) = max(Diag%spfhmax(:),Sfcprop%q2m(:))
          Diag%spfhmin(:) = min(Diag%spfhmin(:),Sfcprop%q2m(:))
        endif

      end subroutine GFS_surface_generic_post_run

      end module GFS_surface_generic_post
