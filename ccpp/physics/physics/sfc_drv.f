!>  \file sfc_drv.f
!!  This file contains the Noah land surface scheme driver.

      module lsm_noah_pre
      contains

      subroutine lsm_noah_pre_init
      end subroutine lsm_noah_pre_init

      subroutine lsm_noah_pre_finalize
      end subroutine lsm_noah_pre_finalize

!! \brief Brief description of the subroutine
!!
!! \section arg_table_lsm_noah_pre_run Argument Table
!! | local_name     | standard_name                                               | long_name                                  | units      | rank | type      |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------------|--------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                      | horizontal loop extent                     | count      |    0 | integer   |           | in     | F        |
!! | km             | soil_vertical_dimension                                     | soil vertical layer dimension              | count      |    0 | integer   |           | in     | F        |
!! | drain          | subsurface_runoff_flux                                      | subsurface runoff flux                     | g m-2 s-1  |    1 | real      | kind_phys | out    | F        |
!! | runof          | surface_runoff_flux                                         | surface runoff flux                        | g m-2 s-1  |    1 | real      | kind_phys | out    | F        |
!! | evbs           | soil_upward_latent_heat_flux                                | soil upward latent heat flux               | W m-2      |    1 | real      | kind_phys | out    | F        |
!! | evcw           | canopy_upward_latent_heat_flux                              | canopy upward latent heat flux             | W m-2      |    1 | real      | kind_phys | out    | F        |
!! | trans          | transpiration_flux                                          | total plant transpiration rate             | kg m-2 s-1 |    1 | real      | kind_phys | out    | F        |
!! | sbsno          | snow_deposition_sublimation_upward_latent_heat_flux         | latent heat flux from snow depo/subl       | W m-2      |    1 | real      | kind_phys | out    | F        |
!! | snowc          | surface_snow_area_fraction                                  | surface snow area fraction                 | frac       |    1 | real      | kind_phys | out    | F        |
!! | snohf          | snow_freezing_rain_upward_latent_heat_flux                  | latent heat flux due to snow and frz rain  | W m-2      |    1 | real      | kind_phys | out    | F        |
!! | smcwlt2        | volume_fraction_of_condensed_water_in_soil_at_wilting_point | soil water fraction at wilting point       | frac       |    1 | real      | kind_phys | out    | F        |
!! | smcref2        | threshold_volume_fraction_of_condensed_water_in_soil        | soil moisture threshold                    | frac       |    1 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                               | error message for error handling in CCPP   | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                  | error flag for error handling in CCPP      | flag       |    0 | integer   |           | out    | F        |
!!
      subroutine lsm_noah_pre_run                                       &
     &  (im,km,drain,runof,evbs,evcw,trans,sbsno,snowc,snohf,smcwlt2,   &
     &   smcref2,errmsg,errflg                                          &
     &  )

      use machine,           only: kind_phys

      implicit none

!  ---  interface variables
      integer, intent(in) :: im, km

      real(kind=kind_phys), dimension(im), intent(inout)  ::            &
     &    drain,runof,evbs,evcw,trans,sbsno,snowc,snohf,smcwlt2,smcref2

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      drain(:)   = 0.0
      runof(:)   = 0.0
      evbs(:)    = 0.0
      evcw(:)    = 0.0
      trans(:)   = 0.0
      sbsno(:)   = 0.0
      snowc(:)   = 0.0
      snohf(:)   = 0.0
      smcwlt2(:) = 0.0
      smcref2(:) = 0.0

      end subroutine lsm_noah_pre_run

      end module lsm_noah_pre


      module lsm_noah_post
      contains

      subroutine lsm_noah_post_init
      end subroutine lsm_noah_post_init

      subroutine lsm_noah_post_finalize
      end subroutine lsm_noah_post_finalize

!> \brief Brief description of the subroutine
!!
!! \section arg_table_lsm_noah_post_run Argument Table
!! | local_name     | standard_name                                               | long_name                                  | units      | rank | type      |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------------|--------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                      | horizontal loop extent                     | count      |    0 | integer   |           | in     | F        |
!! | km             | soil_vertical_dimension                                     | soil vertical layer dimension              | count      |    0 | integer   |           | in     | F        |
!! | flag_lssav     | flag_diagnostics                                            | flag for calculating diagnostic fields     | flag       |    0 | logical   |           | in     | F        |
!! | dtf            | time_step_for_dynamics                                      | dynamics time step                         | s          |    0 | real      | kind_phys | in     | F        |
!! | drain          | subsurface_runoff_flux                                      | subsurface runoff flux                     | g m-2 s-1  |    1 | real      | kind_phys | in     | F        |
!! | runof          | surface_runoff_flux                                         | surface runoff flux                        | g m-2 s-1  |    1 | real      | kind_phys | in     | F        |
!! | runoff         | total_runoff                                                | total runoff                               | kg m-2     |    1 | real      | kind_phys | inout  | F        |
!! | srunoff        | surface_runoff                                              | surface runoff                             | kg m-2     |    1 | real      | kind_phys | inout  | F        |
!! | errmsg         | error_message                                               | error message for error handling in CCPP   | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                  | error flag for error handling in CCPP      | flag       |    0 | integer   |           | out    | F        |
!!
!!  \section lsm_post_general General Algorithm
!!  \section lsm_post_detailed Detailed Algorithm
!!  @{

      subroutine lsm_noah_post_run                                      &
     &  (im,km, flag_lssav,dtf,drain,runof,runoff,srunoff,errmsg,errflg &
     &  )
      use machine,           only: kind_phys
      implicit none

!  ---  interface variables
      integer, intent(in) :: im, km

      logical, intent(in) :: flag_lssav
      real(kind=kind_phys), intent (in)   :: dtf

      real(kind=kind_phys), dimension(im), intent(in   )  ::            &
     &    drain, runof

      real(kind=kind_phys), dimension(im), intent(inout)  ::            &
     &    runoff, srunoff

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if(flag_lssav) then
        runoff(:)  = runoff(:)  + (drain(:)+runof(:)) * dtf * 0.001
        srunoff(:) = srunoff(:) + runof(:) * dtf * 0.001
      end if

      end subroutine lsm_noah_post_run

!! @}
      end module lsm_noah_post
!! @}

      module lsm_noah
      contains

      subroutine lsm_noah_init
      end subroutine lsm_noah_init

      subroutine lsm_noah_finalize
      end subroutine lsm_noah_finalize

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!      call sfc_drv                                                     !
!  ---  inputs:                                                         !
!          ( im, km, ps, u1, v1, t1, q1, soiltyp, vegtype, sigmaf,      !
!            sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          !
!            prsl1, prslki, zf, islimsk, ddvel, slopetyp,               !
!            shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      !
!            isot, ivegsrc,                                             !
!  ---  in/outs:                                                        !
!            weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        !
!            canopy, trans, tsurf, zorl,                                !
!  ---  outputs:                                                        !
!            sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      !
!            cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            !
!            smcwlt2, smcref2, wet1 )                                   !
!                                                                       !
!                                                                       !
!  subprogram called:  sflx                                             !
!                                                                       !
!  program history log:                                                 !
!         xxxx  --             created                                  !
!         200x  -- sarah lu    modified                                 !
!    oct  2006  -- h. wei      modified                                 !
!    apr  2009  -- y.-t. hou   modified to include surface emissivity   !
!                     effect on lw radiation. replaced the comfussing   !
!                     slrad (net sw + dlw) with sfc net sw snet=dsw-usw !
!    sep  2009  -- s. moorthi modification to remove rcl and unit change!
!    nov  2011  -- sarah lu    corrected wet1 calculation
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horiz dimention and num of used pts      1    !
!     km       - integer, vertical soil layer dimension            1    !
!     ps       - real, surface pressure (pa)                       im   !
!     u1, v1   - real, u/v component of surface layer wind         im   !
!     t1       - real, surface layer mean temperature (k)          im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     soiltyp  - integer, soil type (integer index)                im   !
!     vegtype  - integer, vegetation type (integer index)          im   !
!     sigmaf   - real, areal fractional cover of green vegetation  im   !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     dswflx   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     snet     - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     delt     - real, time interval (second)                      1    !
!     tg3      - real, deep soil temperature (k)                   im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, sfc layer 1 mean pressure (pa)              im   !
!     prslki   - real, dimensionless exner function at layer 1     im   !
!     zf       - real, height of bottom layer (m)                  im   !
!     islimsk  - integer, sea/land/ice mask (=0/1/2)               im   !
!     ddvel    - real,                                             im   !
!     slopetyp - integer, class of sfc slope (integer index)       im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, upper bound on max albedo over deep snow    im   !
!     sfalb    - real, mean sfc diffused sw albedo (fractional)    im   !
!     flag_iter- logical,                                          im   !
!     flag_guess-logical,                                          im   !
!     isot     - integer, sfc soil type data source zobler or statsgo   !
!     ivegsrc  - integer, sfc veg type data source umd or igbp          !
!                                                                       !
!  input/outputs:                                                       !
!     weasd    - real, water equivalent accumulated snow depth (mm) im  !
!     snwdph   - real, snow depth (water equiv) over land          im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     smc      - real, total soil moisture content (fractional)   im,km !
!     stc      - real, soil temp (k)                              im,km !
!     slc      - real, liquid soil moisture                       im,km !
!     canopy   - real, canopy moisture content (m)                 im   !
!     trans    - real, total plant transpiration (m/s)             im   !
!     tsurf    - real, surface skin temperature (after iteration)  im   !
!                                                                       !
!  outputs:                                                             !
!     sncovr1  - real, snow cover over land (fractional)           im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     drain    - real, subsurface runoff (mm/s)                    im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!     ep       - real, potential evaporation                       im   !
!     runoff   - real, surface runoff (m/s)                        im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     evbs     - real, direct soil evaporation (m/s)               im   !
!     evcw     - real, canopy water evaporation (m/s)              im   !
!     sbsno    - real, sublimation/deposit from snopack (m/s)      im   !
!     snowc    - real, fractional snow cover                       im   !
!     stm      - real, total soil column moisture content (m)      im   !
!     snohf    - real, snow/freezing-rain latent heat flux (w/m**2)im   !
!     smcwlt2  - real, dry soil moisture threshold                 im   !
!     smcref2  - real, soil moisture threshold                     im   !
!     zorl     - real, surface roughness                           im   !
!     wet1     - real, normalized soil wetness                     im   !
!                                                                       !
!  ====================    end of description    =====================  !

!-----------------------------------
!      subroutine sfc_drv                                                &
! \defgroup Noah_Main GFS Noah Land Surface Model 
!> \defgroup Noah_drv GFS Noah LSM Driver
!!  \brief This is Noah LSM driver module, with the functionality of 
!! preparing variables to run Noah LSM gfssflx(), calling Noah LSM and post-processing
!! variables for return to the parent model suite including unit conversion, as well 
!! as diagnotics calculation. 
!! \section arg_table_lsm_noah_run Argument Table
!! | local_name     | standard_name                                                                | long_name                                                       | units         | rank | type      |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                          | count         |    0 | integer   |           | in     | F        |
!! | km             | soil_vertical_dimension                                                      | soil vertical layer dimension                                   | count         |    0 | integer   |           | in     | F        |
!! | ps             | surface_air_pressure                                                         | surface pressure                                                | Pa            |    1 | real      | kind_phys | in     | F        |
!! | u1             | x_wind_at_lowest_model_layer                                                 | x component of 1st model layer wind                             | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1             | y_wind_at_lowest_model_layer                                                 | y component of 1st model layer wind                             | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature_at_lowest_model_layer                                        | 1st model layer air temperature                                 | K             |    1 | real      | kind_phys | in     | F        |
!! | q1             | specific_humidity_at_lowest_model_layer                                      | 1st model layer specific humidity                               | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | soiltyp        | cell_soil_type                                                               | soil type at each grid cell                                     | index         |    1 | integer   |           | in     | F        |
!! | vegtype        | cell_vegetation_type                                                         | vegetation type at each grid cell                               | index         |    1 | integer   |           | in     | F        |
!! | sigmaf         | vegetation_area_fraction                                                     | areal fractional cover of green vegetation                      | frac          |    1 | real      | kind_phys | in     | F        |
!! | sfcemis        | surface_longwave_emissivity                                                  | surface longwave emissivity                                     | frac          |    1 | real      | kind_phys | in     | F        |
!! | dlwflx         | surface_downwelling_longwave_flux_absorbed_by_ground                         | total sky surface downward longwave flux absorbed by the ground | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | dswsfc         | surface_downwelling_shortwave_flux                                           | total sky surface downward shortwave flux                       | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | snet           | surface_net_downwelling_shortwave_flux                                       | total sky surface net shortwave flux                            | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | delt           | time_step_for_dynamics                                                       | dynamics time step                                              | s             |    0 | real      | kind_phys | in     | F        |
!! | tg3            | deep_soil_temperature                                                        | bottom soil temperature                                         | K             |    1 | real      | kind_phys | in     | F        |
!! | cm             | surface_drag_coefficient_for_momentum_in_air                                 | surface exchange coeff for momentum                             | none          |    1 | real      | kind_phys | in     | F        |
!! | ch             | surface_drag_coefficient_for_heat_and_moisture_in_air                        | surface exchange coeff heat & moisture                          | none          |    1 | real      | kind_phys | in     | F        |
!! | prsl1          | air_pressure_at_lowest_model_layer                                           | Model layer 1 mean pressure                                     | Pa            |    1 | real      | kind_phys | in     | F        |
!! | prslki         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer     | ratio         |    1 | real      | kind_phys | in     | F        |
!! | zf             | height_above_mean_sea_level_at_lowest_model_layer                            | height above MSL at 1st model layer                             | m             |    1 | real      | kind_phys | in     | F        |
!! | islimsk        | sea_land_ice_mask                                                            | landmask: sea/land/ice=0/1/2                                    | flag          |    1 | integer   |           | in     | F        |
!! | ddvel          | surface_wind_enhancement_due_to_convection                                   | surface wind enhancement due to convection                      | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | slopetyp       | surface_slope_classification                                                 | class of sfc slope                                              | index         |    1 | integer   |           | in     | F        |
!! | shdmin         | minimum_vegetation_area_fraction                                             | min fractional coverage of green veg                            | frac          |    1 | real      | kind_phys | in     | F        |
!! | shdmax         | maximum_vegetation_area_fraction                                             | max fractnl cover of green veg (not used)                       | frac          |    1 | real      | kind_phys | in     | F        |
!! | snoalb         | upper_bound_on_max_albedo_over_deep_snow                                     | upper bound on max albedo over deep snow                        | frac          |    1 | real      | kind_phys | in     | F        |
!! | sfalb          | surface_diffused_shortwave_albedo                                            | mean surface diffused shortwave albedo                          | frac          |    1 | real      | kind_phys | in     | F        |
!! | flag_iter      | flag_for_iteration                                                           | flag for iteration                                              | flag          |    1 | logical   |           | in     | F        |
!! | flag_guess     | flag_for_guess_run                                                           | flag for guess run                                              | flag          |    1 | logical   |           | in     | F        |
!! | isot           | soil_type                                                                    | soil type (not used)                                            | index         |    0 | integer   |           | in     | F        |
!! | ivegsrc        | vegetation_type                                                              | vegetation type data source umd or igbp                         | index         |    0 | integer   |           | in     | F        |
!! | weasd          | water_equivalent_accumulated_snow_depth                                      | water equivalent accumulated snow depth                         | mm            |    1 | real      | kind_phys | inout  | F        |
!! | snwdph         | surface_snow_thickness_water_equivalent                                      | water equivalent snow depth over land                           | mm            |    1 | real      | kind_phys | inout  | F        |
!! | tskin          | surface_skin_temperature                                                     | surface skin temperature                                        | K             |    1 | real      | kind_phys | inout  | F        |
!! | tprcp          | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep       | nonnegative precipitation amount in one dynamics time step      | m             |    1 | real      | kind_phys | inout  | F        |
!! | srflag         | flag_for_precipitation_type                                                  | flag for snow or rain precipitation                             | flag          |    1 | real      | kind_phys | inout  | F        |
!! | smc            | volume_fraction_of_soil_moisture                                             | volumetric fraction of soil moisture                            | frac          |    2 | real      | kind_phys | inout  | F        |
!! | stc            | soil_temperature                                                             | soil temperature                                                | K             |    2 | real      | kind_phys | inout  | F        |
!! | slc            | volume_fraction_of_unfrozen_soil_moisture                                    | volume fraction of unfrozen soil moisture                       | frac          |    2 | real      | kind_phys | inout  | F        |
!! | canopy         | canopy_water_amount                                                          | canopy moisture content                                         | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | trans          | transpiration_flux                                                           | total plant transpiration rate                                  | kg m-2 s-1    |    1 | real      | kind_phys | inout  | F        |
!! | tsurf          | surface_skin_temperature_after_iteration                                     | surface skin temperature after iteration                        | K             |    1 | real      | kind_phys | inout  | F        |
!! | zorl           | surface_roughness_length                                                     | surface roughness length                                        | cm            |    1 | real      | kind_phys | inout  | F        |
!! | sncovr1        | surface_snow_area_fraction_for_diagnostics                                   | surface snow area fraction                                      | frac          |    1 | real      | kind_phys | inout  | F        |
!! | qsurf          | surface_specific_humidity                                                    | surface specific humidity                                       | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | gflux          | upward_heat_flux_in_soil                                                     | upward soil heat flux                                           | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | drain          | subsurface_runoff_flux                                                       | subsurface runoff flux                                          | g m-2 s-1     |    1 | real      | kind_phys | inout  | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux                                    | surface upward evaporation flux                                 | kg kg-1 m s-1 |    1 | real      | kind_phys | inout  | F        |
!! | hflx           | kinematic_surface_upward_sensible_heat_flux                                  | surface upward sensible heat flux                               | K m s-1       |    1 | real      | kind_phys | inout  | F        |
!! | ep             | surface_upward_potential_latent_heat_flux                                    | surface upward potential latent heat flux                       | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | runoff         | surface_runoff_flux                                                          | surface runoff flux                                             | g m-2 s-1     |    1 | real      | kind_phys | inout  | F        |
!! | cmm            | surface_drag_wind_speed_for_momentum_in_air                                  | surf mom exch coef time mean surf wind                          | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | chh            | surface_drag_mass_flux_for_heat_and_moisture_in_air                          | surf h&m exch coef time surf wind & density                     | kg m-2 s-1    |    1 | real      | kind_phys | inout  | F        |
!! | evbs           | soil_upward_latent_heat_flux                                                 | soil upward latent heat flux                                    | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | evcw           | canopy_upward_latent_heat_flux                                               | canopy upward latent heat flux                                  | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | sbsno          | snow_deposition_sublimation_upward_latent_heat_flux                          | latent heat flux from snow depo/subl                            | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | snowc          | surface_snow_area_fraction                                                   | surface snow area fraction                                      | frac          |    1 | real      | kind_phys | inout  | F        |
!! | stm            | soil_moisture_content                                                        | soil moisture content                                           | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | snohf          | snow_freezing_rain_upward_latent_heat_flux                                   | latent heat flux due to snow and frz rain                       | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | smcwlt2        | volume_fraction_of_condensed_water_in_soil_at_wilting_point                  | soil water fraction at wilting point                            | frac          |    1 | real      | kind_phys | inout  | F        |
!! | smcref2        | threshold_volume_fraction_of_condensed_water_in_soil                         | soil moisture threshold                                         | frac          |    1 | real      | kind_phys | inout  | F        |
!! | wet1           | normalized_soil_wetness                                                      | normalized soil wetness                                         | frac          |    1 | real      | kind_phys | inout  | F        |
!! | errmsg         | error_message                                                                | error message for error handling in CCPP                        | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                                   | error flag for error handling in CCPP                           | flag          |    0 | integer   |           | out    | F        |
!!
!!  \section general_noah_drv GFS Noah Driver General Algorithm
!!  @{
!  \section detailed_noah GFS Noah Driver Detailed Algorithm
!  @{
      subroutine lsm_noah_run                                            &
     &     ( im, km, ps, u1, v1, t1, q1, soiltyp, vegtype, sigmaf,      &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
     &       prsl1, prslki, zf, islimsk, ddvel, slopetyp,               &
     &       shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
     &       isot, ivegsrc,                                             & !  ---  inputs from here and above
     &       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
     &       canopy, trans, tsurf, zorl,                                & ! --- in/outs from here and above
     &       sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
     &       cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
     &       smcwlt2, smcref2, wet1, errmsg, errflg                     & ! -- outputs from here and above
     &     )

      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, only : grav   => con_g,    cp   => con_cp,          &
     &                     hvap   => con_hvap, rd   => con_rd,          &
     &                     eps    => con_eps, epsm1 => con_epsm1,       &
     &                     rvrdm1 => con_fvirt

      implicit none

!  ---  constant parameters:
      real(kind=kind_phys), parameter :: cpinv   = 1.0/cp
      real(kind=kind_phys), parameter :: hvapi   = 1.0/hvap
      real(kind=kind_phys), parameter :: elocp   = hvap/cp
      real(kind=kind_phys), parameter :: rhoh2o  = 1000.0
      real(kind=kind_phys), parameter :: a2      = 17.2693882
      real(kind=kind_phys), parameter :: a3      = 273.16
      real(kind=kind_phys), parameter :: a4      = 35.86
      real(kind=kind_phys), parameter :: a23m4   = a2*(a3-a4)

      real(kind=kind_phys), save         :: zsoil_noah(4)
      data zsoil_noah / -0.1, -0.4, -1.0, -2.0 /

!  ---  input:
      integer, intent(in) :: im, km, isot, ivegsrc

      integer, dimension(im), intent(in) :: soiltyp, vegtype, slopetyp

      real (kind=kind_phys), dimension(im), intent(in) :: ps, u1, v1,   &
     &       t1, q1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,    &
     &       ch, prsl1, prslki, ddvel, shdmin, shdmax,                  &
     &       snoalb, sfalb, zf

      integer, dimension(im), intent(in) :: islimsk
      real (kind=kind_phys),  intent(in) :: delt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess

!  ---  in/out:
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
     &       snwdph, tskin, tprcp, srflag, canopy, trans, tsurf, zorl

      real (kind=kind_phys), dimension(im,km), intent(inout) ::         &
     &       smc, stc, slc

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(inout) :: sncovr1,   &
     &       qsurf, gflux, drain, evap, hflx, ep, runoff, cmm, chh,     &
     &       evbs, evcw, sbsno, snowc, stm, snohf, smcwlt2, smcref2,    &
     &       wet1

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
     &       q0, qs1, theta1, wind, weasd_old, snwdph_old,              &
     &       tprcp_old, srflag_old, tskin_old, canopy_old

      real (kind=kind_phys), dimension(km) :: et, sldpth, stsoil,       &
     &       smsoil, slsoil

      real (kind=kind_phys), dimension(im,km) :: zsoil, smc_old,        &
     &       stc_old, slc_old

      real (kind=kind_phys) :: alb, albedo, beta, chx, cmx, cmc,        &
     &       dew, drip, dqsdt2, ec, edir, ett, eta, esnow, etp,         &
     &       flx1, flx2, flx3, ffrozp, lwdn, pc, prcp, ptu, q2,         &
     &       q2sat, solnet, rc, rcs, rct, rcq, rcsoil, rsmin,           &
     &       runoff1, runoff2, runoff3, sfcspd, sfcprs, sfctmp,         &
     &       sfcems, sheat, shdfac, shdmin1d, shdmax1d, smcwlt,         &
     &       smcdry, smcref, smcmax, sneqv, snoalb1d, snowh,            &
     &       snomlt, sncovr, soilw, soilm, ssoil, tsea, th2, tbot,      &
     &       xlai, zlvl, swdn, tem,z0

      integer :: couple, ice, nsoil, nroot, slope, stype, vtype
      integer :: i, k

      logical :: flag(im)
!
!===> ...  begin here
!

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!> - Set flag for land points.

      do i = 1, im
        flag(i) = (islimsk(i) == 1)
      enddo

!> - Save land-related prognostic fields for guess run.

      do i = 1, im
        if (flag(i) .and. flag_guess(i)) then
          weasd_old(i)  = weasd(i)
          snwdph_old(i) = snwdph(i)
          tskin_old(i)  = tskin(i)
          canopy_old(i) = canopy(i)
          tprcp_old(i)  = tprcp(i)
          srflag_old(i) = srflag(i)

          do k = 1, km
            smc_old(i,k) = smc(i,k)
            stc_old(i,k) = stc(i,k)
            slc_old(i,k) = slc(i,k)
          enddo
        endif
      enddo

!  --- ...  initialization block

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          ep(i)     = 0.0
          evap (i)  = 0.0
          hflx (i)  = 0.0
          gflux(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0
          snowc(i)  = 0.0
          snohf(i)  = 0.0
        endif
      enddo

!  --- ...  initialize variables

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          wind(i) = max(sqrt( u1(i)*u1(i) + v1(i)*v1(i) )               &
     &                + max(0.0, min(ddvel(i), 30.0)), 1.0)

          q0(i)   = max(q1(i), 1.e-8)   !* q1=specific humidity at level 1 (kg/kg)
          theta1(i) = t1(i) * prslki(i) !* adiabatic temp at level 1 (k)

          rho(i) = prsl1(i) / (rd*t1(i)*(1.0+rvrdm1*q0(i)))
          qs1(i) = fpvs( t1(i) )        !* qs1=sat. humidity at level 1 (kg/kg)
          qs1(i) = max(eps*qs1(i) / (prsl1(i)+epsm1*qs1(i)), 1.e-8)
          q0 (i) = min(qs1(i), q0(i))
        endif
      enddo

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          do k = 1, km
            zsoil(i,k) = zsoil_noah(k)
          enddo
        endif
      enddo

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then

!> - Prepare variables to run Noah LSM: 
!!  -   1. configuration information (c):
!!\n  ----------------------------------------
!!\n  \a couple  - couple-uncouple flag (=1: coupled, =0: uncoupled)
!!\n  \a ffrozp  - flag for snow-rain detection (1.=snow, 0.=rain)
!!\n  \a ice     - sea-ice flag (=1: sea-ice, =0: land)
!!\n  \a dt      - timestep (sec) (dt should not exceed 3600 secs) = delt
!!\n  \a zlvl    - height (\f$m\f$) above ground of atmospheric forcing variables
!!\n  \a nsoil   - number of soil layers (at least 2)
!!\n  \a sldpth  - the thickness of each soil layer (\f$m\f$)

          couple = 1                      ! run noah lsm in 'couple' mode

          if     (srflag(i) == 1.0) then  ! snow phase
            ffrozp = 1.0
          elseif (srflag(i) == 0.0) then  ! rain phase
            ffrozp = 0.0
          endif
          ice = 0

          zlvl = zf(i)

          nsoil = km
          sldpth(1) = - zsoil(i,1)
          do k = 2, km
            sldpth(k) = zsoil(i,k-1) - zsoil(i,k)
          enddo

!>  -   2. forcing data (f):
!!\n  ---------------------------------------
!!\n  \a lwdn    - lw dw radiation flux (\f$W m^{-2}\f$)
!!\n  \a solnet  - net sw radiation flux (dn-up) (\f$W m^{-2}\f$)
!!\n  \a sfcprs  - pressure at height zlvl above ground (pascals)
!!\n  \a prcp    - precip rate (\f$kg m^{-2} s^{-1}\f$)
!!\n  \a sfctmp  - air temperature (\f$K\f$) at height zlvl above ground
!!\n  \a th2     - air potential temperature (\f$K\f$) at height zlvl above ground
!!\n  \a q2      - mixing ratio at height zlvl above ground (\f$kg kg^{-1}\f$)

          lwdn   = dlwflx(i)         !..downward lw flux at sfc in w/m2
          swdn   = dswsfc(i)         !..downward sw flux at sfc in w/m2
          solnet = snet(i)           !..net sw rad flx (dn-up) at sfc in w/m2
          sfcems = sfcemis(i)

          sfcprs = prsl1(i)
          prcp   = rhoh2o * tprcp(i) / delt
          sfctmp = t1(i)
          th2    = theta1(i)
          q2     = q0(i)

!>  -   3. other forcing (input) data (i):
!!\n   ---------------------------------------
!!\n  \a sfcspd  - wind speed (\f$m s^{-1}\f$) at height zlvl above ground
!!\n  \a q2sat   - sat mixing ratio at height zlvl above ground (\f$kg kg^{-1}\f$)
!!\n  \a dqsdt2  - slope of sat specific humidity curve at t=sfctmp (\f$kg kg^{-1} k^{-1}\f$)

          sfcspd = wind(i)
          q2sat =  qs1(i)
          dqsdt2 = q2sat * a23m4/(sfctmp-a4)**2

!>  -   4. canopy/soil characteristics (s):
!!\n      ------------------------------------
!!\n \a vegtyp  - vegetation type (integer index)                   -> vtype
!!\n \a soiltyp - soil type (integer index)                         -> stype
!!\n \a slopetyp- class of sfc slope (integer index)                -> slope
!!\n \a shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!!\n \a shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
!!\n \a ptu     - photo thermal unit (plant phenology for annuals/crops)
!!\n \a alb     - backround snow-free surface albedo (fraction)
!!\n \a snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
!!\n \a tbot    - bottom soil temperature (local yearly-mean sfc air temp)

          vtype = vegtype(i)
          stype = soiltyp(i)
          slope = slopetyp(i)
          shdfac= sigmaf(i)

          shdmin1d = shdmin(i)
          shdmax1d = shdmax(i)
          snoalb1d = snoalb(i)

          ptu  = 0.0
          alb  = sfalb(i)
          tbot = tg3(i)

!>  -   5. history (state) variables (h):
!!\n      ------------------------------
!!\n \a cmc        - canopy moisture content (\f$m\f$)
!!\n \a t1         - ground/canopy/snowpack effective skin temperature (\f$K\f$)   -> tsea
!!\n \a stc(nsoil) - soil temp (\f$K\f$)                                         -> stsoil
!!\n \a smc(nsoil) - total soil moisture content (volumetric fraction)     -> smsoil
!!\n \a sh2o(nsoil)- unfrozen soil moisture content (volumetric fraction)  -> slsoil
!!\n \a snowh      - actual snow depth (\f$m\f$)
!!\n \a sneqv      - liquid water-equivalent snow depth (\f$m\f$)
!!\n \a albedo     - surface albedo including snow effect (unitless fraction)
!!\n \a ch         - surface exchange coefficient for heat and moisture (\f$m s^{-1}\f$) -> chx
!!\n \a cm         - surface exchange coefficient for momentum (\f$m s^{-1}\f$)          -> cmx
!!\n \a z0         - surface roughness (\f$m\f$)     -> zorl(\f$cm\f$)

          cmc = canopy(i) * 0.001            ! convert from mm to m
          tsea = tsurf(i)                    ! clu_q2m_iter

          do k = 1, km
            stsoil(k) = stc(i,k)
            smsoil(k) = smc(i,k)
            slsoil(k) = slc(i,k)
          enddo

          snowh = snwdph(i) * 0.001         ! convert from mm to m
          sneqv = weasd(i)  * 0.001         ! convert from mm to m
          if (sneqv /= 0.0 .and. snowh == 0.0) then
            snowh = 10.0 * sneqv
          endif

          chx    = ch(i)  * wind(i)              ! compute conductance
          cmx    = cm(i)  * wind(i)
          chh(i) = chx * rho(i)
          cmm(i) = cmx

!  ---- ... outside sflx, roughness uses cm as unit
          z0 = zorl(i)/100.

!> - Call Noah LSM gfssflx(). 

!          call sflx                                                     &
          call gfssflx                                                  & ! ccppdox: these is sflx in mpbl
!  ---  inputs:
     &     ( nsoil, couple, ice, ffrozp, delt, zlvl, sldpth,            &
     &       swdn, solnet, lwdn, sfcems, sfcprs, sfctmp,                &
     &       sfcspd, prcp, q2, q2sat, dqsdt2, th2, ivegsrc,             &
     &       vtype, stype, slope, shdmin1d, alb, snoalb1d,              &
!  ---  input/outputs:
     &       tbot, cmc, tsea, stsoil, smsoil, slsoil, sneqv, chx, cmx,  &
     &       z0,                                                        &
!  ---  outputs:
     &       nroot, shdfac, snowh, albedo, eta, sheat, ec,              &
     &       edir, et, ett, esnow, drip, dew, beta, etp, ssoil,         &
     &       flx1, flx2, flx3, runoff1, runoff2, runoff3,               &
     &       snomlt, sncovr, rc, pc, rsmin, xlai, rcs, rct, rcq,        &
     &       rcsoil, soilw, soilm, smcwlt, smcdry, smcref, smcmax)

!> - Noah LSM: prepare variables for return to parent model and unit conversion.
!>  -   6. output (o):
!!\n  ------------------------------
!!\n \a eta     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a sheat   - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a beta    - ratio of actual/potential evap (dimensionless)
!!\n \a etp     - potential evaporation (\f$W m^{-2}\f$)
!!\n \a ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!!\n \a runoff1 - surface runoff (\f$m s^{-1}\f$), not infiltrating the surface
!!\n \a runoff2 - subsurface runoff (\f$m s^{-1}\f$), drainage out bottom

          evap(i)  = eta
          hflx(i)  = sheat
          gflux(i) = ssoil

          evbs(i)  = edir
          evcw(i)  = ec
          trans(i) = ett
          sbsno(i) = esnow
          snowc(i) = sncovr
          stm(i)   = soilm
          snohf(i) = flx1 + flx2 + flx3

          smcwlt2(i) = smcwlt
          smcref2(i) = smcref

          ep(i)      = etp
          tsurf(i)   = tsea

          do k = 1, km
            stc(i,k) = stsoil(k)
            smc(i,k) = smsoil(k)
            slc(i,k) = slsoil(k)
          enddo
          wet1(i) = smsoil(1) / smcmax !Sarah Lu added 09/09/2010 (for GOCART)

!  --- ...  unit conversion (from m s-1 to mm s-1)
          runoff(i)  = runoff1 * 1000.0
          drain (i)  = runoff2 * 1000.0

!  --- ...  unit conversion (from m to mm)
          canopy(i)  = cmc   * 1000.0
          snwdph(i)  = snowh * 1000.0
          weasd(i)   = sneqv * 1000.0
          sncovr1(i) = sncovr
!  ---- ... outside sflx, roughness uses cm as unit (update after snow's
!  effect)
          zorl(i) = z0*100.

!  --- ...  do not return the following output fields to parent model
!    ec      - canopy water evaporation (m s-1)
!    edir    - direct soil evaporation (m s-1)
!    et(nsoil)-plant transpiration from a particular root layer (m s-1)
!    ett     - total plant transpiration (m s-1)
!    esnow   - sublimation from (or deposition to if <0) snowpack (m s-1)
!    drip    - through-fall of precip and/or dew in excess of canopy
!              water-holding capacity (m)
!    dew     - dewfall (or frostfall for t<273.15) (m)
!    beta    - ratio of actual/potential evap (dimensionless)
!    flx1    - precip-snow sfc (w m-2)
!    flx2    - freezing rain latent heat flux (w m-2)
!    flx3    - phase-change heat flux from snowmelt (w m-2)
!    snomlt  - snow melt (m) (water equivalent)
!    sncovr  - fractional snow cover (unitless fraction, 0-1)
!    runoff3 - numerical trunctation in excess of porosity (smcmax)
!              for a given soil layer at the end of a time step
!    rc      - canopy resistance (s m-1)
!    pc      - plant coefficient (unitless fraction, 0-1) where pc*etp
!              = actual transp
!    xlai    - leaf area index (dimensionless)
!    rsmin   - minimum canopy resistance (s m-1)
!    rcs     - incoming solar rc factor (dimensionless)
!    rct     - air temperature rc factor (dimensionless)
!    rcq     - atmos vapor pressure deficit rc factor (dimensionless)
!    rcsoil  - soil moisture rc factor (dimensionless)
!    soilw   - available soil moisture in root zone (unitless fraction
!              between smcwlt and smcmax)
!    soilm   - total soil column moisture content (frozen+unfrozen) (m)
!    smcwlt  - wilting point (volumetric)
!    smcdry  - dry soil moisture threshold where direct evap frm top
!              layer ends (volumetric)
!    smcref  - soil moisture threshold where transpiration begins to
!              stress (volumetric)
!    smcmax  - porosity, i.e. saturated value of soil moisture
!              (volumetric)
!    nroot   - number of root layers, a function of veg type, determined
!              in subroutine redprm.

        endif   ! end if_flag_iter_and_flag_block
      enddo   ! end do_i_loop

!> - Compute specific humidity at surface (\a qsurf).

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          rch(i)   = rho(i) * cp * ch(i) * wind(i)
          qsurf(i) = q1(i)  + evap(i) / (elocp * rch(i))
        endif
      enddo

!> - Compute surface upward sensible heat flux (\a hflx) and evaporation
!! flux (\a evap). 
      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          tem     = 1.0 / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif
      enddo

!> - Restore land-related prognostic fields for guess run.

      do i = 1, im
        if (flag(i)) then
          if (flag_guess(i)) then
            weasd(i)  = weasd_old(i)
            snwdph(i) = snwdph_old(i)
            tskin(i)  = tskin_old(i)
            canopy(i) = canopy_old(i)
            tprcp(i)  = tprcp_old(i)
            srflag(i) = srflag_old(i)

            do k = 1, km
              smc(i,k) = smc_old(i,k)
              stc(i,k) = stc_old(i,k)
              slc(i,k) = slc_old(i,k)
            enddo
          else
            tskin(i) = tsurf(i)
          endif
        endif
      enddo
!
      return
!...................................
!      end subroutine sfc_drv
      end subroutine lsm_noah_run
!-----------------------------------
!! @}

      end module lsm_noah
