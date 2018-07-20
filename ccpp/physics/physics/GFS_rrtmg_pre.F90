!> \file GFS_rrtmg_pre.f90
!! This file contains
      module GFS_rrtmg_pre

      public GFS_rrtmg_pre_run

      contains

!> \defgroup GFS_rrtmg_pre GFS RRTMG Scheme Pre
!! @{
!! \section arg_table_GFS_rrtmg_pre_init Argument Table
!!
      subroutine GFS_rrtmg_pre_init ()
      end subroutine GFS_rrtmg_pre_init

!> \section arg_table_GFS_rrtmg_pre_run Argument Table
!! | local_name        | standard_name                                                 | long_name                                                                     | units    | rank |  type            |   kind    | intent | optional |
!! |-------------------|---------------------------------------------------------------|-------------------------------------------------------------------------------|----------|------|------------------|-----------|--------|----------|
!! | Model             | FV3-GFS_Control_type                                          | Fortran DDT containing FV3-GFS model control parameters                       | DDT      |    0 | GFS_control_type |           | in     | F        |
!! | Grid              | FV3-GFS_Grid_type                                             | Fortran DDT containing FV3-GFS grid and interpolation related data            | DDT      |    0 | GFS_grid_type    |           | in     | F        |
!! | Sfcprop           | FV3-GFS_Sfcprop_type                                          | Fortran DDT containing FV3-GFS surface fields                                 | DDT      |    0 | GFS_sfcprop_type |           | in     | F        |
!! | Statein           | FV3-GFS_Statein_type                                          | Fortran DDT containing FV3-GFS prognostic state data in from dycore           | DDT      |    0 | GFS_statein_type |           | in     | F        |
!! | Tbd               | FV3-GFS_Tbd_type                                              | Fortran DDT containing FV3-GFS data not yet assigned to a defined container   | DDT      |    0 | GFS_tbd_type     |           | in     | F        |
!! | Cldprop           | FV3-GFS_Cldprop_type                                          | Fortran DDT containing FV3-GFS cloud fields needed by radiation from physics  | DDT      |    0 | GFS_cldprop_type |           | in     | F        |
!! | Radtend           | FV3-GFS_Radtend_type                                          | Fortran DDT containing FV3-GFS radiation tendencies                           | DDT      |    0 | GFS_radtend_type |           | inout  | F        |
!! | lm                | vertical_layer_dimension_for_radiation                        | number of vertical layers for radiation calculation                           | count    |    0 | integer          |           | in     | F        |
!! | im                | horizontal_loop_extent                                        | horizontal loop extent                                                        | count    |    0 | integer          |           | in     | F        |
!! | lmk               | adjusted_vertical_layer_dimension_for_radiation               | number of vertical layers for radiation                                       | count    |    0 | integer          |           | in     | F        |
!! | lmp               | adjusted_vertical_level_dimension_for_radiation               | number of vertical levels for radiation                                       | count    |    0 | integer          |           | in     | F        |
!! | kd                | vertical_index_difference_between_inout_and_local             | vertical index difference between in/out and local                            | index    |    0 | integer          |           | out    | F        |
!! | kt                | vertical_index_difference_between_layer_and_upper_bound       | vertical index difference between layer and upper bound                       | index    |    0 | integer          |           | out    | F        |
!! | kb                | vertical_index_difference_between_layer_and_lower_bound       | vertical index difference between layer and lower bound                       | index    |    0 | integer          |           | out    | F        |
!! | raddt             | time_step_for_radiation                                       | radiation time step                                                           | s        |    0 | real             | kind_phys | out    | F        |
!! | plvl              | air_pressure_at_interface_for_radiation_in_hPa                | air pressure at vertical interface for radiation calculation                  | hPa      |    2 | real             | kind_phys | out    | F        |
!! | plyr              | air_pressure_at_layer_for_radiation_in_hPa                    | air pressure at vertical layer for radiation calculation                      | hPa      |    2 | real             | kind_phys | out    | F        |
!! | tlvl              | air_temperature_at_interface_for_radiation                    | air temperature at vertical interface for radiation calculation               | K        |    2 | real             | kind_phys | out    | F        |
!! | tlyr              | air_temperature_at_layer_for_radiation                        | air temperature at vertical layer for radiation calculation                   | K        |    2 | real             | kind_phys | out    | F        |
!! | tsfg              | surface_ground_temperature_for_radiation                      | surface ground temperature for radiation                                      | K        |    1 | real             | kind_phys | out    | F        |
!! | tsfa              | surface_air_temperature_for_radiation                         | lowest model layer air temperature for radiation                              | K        |    1 | real             | kind_phys | out    | F        |
!! | qlyr              | water_vapor_specific_humidity_at_layer_for_radiation          | water vapor specific humidity at vertical layer for radiation calculation     | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | olyr              | ozone_concentration_at_layer_for_radiation                    | ozone concentration                                                           | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_co2        | volume_mixing_ratio_co2                                       | CO2 volume mixing ratio                                                       | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_n2o        | volume_mixing_ratio_n2o                                       | N2O volume mixing ratio                                                       | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_ch4        | volume_mixing_ratio_ch4                                       | CH4 volume mixing ratio                                                       | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_o2         | volume_mixing_ratio_o2                                        | O2 volume mixing ratio                                                        | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_co         | volume_mixing_ratio_co                                        | CO volume mixing ratio                                                        | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_cfc11      | volume_mixing_ratio_cfc11                                     | CFC11 volume mixing ratio                                                     | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_cfc12      | volume_mixing_ratio_cfc12                                     | CFC12 volume mixing ratio                                                     | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_cfc22      | volume_mixing_ratio_cfc22                                     | CFC22 volume mixing ratio                                                     | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_ccl4       | volume_mixing_ratio_ccl4                                      | CCL4 volume mixing ratio                                                      | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | gasvmr_cfc113     | volume_mixing_ratio_cfc113                                    | CFC113 volume mixing ratio                                                    | kg kg-1  |    2 | real             | kind_phys | out    | F        |
!! | faersw1           | aerosol_optical_depth_for_shortwave_bands_01-16               | aerosol optical depth for shortwave bands 01-16                               | none     |    3 | real             | kind_phys | out    | F        |
!! | faersw2           | aerosol_single_scattering_albedo_for_shortwave_bands_01-16    | aerosol single scattering albedo for shortwave bands 01-16                    | frac     |    3 | real             | kind_phys | out    | F        |
!! | faersw3           | aerosol_asymmetry_parameter_for_shortwave_bands_01-16         | aerosol asymmetry parameter for shortwave bands 01-16                         | none     |    3 | real             | kind_phys | out    | F        |
!! | faerlw1           | aerosol_optical_depth_for_longwave_bands_01-16                | aerosol optical depth for longwave bands 01-16                                | none     |    3 | real             | kind_phys | out    | F        |
!! | faerlw2           | aerosol_single_scattering_albedo_for_longwave_bands_01-16     | aerosol single scattering albedo for longwave bands 01-16                     | frac     |    3 | real             | kind_phys | out    | F        |
!! | faerlw3           | aerosol_asymmetry_parameter_for_longwave_bands_01-16          | aerosol asymmetry parameter for longwave bands 01-16                          | none     |    3 | real             | kind_phys | out    | F        |
!! | aerodp            | atmosphere_optical_thickness_due_to_ambient_aerosol_particles | vertical integrated optical depth for various aerosol species                 | none     |    2 | real             | kind_phys | out    | F        |
!! | clouds1           | total_cloud_fraction                                          | layer total cloud fraction                                                    | frac     |    2 | real             | kind_phys | out    | F        |
!! | clouds2           | cloud_liquid_water_path                                       | layer cloud liquid water path                                                 | g m-2    |    2 | real             | kind_phys | out    | F        |
!! | clouds3           | mean_effective_radius_for_liquid_cloud                        | mean effective radius for liquid cloud                                        | micron   |    2 | real             | kind_phys | out    | F        |
!! | clouds4           | cloud_ice_water_path                                          | layer cloud ice water path                                                    | g m-2    |    2 | real             | kind_phys | out    | F        |
!! | clouds5           | mean_effective_radius_for_ice_cloud                           | mean effective radius for ice cloud                                           | micron   |    2 | real             | kind_phys | out    | F        |
!! | clouds6           | cloud_rain_water_path                                         | cloud rain water path                                                         | g m-2    |    2 | real             | kind_phys | out    | F        |
!! | clouds7           | mean_effective_radius_for_rain_drop                           | mean effective radius for rain drop                                           | micron   |    2 | real             | kind_phys | out    | F        |
!! | clouds8           | cloud_snow_water_path                                         | cloud snow water path                                                         | g m-2    |    2 | real             | kind_phys | out    | F        |
!! | clouds9           | mean_effective_radius_for_snow_flake                          | mean effective radius for snow flake                                          | micron   |    2 | real             | kind_phys | out    | F        |
!! | cldsa             | cloud_area_fraction_for_radiation                             | fraction of clouds for low, middle,high, total and BL                         | frac     |    2 | real             | kind_phys | out    | F        |
!! | mtopa             | model_layer_number_at_cloud_top                               | vertical indices for low, middle and high cloud tops                          | index    |    2 | integer          |           | out    | F        |
!! | mbota             | model_layer_number_at_cloud_base                              | vertical indices for low, middle and high cloud bases                         | index    |    2 | integer          |           | out    | F        |
!! | errmsg            | error_message                                                 | error message for error handling in CCPP                                      | none     |    0 | character        | len=*     | out    | F        |
!! | errflg            | error_flag                                                    | error flag for error handling in CCPP                                         | flag     |    0 | integer          |           | out    | F        |
!!
      ! Attention - the output arguments lm, im, lmk, lmp must not be set
      ! in the CCPP version - they are defined in the interstitial_create routine
      subroutine GFS_rrtmg_pre_run (Model, Grid, Sfcprop, Statein,   & ! input
          Tbd, Cldprop,                                              &
          Radtend,                                                   & ! input/output
          lm, im, lmk, lmp, kd, kt, kb, raddt, plvl, plyr,           & ! output
          tlvl, tlyr, tsfg, tsfa, qlyr, olyr,                        &
          gasvmr_co2,   gasvmr_n2o,   gasvmr_ch4,   gasvmr_o2,       &
          gasvmr_co,    gasvmr_cfc11, gasvmr_cfc12,                  &
          gasvmr_cfc22, gasvmr_ccl4,  gasvmr_cfc113,                 &
          faersw1,  faersw2,  faersw3,                               &
          faerlw1, faerlw2, faerlw3, aerodp,                         &
          clouds1, clouds2, clouds3, clouds4, clouds5, clouds6,      &
          clouds7, clouds8, clouds9, cldsa, mtopa, mbota, errmsg, errflg)

      use machine,                   only: kind_phys
      use GFS_typedefs,              only: GFS_statein_type,   &
                                           GFS_stateout_type,  &
                                           GFS_sfcprop_type,   &
                                           GFS_coupling_type,  &
                                           GFS_control_type,   &
                                           GFS_grid_type,      &
                                           GFS_tbd_type,       &
                                           GFS_cldprop_type,   &
                                           GFS_radtend_type,   &
                                           GFS_diag_type
      use physparam
      use physcons,                  only: eps   => con_eps,         &
     &                                     epsm1 => con_epsm1,       &
     &                                     fvirt => con_fvirt        &
     &,                                    rocp  => con_rocp
      use radcons,                   only: itsfc,ltp, lextop, qmin,  &
                                           qme5, qme6, epsq, prsmin
      use funcphys,                  only: fpvs

      use module_radiation_astronomy,only: coszmn                         ! sol_init, sol_update
      use module_radiation_gases,    only: NF_VGAS, getgases, getozn      ! gas_init, gas_update,
      use module_radiation_aerosols, only: NF_AESW, NF_AELW, setaer,   &  ! aer_init, aer_update,
     &                                     NSPC1
      use module_radiation_clouds,   only: NF_CLDS,                    &  ! cld_init
     &                                     progcld1, progcld2,progcld3,&
     &                                     progclduni, diagcld1
      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type,   &
     &                                     profsw_type, NBDSW
      use module_radlw_parameters,   only: topflw_type, sfcflw_type,    &
     &                                     proflw_type, NBDLW

      implicit none

      type(GFS_control_type),              intent(in)    :: Model
      type(GFS_grid_type),                 intent(in)    :: Grid
      type(GFS_sfcprop_type),              intent(in)    :: Sfcprop
      type(GFS_statein_type),              intent(in)    :: Statein
      type(GFS_radtend_type),              intent(inout) :: Radtend
      type(GFS_tbd_type),                  intent(in)    :: Tbd
      type(GFS_cldprop_type),              intent(in)    :: Cldprop

#ifdef CCPP
      integer, intent(in) :: im, lm, lmk, lmp
#else
      integer, intent(out) :: im, lm, lmk, lmp
#endif
      integer,              intent(out) :: kd, kt, kb
      real(kind=kind_phys), intent(out) :: raddt

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+1+LTP), intent(out) :: plvl
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: plyr
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+1+LTP), intent(out) :: tlvl
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: tlyr
      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                  intent(out) :: tsfg
      real(kind=kind_phys), dimension(size(Grid%xlon,1)),                  intent(out) :: tsfa
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: qlyr
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: olyr

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_co2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_n2o
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_ch4
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_o2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_co
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_cfc11
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_cfc12
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_cfc22
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_ccl4
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: gasvmr_cfc113

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW), intent(out) :: faersw1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW), intent(out) :: faersw2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW), intent(out) :: faersw3
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW), intent(out) :: faerlw1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW), intent(out) :: faerlw2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW), intent(out) :: faerlw3

      real(kind=kind_phys), dimension(size(Grid%xlon,1),NSPC1),            intent(out) :: aerodp
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds2
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds3
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds4
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds5
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds6
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds7
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds8
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP),   intent(out) :: clouds9
      real(kind=kind_phys), dimension(size(Grid%xlon,1),5),                intent(out) :: cldsa
      integer,              dimension(size(Grid%xlon,1),3),                intent(out) :: mbota
      integer,              dimension(size(Grid%xlon,1),3),                intent(out) :: mtopa

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Local variables
      integer :: me, nfxr, ntrac

      integer :: i, j, k, k1, lv, itop, ibtc, LP1, lla, llb, lya, lyb

      real(kind=kind_phys) :: es, qs, delt, tem0d

      real(kind=kind_phys), dimension(size(Grid%xlon,1)) :: cvt1, cvb1, tem1d, tskn

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP) :: &
                          htswc, htlwc, gcice, grain, grime, htsw0, htlw0, &
                          rhly, tvly,qstl, vvel, clw, ciw, prslk1, tem2da, &
                          tem2db, cldcov, deltaq, cnvc, cnvw

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,2:Model%ntrac) :: tracer1
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NF_CLDS) :: clouds
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NF_VGAS) :: gasvmr
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDSW,NF_AESW)::faersw
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levr+LTP,NBDLW,NF_AELW)::faerlw
!
!===> ...  begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !--- set commonly used integers
      me = Model%me
#ifndef CCPP
      LM = Model%levr
      IM = size(Grid%xlon,1)
#endif
      NFXR = Model%nfxr
      NTRAC = Model%ntrac        ! tracers in grrad strip off sphum - start tracer1(2:NTRAC)

      LP1 = LM + 1               ! num of in/out levels


!  --- ...  set local /level/layer indexes corresponding to in/out
!  variables

#ifndef CCPP
      LMK = LM + LTP             ! num of local layers
      LMP = LMK + 1              ! num of local levels
#endif

      if ( lextop ) then
        if ( ivflip == 1 ) then    ! vertical from sfc upward
          kd = 0                   ! index diff between in/out and local
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
          lla = LMK                ! local index at the 2nd level from top
          llb = LMP                ! local index at toa level
          lya = LM                 ! local index for the 2nd layer from top
          lyb = LP1                ! local index for the top layer
        else                       ! vertical from toa downward
          kd = 1                   ! index diff between in/out and local
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
          lla = 2                  ! local index at the 2nd level from top
          llb = 1                  ! local index at toa level
          lya = 2                  ! local index for the 2nd layer from top
          lyb = 1                  ! local index for the top layer
        endif                    ! end if_ivflip_block
      else
        kd = 0
        if ( ivflip == 1 ) then  ! vertical from sfc upward
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
        else                     ! vertical from toa downward
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
        endif                    ! end if_ivflip_block
      endif   ! end if_lextop_block

      raddt = min(Model%fhswr, Model%fhlwr)
!     print *,' in grrad : raddt=',raddt


!> -# Setup surface ground temperature and ground/air skin temperature
!! if required.

      if ( itsfc == 0 ) then            ! use same sfc skin-air/ground temp
        do i = 1, IM
          tskn(i) = Sfcprop%tsfc(i)
          tsfg(i) = Sfcprop%tsfc(i)
        enddo
      else                              ! use diff sfc skin-air/ground temp
        do i = 1, IM
          tskn(i) = Sfcprop%tsfc(i)
          tsfg(i) = Sfcprop%tsfc(i)
        enddo
      endif


!> -# Prepare atmospheric profiles for radiation input.
!
!           convert pressure unit from pa to mb
      do k = 1, LM
        k1 = k + kd
        do i = 1, IM
          plvl(i,k1)   = 0.01 * Statein%prsi(i,k)   ! pa to mb (hpa)
          plyr(i,k1)   = 0.01 * Statein%prsl(i,k)   ! pa to mb (hpa)
          tlyr(i,k1)   = Statein%tgrs(i,k)
          prslk1(i,k1) = Statein%prslk(i,k)

!>  - Compute relative humidity.
!         es  = min( Statein%prsl(i,k), 0.001 * fpvs( Statein%tgrs(i,k)   ) )   ! fpvs in pa
          es  = min( Statein%prsl(i,k),  fpvs( Statein%tgrs(i,k) ) )  ! fpvs and prsl in pa
          qs  = max( QMIN, eps * es / (Statein%prsl(i,k) + epsm1*es) )
          rhly(i,k1) = max( 0.0, min( 1.0, max(QMIN, Statein%qgrs(i,k,1))/qs ) )
          qstl(i,k1) = qs
        enddo
      enddo

      !--- recast remaining all tracers (except sphum) forcing them all
      !to be positive
      do j = 2, NTRAC
        do k = 1, LM
          k1 = k + kd
          tracer1(:,k1,j) = max(0.0,Statein%qgrs(:,k,j))
        enddo
      enddo

      do i = 1, IM
        plvl(i,LP1+kd) = 0.01 * Statein%prsi(i,LP1)  ! pa to mb (hpa)
      enddo

      if ( lextop ) then                 ! values for extra top layer
        do i = 1, IM
          plvl(i,llb) = prsmin
          if ( plvl(i,lla) <= prsmin ) plvl(i,lla) = 2.0*prsmin
          plyr(i,lyb)   = 0.5 * plvl(i,lla)
          tlyr(i,lyb)   = tlyr(i,lya)
          prslk1(i,lyb) = (plyr(i,lyb)*0.00001) ** rocp ! plyr in Pa
          rhly(i,lyb)   = rhly(i,lya)
          qstl(i,lyb)   = qstl(i,lya)
        enddo

!  ---  note: may need to take care the top layer amount
       tracer1(:,lyb,:) = tracer1(:,lya,:)
      endif


!>  - Get layer ozone mass mixing ratio (if use ozone climatology data,
!!    call getozn()).

      if (Model%ntoz > 0) then            ! interactive ozone generation
        olyr(:,:) = max( QMIN, tracer1(:,1:LMK,Model%ntoz) )
      else                                ! climatological ozone
        call getozn (prslk1, Grid%xlat, IM, LMK,    &     !  ---  inputs
                     olyr)                                !  ---  outputs
      endif                               ! end_if_ntoz

!>  - Call coszmn(), to compute cosine of zenith angle.
      call coszmn (Grid%xlon,Grid%sinlat,           &     !  ---  inputs
                   Grid%coslat,Model%solhr, IM, me, &
                   Radtend%coszen, Radtend%coszdg)        !  --- outputs


!>  - Call getgases(), to set up non-prognostic gas volume mixing
!!    ratioes (gasvmr).
!  - gasvmr(:,:,1)  -  co2 volume mixing ratio
!  - gasvmr(:,:,2)  -  n2o volume mixing ratio
!  - gasvmr(:,:,3)  -  ch4 volume mixing ratio
!  - gasvmr(:,:,4)  -  o2  volume mixing ratio
!  - gasvmr(:,:,5)  -  co  volume mixing ratio
!  - gasvmr(:,:,6)  -  cf11 volume mixing ratio
!  - gasvmr(:,:,7)  -  cf12 volume mixing ratio
!  - gasvmr(:,:,8)  -  cf22 volume mixing ratio
!  - gasvmr(:,:,9)  -  ccl4 volume mixing ratio
!  - gasvmr(:,:,10) -  cfc113 volumne mixing ratio

!  --- ...  set up non-prognostic gas volume mixing ratioes

      call getgases (plvl, Grid%xlon, Grid%xlat, IM, LMK,  & !  --- inputs
                     gasvmr)                                 !  --- outputs

!CCPP: re-assign gasvmr(:,:,NF_VGAS) to gasvmr_X(:,:)
      do k = 1, LMK
        do i = 1, IM
           gasvmr_co2    (i,k)  = gasvmr(i,k,1)
           gasvmr_n2o    (i,k)  = gasvmr(i,k,2)
           gasvmr_ch4    (i,k)  = gasvmr(i,k,3)
           gasvmr_o2     (i,k)  = gasvmr(i,k,4)
           gasvmr_co     (i,k)  = gasvmr(i,k,5)
           gasvmr_cfc11  (i,k)  = gasvmr(i,k,6)
           gasvmr_cfc12  (i,k)  = gasvmr(i,k,7)   
           gasvmr_cfc22  (i,k)  = gasvmr(i,k,8)    
           gasvmr_ccl4   (i,k)  = gasvmr(i,k,9)  
           gasvmr_cfc113 (i,k)  = gasvmr(i,k,10) 
         enddo
      enddo

!>  - Get temperature at layer interface, and layer moisture.
      do k = 2, LMK
        do i = 1, IM
          tem2da(i,k) = log( plyr(i,k) )
          tem2db(i,k) = log( plvl(i,k) )
        enddo
      enddo

      if (ivflip == 0) then              ! input data from toa to sfc
        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = 1.0
          tsfa  (i)   = tlyr(i,LMK)                  ! sfc layer air temp
          tlvl(i,1)   = tlyr(i,1)
          tlvl(i,LMP) = tskn(i)
        enddo

        do k = 1, LM
          k1 = k + kd
          do i = 1, IM
            qlyr(i,k1) = max( tem1d(i), Statein%qgrs(i,k,1) )
            tem1d(i)   = min( QME5, qlyr(i,k1) )
            tvly(i,k1) = Statein%tgrs(i,k) * (1.0 + fvirt*qlyr(i,k1)) ! virtual T (K)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        do k = 2, LMK
          do i = 1, IM
            tlvl(i,k) = tlyr(i,k) + (tlyr(i,k-1) - tlyr(i,k))           &
     &                * (tem2db(i,k)   - tem2da(i,k))                   &
     &                / (tem2da(i,k-1) - tem2da(i,k))
          enddo
        enddo

      else                               ! input data from sfc to toa

        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = log( plvl(i,1) )
          tsfa  (i)   = tlyr(i,1)                    ! sfc layer air temp
          tlvl(i,1)   = tskn(i)
          tlvl(i,LMP) = tlyr(i,LMK)
        enddo

        do k = LM, 1, -1
          do i = 1, IM
            qlyr(i,k) = max( tem1d(i), Statein%qgrs(i,k,1) )
            tem1d(i)  = min( QME5, qlyr(i,k) )
            tvly(i,k) = Statein%tgrs(i,k) * (1.0 + fvirt*qlyr(i,k)) ! virtual T (K)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
          enddo
        endif

        do k = 1, LMK-1
          do i = 1, IM
            tlvl(i,k+1) = tlyr(i,k) + (tlyr(i,k+1) - tlyr(i,k))         &
     &                  * (tem2db(i,k+1) - tem2da(i,k))                 &
     &                  / (tem2da(i,k+1) - tem2da(i,k))
          enddo
        enddo

      endif                              ! end_if_ivflip

!>  - Call module_radiation_aerosols::setaer(),to setup aerosols
!! property profile for radiation.

!check  print *,' in grrad : calling setaer '

      call setaer (plvl, plyr, prslk1, tvly, rhly, Sfcprop%slmsk, & !  ---  inputs
                   tracer1, Grid%xlon, Grid%xlat, IM, LMK, LMP,   &
                   Model%lsswr, Model%lslwr,                      &
                   faersw, faerlw, aerodp)                          !  ---  outputs

! CCPP
      do j = 1,NBDSW
        do k = 1, LMK
          do i = 1, IM 
            ! NF_AESW = 3
            faersw1(i,k,j) = faersw(i,k,j,1)
            faersw2(i,k,j) = faersw(i,k,j,2)
            faersw3(i,k,j) = faersw(i,k,j,3)
          enddo 
        enddo
       enddo

      do j = 1,NBDLW
        do k = 1, LMK
          do i = 1, IM
            ! NF_AELW = 3
            faerlw1(i,k,j) = faerlw(i,k,j,1)
            faerlw2(i,k,j) = faerlw(i,k,j,2)
            faerlw3(i,k,j) = faerlw(i,k,j,3)
          enddo
        enddo
       enddo

!>  - Obtain cloud information for radiation calculations
!!    (clouds,cldsa,mtopa,mbota)
!!\n   for  prognostic cloud:
!!    - For Zhao/Moorthi's prognostic cloud scheme,
!!      call module_radiation_clouds::progcld1()
!!    - For Zhao/Moorthi's prognostic cloud+pdfcld,
!!      call module_radiation_clouds::progcld3()
!!      call module_radiation_clouds::progclduni() for unified cloud and ncld=2
!>  - If cloud condensate is not computed (ntcw=0), using the legacy
!!   cloud scheme, compute cloud information based on Slingo's
!!   diagnostic cloud scheme (call module_radiation_clouds::diagcld1())

!  --- ...  obtain cloud information for radiation calculations

      if (Model%ntcw > 0) then                   ! prognostic cloud scheme
        if (Model%uni_cld .and. Model%ncld >= 2) then
          clw(:,:) = tracer1(:,1:LMK,Model%ntcw)              ! cloud water amount
          ciw(:,:) = 0.0
          do j = 2, Model%ncld
            ciw(:,:) = ciw(:,:) + tracer1(:,1:LMK,Model%ntcw+j-1)   ! cloud ice amount
          enddo

          do k = 1, LMK
            do i = 1, IM
              if ( clw(i,k) < EPSQ ) clw(i,k) = 0.0
              if ( ciw(i,k) < EPSQ ) ciw(i,k) = 0.0
            enddo
          enddo
        else
          clw(:,:) = 0.0
          do j = 1, Model%ncld
            clw(:,:) = clw(:,:) + tracer1(:,1:LMK,Model%ntcw+j-1)   ! cloud condensate amount
          enddo

          do k = 1, LMK
            do i = 1, IM
              if ( clw(i,k) < EPSQ ) clw(i,k) = 0.0
            enddo
          enddo
        endif
!
!  --- add suspended convective cloud water to grid-scale cloud water
!      only for cloud fraction & radiation computation
!      it is to enhance cloudiness due to suspended convec cloud water
!      for zhao/moorthi's (icmphys=1) &
!          ferrier's (icmphys=2) microphysics schemes
!
        if (Model%shoc_cld) then                                       ! all but MG microphys
          cldcov(:,1:LM) = Tbd%phy_f3d(:,1:LM,Model%ntot3d-2)
        elseif (Model%ncld == 2) then                                  ! MG microphys (icmphys = 1)
          cldcov(:,1:LM) = Tbd%phy_f3d(:,1:LM,1)
        else                                                           ! neither of the other two cases
          cldcov = 0
        endif

        if ((Model%num_p3d == 4) .and. (Model%npdf3d == 3)) then       ! icmphys = 3
          deltaq(:,1:LM) = Tbd%phy_f3d(:,1:LM,5)
          cnvw  (:,1:LM) = Tbd%phy_f3d(:,1:LM,6)
          cnvc  (:,1:LM) = Tbd%phy_f3d(:,1:LM,7)
        elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then  ! icmphys = 1
          deltaq(:,1:LM) = 0.
          cnvw  (:,1:LM) = Tbd%phy_f3d(:,1:LM,Model%num_p3d+1)
          cnvc  (:,1:LM) = 0.
        else                                                           !  icmphys = 1 (ncld=2)
          deltaq = 0.0
          cnvw   = 0.0
          cnvc   = 0.0
        endif

        if (lextop) then
          cldcov(:,lyb) = cldcov(:,lya)
          deltaq(:,lyb) = deltaq(:,lya)
          cnvw  (:,lyb) = cnvw  (:,lya)
          cnvc  (:,lyb) = cnvc  (:,lya)
        endif

        if (icmphys == 1) then
          clw(:,1:LMK) = clw(:,1:LMK) + cnvw(:,1:LMK)
        endif
!

        if (icmphys == 1) then           ! zhao/moorthi's prognostic cloud scheme
                                         ! or unified cloud and/or with MG microphysics

          if (Model%uni_cld .and. Model%ncld >= 2) then
            call progclduni (plyr, plvl, tlyr, tvly, clw, ciw,    &    !  ---  inputs
                             Grid%xlat, Grid%xlon, Sfcprop%slmsk, &
                             IM, LMK, LMP, cldcov(:,1:LMK),       &
                             clouds, cldsa, mtopa, mbota)              !  ---  outputs
          else
            call progcld1 (plyr ,plvl, tlyr, tvly, qlyr, qstl,    &    !  ---  inputs
                           rhly, clw, Grid%xlat,Grid%xlon,        &
                           Sfcprop%slmsk, IM, LMK, LMP,           &
                           Model%uni_cld, Model%lmfshal,          &
                           Model%lmfdeep2, cldcov(:,1:LMK),       &
                           clouds, cldsa, mtopa, mbota)                !  ---  outputs
          endif

        elseif(icmphys == 3) then      ! zhao/moorthi's prognostic cloud+pdfcld

          call progcld3 (plyr, plvl, tlyr, tvly, qlyr, qstl, rhly,&    ! ---  inputs
                         clw, cnvw, cnvc, Grid%xlat, Grid%xlon,   &
                         Sfcprop%slmsk,im, lmk, lmp, deltaq,      &
                         Model%sup, Model%kdt, me,                &
                         clouds, cldsa, mtopa, mbota)                  ! ---  outputs

        endif                            ! end if_icmphys

      else                               ! diagnostic cloud scheme

        cvt1(:) = 0.01 * Cldprop%cvt(:)
        cvb1(:) = 0.01 * Cldprop%cvb(:)

        do k = 1, LM
          k1 = k + kd
          vvel(:,k1) = 0.01 * Statein%vvl(:,k)
        enddo
        if (lextop) then
          vvel(:,lyb) = vvel(:,lya)
        endif

!  ---  compute diagnostic cloud related quantities

        call diagcld1 (plyr, plvl, tlyr, rhly, vvel, Cldprop%cv,  &    !  ---  inputs
                       cvt1, cvb1, Grid%xlat, Grid%xlon,          &
                       Sfcprop%slmsk, IM, LMK, LMP,               &
                       clouds, cldsa, mtopa, mbota)                    !  ---  outputs

      endif                                ! end_if_ntcw

! CCPP
       do k = 1, LMK
         do i = 1, IM
            clouds1(i,k) = clouds(i,k,1)
            clouds2(i,k) = clouds(i,k,2)
            clouds3(i,k) = clouds(i,k,3)
            clouds4(i,k) = clouds(i,k,4)
            clouds5(i,k) = clouds(i,k,5)
            clouds6(i,k) = clouds(i,k,6)
            clouds7(i,k) = clouds(i,k,7)
            clouds8(i,k) = clouds(i,k,8)
            clouds9(i,k) = clouds(i,k,9)
         enddo
       enddo

      end subroutine GFS_rrtmg_pre_run
   
!> \section arg_table_GFS_rrtmg_pre_finalize Argument Table
!!
      subroutine GFS_rrtmg_pre_finalize ()
      end subroutine GFS_rrtmg_pre_finalize

!! @}
      end module GFS_rrtmg_pre


