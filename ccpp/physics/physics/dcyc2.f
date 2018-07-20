      module dcyc2t3_pre

      contains

!! \section arg_table_dcyc2t3_pre_init Argument Table
!!
      subroutine dcyc2t3_pre_init()
      end subroutine dcyc2t3_pre_init


!! \section arg_table_dcyc2t3_pre_run Argument Table
!!
      subroutine dcyc2t3_pre_run()
      end subroutine dcyc2t3_pre_run


!! \section arg_table_dcyc2t3_pre_finalize Argument Table
!!
      subroutine dcyc2t3_pre_finalize()
      end subroutine dcyc2t3_pre_finalize

      end module dcyc2t3_pre


      module dcyc2t3

      contains

!! \section arg_table_dcyc2t3_init Argument Table
!!
      subroutine dcyc2t3_init()
      end subroutine dcyc2t3_init


! ===================================================================== !
!  description:                                                         !
!                                                                       !
!    dcyc2t3 fits radiative fluxes and heating rates from a coarse      !
!    radiation calc time interval into model's more frequent time steps.!
!    solar heating rates and fluxes are scaled by the ratio of cosine   !
!    of zenith angle at the current time to the mean value used in      !
!    radiation calc.  surface downward lw flux is scaled by the ratio   !
!    of current surface air temperature (temp**4) to the corresponding  !
!    temperature saved during lw radiation calculation. upward lw flux  !
!    at the surface is computed by current ground surface temperature.  !
!    surface emissivity effect will be taken in other part of the model.!
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call dcyc2t3                                                       !
!      inputs:                                                          !
!          ( solhr,slag,sdec,cdec,sinlat,coslat,                        !
!            xlon,coszen,tsea,tf,tsflw,sfcemis,                         !
!            sfcdsw,sfcnsw,sfcdlw,swh,swhc,hlw,hlwc,                    !
!            sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   !
!            sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   !
!            ix, im, levs,                                              !
!      input/output:                                                    !
!            dtdt,dtdtc,                                                !
!      outputs:                                                         !
!            adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,xmu,xcosz,         !
!            adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   !
!            adjdnnbmd,adjdnndfd,adjdnvbmd,adjdnvdfd)                   !
!                                                                       !
!                                                                       !
!  program history:                                                     !
!          198?  nmc mrf    - created, similar as treatment in gfdl     !
!                             radiation treatment                       !
!          1994  y. hou     - modified solar zenith angle calculation   !
!     nov  2004  x. wu      - add sfc sw downward flux to the variable  !
!                             list for sea-ice model                    !
!     mar  2008  y. hou     - add cosine of zenith angle as output for  !
!                             sunshine duration time calc.              !
!     sep  2008  y. hou     - separate net sw and downward lw in slrad, !
!                 changed the sign of sfc net sw to consistent with     !
!                 other parts of the mdl (positive value defines from   !
!                 atmos to the ground). rename output fluxes as adjusted!
!                 fluxes. other minor changes such as renaming some of  !
!                 passing argument names to be consistent with calling  !
!                 program.                                              !
!     apr  2009  y. hou     - integrated with the new parallel model    !
!                 along with other modifications                        !
!     mar  2011  y. hou     - minor modification including rearrange    !
!                 loop orders and loop structures to improve efficiency !
!     mar  2014  x. wu      - add sfc nir/vis bm/df to the variable     !
!                             list for the coupled model input          !
!     jul  2014  s moorthi  - merge gfs and nems versions               !
!     jun  2014  y. hou     - revised to include both up and down sw    !
!                 spectral component fluxes                             !
!     Oct  2014  y. hous s. moorthi - add emissivity contribution to    !
!                             upward longwave flux                      !
!                                                                       !
!  subprograms called:  none                                            !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                              !
!     solhr        - real, forecast time in 24-hour form (hr)           !
!     slag         - real, equation of time in radians                  !
!     sdec, cdec   - real, sin and cos of the solar declination angle   !
!     sinlat(im), coslat(im):                                           !
!                  - real, sin and cos of latitude                      !
!     xlon   (im)  - real, longitude in radians                         !
!     coszen (im)  - real, avg of cosz over daytime sw call interval    !
!     tsea   (im)  - real, ground surface temperature (k)               !
!     tf     (im)  - real, surface air (layer 1) temperature (k)        !
!     sfcemis(im)  - real, surface emissivity (fraction)                !
!     tsflw  (im)  - real, sfc air (layer 1) temp in k saved in lw call !
!     sfcdsw (im)  - real, total sky sfc downward sw flux ( w/m**2 )    !
!     sfcnsw (im)  - real, total sky sfc net sw into ground (w/m**2)    !
!     sfcdlw (im)  - real, total sky sfc downward lw flux ( w/m**2 )    !
!     swh(ix,levs) - real, total sky sw heating rates ( k/s )           !
!     swhc(ix,levs) - real, clear sky sw heating rates ( k/s )          !
!     hlw(ix,levs) - real, total sky lw heating rates ( k/s )           !
!     hlwc(ix,levs) - real, clear sky lw heating rates ( k/s )          !
!     sfcnirbmu(im)- real, tot sky sfc nir-beam sw upward flux (w/m2)   !
!     sfcnirdfu(im)- real, tot sky sfc nir-diff sw upward flux (w/m2)   !
!     sfcvisbmu(im)- real, tot sky sfc uv+vis-beam sw upward flux (w/m2)!
!     sfcvisdfu(im)- real, tot sky sfc uv+vis-diff sw upward flux (w/m2)!
!     sfcnirbmd(im)- real, tot sky sfc nir-beam sw downward flux (w/m2) !
!     sfcnirdfd(im)- real, tot sky sfc nir-diff sw downward flux (w/m2) !
!     sfcvisbmd(im)- real, tot sky sfc uv+vis-beam sw dnward flux (w/m2)!
!     sfcvisdfd(im)- real, tot sky sfc uv+vis-diff sw dnward flux (w/m2)!
!     ix, im       - integer, horiz. dimention and num of used points   !
!     levs         - integer, vertical layer dimension                  !
!                                                                       !
!  input/output:                                                        !
!     dtdt(im,levs)- real, model time step adjusted total radiation     !
!                          heating rates ( k/s )                        !
!     dtdtc(im,levs)- real, model time step adjusted clear sky radiation!
!                          heating rates ( k/s )                        !
!                                                                       !
!  outputs:                                                             !
!     adjsfcdsw(im)- real, time step adjusted sfc dn sw flux (w/m**2)   !
!     adjsfcnsw(im)- real, time step adj sfc net sw into ground (w/m**2)!
!     adjsfcdlw(im)- real, time step adjusted sfc dn lw flux (w/m**2)   !
!     adjsfculw(im)- real, sfc upward lw flux at current time (w/m**2)  !
!     adjnirbmu(im)- real, t adj sfc nir-beam sw upward flux (w/m2)     !
!     adjnirdfu(im)- real, t adj sfc nir-diff sw upward flux (w/m2)     !
!     adjvisbmu(im)- real, t adj sfc uv+vis-beam sw upward flux (w/m2)  !
!     adjvisdfu(im)- real, t adj sfc uv+vis-diff sw upward flux (w/m2)  !
!     adjnirbmd(im)- real, t adj sfc nir-beam sw downward flux (w/m2)   !
!     adjnirdfd(im)- real, t adj sfc nir-diff sw downward flux (w/m2)   !
!     adjvisbmd(im)- real, t adj sfc uv+vis-beam sw dnward flux (w/m2)  !
!     adjvisdfd(im)- real, t adj sfc uv+vis-diff sw dnward flux (w/m2)  !
!     xmu   (im)   - real, time step zenith angle adjust factor for sw  !
!     xcosz (im)   - real, cosine of zenith angle at current time step  !
!                                                                       !
!  ====================    end of description    =====================  !

!-----------------------------------
!! \section arg_table_dcyc2t3_run Argument Table
!! | local_name     | standard_name                                                                                  | long_name                                                                                            | units   | rank | type      | kind      | intent | optional |
!! |----------------|------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | solhr          | forecast_hour                                                                                  | forecast time in 24-hour form                                                                        | h       |    0 | real      | kind_phys | in     | F        |
!! | slag           | equation_of_time                                                                               | equation of time                                                                                     | radians |    0 | real      | kind_phys | in     | F        |
!! | sdec           | sine_of_solar_declination_angle                                                                | sine of solar declination angle                                                                      | none    |    0 | real      | kind_phys | in     | F        |
!! | cdec           | cosine_of_solar_declination_angle                                                              | cosine of solar declination angle                                                                    | none    |    0 | real      | kind_phys | in     | F        |
!! | sinlat         | sine_of_latitude                                                                               | sine of latitude                                                                                     | none    |    1 | real      | kind_phys | in     | F        |
!! | coslat         | cosine_of_latitude                                                                             | cosine of latitude                                                                                   | none    |    1 | real      | kind_phys | in     | F        |
!! | xlon           | longitude                                                                                      | longitude of grid box                                                                                | radians |    1 | real      | kind_phys | in     | F        |
!! | coszen         | cosine_of_zenith_angle                                                                         | average of cosine of zenith angle over daytime shortwave call time interval                          | none    |    1 | real      | kind_phys | in     | F        |
!! | tsea           | surface_skin_temperature                                                                       | surface skin temperature                                                                             | K       |    1 | real      | kind_phys | in     | F        |
!! | tf             | air_temperature_at_lowest_model_layer                                                          | air temperature at lowest model layer                                                                | K       |    1 | real      | kind_phys | in     | F        |
!! | tsflw          | surface_midlayer_air_temperature_in_longwave_radiation                                         | surface (first layer) air temperature saved in longwave radiation call                               | K       |    1 | real      | kind_phys | in     | F        |
!! | sfcemis        | surface_longwave_emissivity                                                                    | surface emissivity                                                                                   | frac    |    1 | real      | kind_phys | in     | F        |
!! | sfcdsw         | surface_downwelling_shortwave_flux_on_radiation_time_step                                      | total sky surface downwelling shortwave flux on radiation time step                                  | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcnsw         | surface_net_downwelling_shortwave_flux_on_radiation_time_step                                  | total sky surface net downwelling shortwave flux on radiation time step                              | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcdlw         | surface_downwelling_longwave_flux_on_radiation_time_step                                       | total sky surface downwelling longwave flux on radiation time step                                   | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | swh            | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | total sky shortwave heating rate on radiation time step                                              | K s-1   |    2 | real      | kind_phys | in     | F        |
!! | swhc           | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | clear sky shortwave heating rate on radiation time step                                              | K s-1   |    2 | real      | kind_phys | in     | F        |
!! | hlw            | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                     | total sky longwave heating rate on radiation time step                                               | K s-1   |    2 | real      | kind_phys | in     | F        |
!! | hlwc           | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step  | clear sky longwave heating rate on radiation time step                                               | K s-1   |    2 | real      | kind_phys | in     | F        |
!! | sfcnirbmu      | surface_upwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step                   | total sky surface upwelling beam near-infrared shortwave flux on radiation time step                 | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcnirdfu      | surface_upwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step                  | total sky surface upwelling diffuse near-infrared shortwave flux on radiation time step              | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcvisbmu      | surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step         | total sky surface upwelling beam ultraviolet plus visible shortwave flux on radiation time step      | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcvisdfu      | surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step        | total sky surface upwelling diffuse ultraviolet plus visible shortwave flux on radiation time step   | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcnirbmd      | surface_downwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step                 | total sky surface downwelling beam near-infrared shortwave flux on radiation time step               | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcnirdfd      | surface_downwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step                | total sky surface downwelling diffuse near-infrared shortwave flux on radiation time step            | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcvisbmd      | surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step       | total sky surface downwelling beam ultraviolet plus visible shortwave flux on radiation time step    | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | sfcvisdfd      | surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step      | total sky surface downwelling diffuse ultraviolet plus visible shortwave flux on radiation time step | W m-2   |    1 | real      | kind_phys | in     | F        |
!! | ix             | horizontal_dimension                                                                           | horizontal dimension                                                                                 | count   |    0 | integer   |           | in     | F        |
!! | im             | horizontal_loop_extent                                                                         | horizontal loop extent                                                                               | count   |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                                             | number of vertical layers                                                                            | count   |    0 | integer   |           | in     | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics                                               | total radiative heating rate at current time                                                         | K s-1   |    2 | real      | kind_phys | inout  | F        |
!! | dtdtc          | tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky                        | clear sky radiative (shortwave + longwave) heating rate at current time                              | K s-1   |    2 | real      | kind_phys | inout  | F        |
!! | adjsfcdsw      | surface_downwelling_shortwave_flux                                                             | surface downwelling shortwave flux at current time                                                   | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjsfcnsw      | surface_net_downwelling_shortwave_flux                                                         | surface net downwelling shortwave flux at current time                                               | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux                                                              | surface downwelling longwave flux at current time                                                    | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjsfculw      | surface_upwelling_longwave_flux                                                                | surface upwelling longwave flux at current time                                                      | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | xmu            | zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes                                   | zenith angle temporal adjustment factor for shortwave fluxes                                         | none    |    1 | real      | kind_phys | out    | F        |
!! | xcosz          | instantaneous_cosine_of_zenith_angle                                                           | cosine of zenith angle at current time                                                               | none    |    1 | real      | kind_phys | out    | F        |
!! | adjnirbmu      | surface_upwelling_direct_near_infrared_shortwave_flux                                          | surface upwelling beam near-infrared shortwave flux at current time                                  | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjnirdfu      | surface_upwelling_diffuse_near_infrared_shortwave_flux                                         | surface upwelling diffuse near-infrared shortwave flux at current time                               | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjvisbmu      | surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux                                | surface upwelling beam ultraviolet plus visible shortwave flux at current time                       | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjvisdfu      | surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux                               | surface upwelling diffuse ultraviolet plus visible shortwave flux at current time                    | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjnirbmd      | surface_downwelling_direct_near_infrared_shortwave_flux                                        | surface downwelling beam near-infrared shortwave flux at current time                                | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjnirdfd      | surface_downwelling_diffuse_near_infrared_shortwave_flux                                       | surface downwelling diffuse near-infrared shortwave flux at current time                             | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjvisbmd      | surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux                              | surface downwelling beam ultraviolet plus visible shortwave flux at current time                     | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | adjvisdfd      | surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux                             | surface downwelling diffuse ultraviolet plus visible shortwave flux at current time                  | W m-2   |    1 | real      | kind_phys | out    | F        |
!! | errmsg         | error_message                                                                                  | error message for error handling in CCPP                                                             | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                                                     | error flag for error handling in CCPP                                                                | flag    |    0 | integer   |           | out    | F        |
!!
      subroutine dcyc2t3_run                                            &
!...................................
!  ---  inputs:
     &     ( solhr,slag,sdec,cdec,sinlat,coslat,                        &
     &       xlon,coszen,tsea,tf,tsflw,sfcemis,                         &
     &       sfcdsw,sfcnsw,sfcdlw,swh,swhc,hlw,hlwc,                    &
     &       sfcnirbmu,sfcnirdfu,sfcvisbmu,sfcvisdfu,                   &
     &       sfcnirbmd,sfcnirdfd,sfcvisbmd,sfcvisdfd,                   &
     &       ix, im, levs,                                              &
!  ---  input/output:
     &       dtdt,dtdtc,                                                &
!  ---  outputs:
     &       adjsfcdsw,adjsfcnsw,adjsfcdlw,adjsfculw,xmu,xcosz,         &
     &       adjnirbmu,adjnirdfu,adjvisbmu,adjvisdfu,                   &
     &       adjnirbmd,adjnirdfd,adjvisbmd,adjvisdfd,                   &
     &       errmsg,errflg
     &     )
!
      use machine,         only : kind_phys
      use physcons,        only : con_pi, con_sbc

      implicit none
!
!  ---  constant parameters:
      real(kind=kind_phys), parameter :: f_eps  = 0.0001, hour12 = 12.0

!  ---  inputs:
      integer, intent(in) :: ix, im, levs

      real(kind=kind_phys), intent(in) :: solhr, slag, cdec, sdec

      real(kind=kind_phys), dimension(im), intent(in) ::                &
     &      sinlat, coslat, xlon, coszen, tsea, tf, tsflw, sfcdlw,      &
     &      sfcdsw, sfcnsw, sfcemis
      real(kind=kind_phys), dimension(im), intent(in) ::                &
     &      sfcnirbmu, sfcnirdfu, sfcvisbmu, sfcvisdfu,                 &
     &      sfcnirbmd, sfcnirdfd, sfcvisbmd, sfcvisdfd

      real(kind=kind_phys), dimension(ix,levs), intent(in) :: swh,  hlw
     &,                                                       swhc, hlwc&

!  ---  input/output:
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: dtdt   &
     &,                                                          dtdtc

!  ---  outputs:
      real(kind=kind_phys), dimension(im), intent(out) ::               &
     &      adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw, xmu, xcosz,     &
     &      adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                 &
     &      adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      integer :: i, k
      real(kind=kind_phys) :: cns, ss, cc, ch, tem1, tem2
!
!===> ...  begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
      cns = con_pi * (solhr - hour12) / hour12 + slag
!
      do i = 1, im

!  --- ...  lw time-step adjustment
!           -----------------------
!  --- ...  adjust sfc downward lw flux to account for t changes in layer 1
!           compute 4th power of the ratio of layer 1 tf over the mean value tsflw

        tem1 = tf(i) / tsflw(i)
        tem2  = tem1 * tem1
        adjsfcdlw(i) = sfcdlw(i) * tem2 * tem2

!  --- ...  compute sfc upward lw flux from current sfc temp,
!      note: sfc emiss effect is not appied here, and will be dealt in other place

        tem2 = tsea(i) * tsea(i)
        adjsfculw(i) =  sfcemis(i) * con_sbc * tem2 * tem2
     &               + (1.0 - sfcemis(i)) * adjsfcdlw(i)
!
!  --- ...  sw time-step adjustment
!           -----------------------

        ss     = sinlat(i) * sdec
        cc     = coslat(i) * cdec
        ch     = cc * cos( xlon(i)+cns )
        xcosz(i) = ch + ss        ! cosine of solar zenith angle at current time

!  --- ...  normalize by average value over radiation period for daytime.

        if ( xcosz(i) > f_eps .and. coszen(i) > f_eps ) then
          xmu(i) = xcosz(i) / coszen(i)
        else
          xmu(i) = 0.0
        endif

!  --- ...  adjust sfc net and downward sw fluxes for zenith angle changes
!      note: sfc emiss effect will not be appied here

        adjsfcnsw(i) = sfcnsw(i)   * xmu(i)
        adjsfcdsw(i) = sfcdsw(i)   * xmu(i)

        adjnirbmu(i)  = sfcnirbmu(i) * xmu(i)
        adjnirdfu(i)  = sfcnirdfu(i) * xmu(i)
        adjvisbmu(i)  = sfcvisbmu(i) * xmu(i)
        adjvisdfu(i)  = sfcvisdfu(i) * xmu(i)

        adjnirbmd(i)  = sfcnirbmd(i) * xmu(i)
        adjnirdfd(i)  = sfcnirdfd(i) * xmu(i)
        adjvisbmd(i)  = sfcvisbmd(i) * xmu(i)
        adjvisdfd(i)  = sfcvisdfd(i) * xmu(i)
      enddo

!  --- ...  adjust sw heating rates with zenith angle change and
!           add with lw heating to temperature tendency

      do k = 1, levs
        do i = 1, im
          dtdt(i,k)  = dtdt(i,k)  + swh(i,k)*xmu(i)  + hlw(i,k)
          dtdtc(i,k) = dtdtc(i,k) + swhc(i,k)*xmu(i) + hlwc(i,k)
        enddo
      enddo
!
      return
!...................................
      end subroutine dcyc2t3_run
!-----------------------------------


!! \section arg_table_dcyc2t3_finalize Argument Table
!!
      subroutine dcyc2t3_finalize()
      end subroutine dcyc2t3_finalize

      end module dcyc2t3





      module dcyc2t3_post

      contains

!! \section arg_table_dcyc2t3_post_init Argument Table
!!
      subroutine dcyc2t3_post_init()
      end subroutine dcyc2t3_post_init


!! \section arg_table_dcyc2t3_post_run Argument Table
!! | local_name     | standard_name                          | long_name                                              | units   | rank | type                  | kind      | intent | optional |
!! |----------------|----------------------------------------|--------------------------------------------------------|---------|------|-----------------------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                 | horizontal loop extent                                 | count   |    0 | integer               |           | in     | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux      | surface downwelling longwave flux at current time      | W m-2   |    1 | real                  | kind_phys | in     | F        |
!! | adjsfculw      | surface_upwelling_longwave_flux        | surface upwelling longwave flux at current time        | W m-2   |    1 | real                  | kind_phys | in     | F        |
!! | adjsfcdsw      | surface_downwelling_shortwave_flux     | surface downwelling shortwave flux at current time     | W m-2   |    1 | real                  | kind_phys | in     | F        |
!! | adjsfcnsw      | surface_net_downwelling_shortwave_flux | surface net downwelling shortwave flux at current time | W m-2   |    1 | real                  | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                      | GFS diagnostics derived data type variable             | DDT     |    0 | GFS_diag_type         |           | inout  | F        |
!! | errmsg         | error_message                          | error message for error handling in CCPP               | none    |    0 | character             | len=*     | out    | F        |
!! | errflg         | error_flag                             | error flag for error handling in CCPP                  | flag    |    0 | integer               |           | out    | F        |
!!
      subroutine dcyc2t3_post_run(                                      &
     &           im, adjsfcdlw, adjsfculw, adjsfcdsw, adjsfcnsw, Diag,  &
     &           errmsg, errflg)

      use GFS_typedefs, only: GFS_diag_type
      use machine,      only: kind_phys

      implicit none

      integer, intent(in) :: im
      real(kind=kind_phys), dimension(im), intent(in) ::                &
     &      adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw
      type(GFS_diag_type), intent(inout) :: Diag
      character(len=*),      intent(out) :: errmsg
      integer,               intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      Diag%dlwsfci(:) = adjsfcdlw(:)
      Diag%ulwsfci(:) = adjsfculw(:)
      Diag%uswsfci(:) = adjsfcdsw(:) - adjsfcnsw(:)
      Diag%dswsfci(:) = adjsfcdsw(:)

      return

      end subroutine dcyc2t3_post_run


!! \section arg_table_dcyc2t3_post_finalize Argument Table
!!
      subroutine dcyc2t3_post_finalize()
      end subroutine dcyc2t3_post_finalize

      end module dcyc2t3_post



