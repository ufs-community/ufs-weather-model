!>  \file sfc_diag.f
!!  This file contains the land surface diagnose calculation scheme.

!> \defgroup Sfc_diag Land Surface Diagnose Calculation
!! @{

      module sfc_diag
      contains
  
      subroutine sfc_diag_init
      end subroutine sfc_diag_init
      
      subroutine sfc_diag_finalize
      end subroutine sfc_diag_finalize
      
!> \brief Brief description of the subroutine
!!
!! \section arg_table_sfc_diag_run Arguments
!! | local_name     | standard_name                                                                | long_name                                                   | units         | rank | type      |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|-------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                      | count         |    0 | integer   |           | in     | F        |
!! | ps             | surface_air_pressure                                                         | surface pressure                                            | Pa            |    1 | real      | kind_phys | in     | F        |
!! | u1             | x_wind_at_lowest_model_layer                                                 | x component of 1st model layer wind                         | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1             | y_wind_at_lowest_model_layer                                                 | y component of 1st model layer wind                         | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature_at_lowest_model_layer                                        | 1st model layer air temperature                             | K             |    1 | real      | kind_phys | in     | F        |
!! | q1             | specific_humidity_at_lowest_model_layer                                      | 1st model layer specific humidity                           | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | tskin          | surface_skin_temperature                                                     | surface skin temperature                                    | K             |    1 | real      | kind_phys | in     | F        |
!! | qsurf          | surface_specific_humidity                                                    | surface specific humidity                                   | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | f10m           | ratio_of_wind_at_lowest_model_layer_and_wind_at_10m                          | ratio of fm10 and fm                                        | ratio         |    1 | real      | kind_phys | out    | F        |
!! | u10m           | x_wind_at_10m                                                                | x component of wind at 10 m                                 | m s-1         |    1 | real      | kind_phys | out    | F        |
!! | v10m           | y_wind_at_10m                                                                | y component of wind at 10 m                                 | m s-1         |    1 | real      | kind_phys | out    | F        |
!! | t2m            | temperature_at_2m                                                            | temperature at 2 m                                          | K             |    1 | real      | kind_phys | out    | F        |
!! | q2m            | specific_humidity_at_2m                                                      | specific humidity at 2 m                                    | kg kg-1       |    1 | real      | kind_phys | out    | F        |
!! | prslki         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer | ratio         |    1 | real      | kind_phys | in     | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux                                    | surface upward evaporation flux                             | kg kg-1 m s-1 |    1 | real      | kind_phys | in     | F        |
!! | fm             | Monin-Obukhov_similarity_function_for_momentum                               | Monin-Obukhov similarity parameter for momentum             | none          |    1 | real      | kind_phys | in     | F        |
!! | fh             | Monin-Obukhov_similarity_function_for_heat                                   | Monin-Obukhov similarity parameter for heat                 | none          |    1 | real      | kind_phys | in     | F        |
!! | fm10           | Monin-Obukhov_similarity_function_for_momentum_at_10m                        | Monin-Obukhov similarity parameter for momentum             | none          |    1 | real      | kind_phys | in     | F        |
!! | fh2            | Monin-Obukhov_similarity_function_for_heat_at_2m                             | Monin-Obukhov similarity parameter for heat                 | none          |    1 | real      | kind_phys | in     | F        |
!! | errmsg         | error_message                                                                | error message for error handling in CCPP                    | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                                   | error flag for error handling in CCPP                       | flag          |    0 | integer   |           | out    | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
! DH* TODO add intent() information for variables
      subroutine sfc_diag_run                                           &
     &                   (im,ps,u1,v1,t1,q1,                            &
     &                    tskin,qsurf,f10m,u10m,v10m,t2m,q2m,           &
     &                    prslki,evap,fm,fh,fm10,fh2,errmsg,errflg      &
     &                   )
!!
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, grav => con_g,  cp => con_cp,                       &
     &              eps => con_eps, epsm1 => con_epsm1
      implicit none
!
      integer              im
      real, dimension(im) :: ps,   u1,   v1,   t1,  q1,  tskin,  qsurf, &
     &                       f10m, u10m, v10m, t2m, q2m, prslki, evap,  &
     &                       fm,   fh,   fm10, fh2
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!     locals
!
      real (kind=kind_phys), parameter :: qmin=1.0e-8
      integer              k,i
!
      real(kind=kind_phys)        fhi, qss, wrk
!     real(kind=kind_phys) sig2k, fhi, qss
!
!     real, parameter :: g=grav
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!     estimate sigma ** k at 2 m
!
!     sig2k = 1. - 4. * g * 2. / (cp * 280.)
!
!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals
!
!!
      do i = 1, im
        f10m(i) = fm10(i) / fm(i)
!       f10m(i) = min(f10m(i),1.)
        u10m(i) = f10m(i) * u1(i)
        v10m(i) = f10m(i) * v1(i)
        fhi     = fh2(i) / fh(i)
!       t2m(i)  = tskin(i)*(1. - fhi) + t1(i) * prslki(i) * fhi
!       sig2k   = 1. - (grav+grav) / (cp * t2m(i))
!       t2m(i)  = t2m(i) * sig2k
        wrk     = 1.0 - fhi

        t2m(i)  = tskin(i)*wrk + t1(i)*prslki(i)*fhi - (grav+grav)/cp

        if(evap(i) >= 0.) then !  for evaporation>0, use inferred qsurf to deduce q2m
          q2m(i) = qsurf(i)*wrk + max(qmin,q1(i))*fhi
        else                   !  for dew formation, use saturated q at tskin
          qss    = fpvs(tskin(i))
          qss    = eps * qss / (ps(i) + epsm1 * qss)
          q2m(i) = qss*wrk + max(qmin,q1(i))*fhi
        endif
        qss    = fpvs(t2m(i))
        qss    = eps * qss / (ps(i) + epsm1 * qss)
        q2m(i) = min(q2m(i),qss)
      enddo

      return
      end subroutine sfc_diag_run
!> @}

      end module sfc_diag
!> @}
