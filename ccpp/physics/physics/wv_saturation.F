!
! Common block and statement functions for saturation vapor pressure
! look-up procedure, J. J. Hack, February 1990
!
! $Id: wv_saturation.F90,v 1.1.12.2.10.1 2014-04-14 16:04:56 dbarahon Exp $
!
      module wv_saturation
#ifdef GEOS5
       use MAPL_ConstantsMod, r8 => MAPL_R8
#endif
#ifdef NEMS_GSM
       use funcphys,  only: fpvsl, fpvsi, fpvs    ! saturation vapor pressure for water & ice
#endif
       use machine,   only: r8 => kind_phys


!++jtb (comm out)



!--jtb

       implicit none
       private
       save
!
! Public interfaces
!
       public gestbl
       public estblf
       public aqsat
       public aqsatd
       public vqsatd
       public fqsatd
       public qsat_water
       public vqsat_water
       public qsat_ice
       public vqsat_ice
       public vqsatd_water
       public aqsat_water
       public vqsatd2_water
       public vqsatd2_water_single
       public vqsatd2_ice_single
       public vqsatd2
       public vqsatd2_single
       public polysvp
!
! Data used by cldwat
!
       public hlatv, tmin, hlatf, rgasv, pcf, cp, epsqs, ttrice
!
! Data
!
       integer plenest
       parameter (plenest=250)
!
! Table of saturation vapor pressure values es from tmin degrees
! to tmax+1 degrees k in one degree increments.  ttrice defines the
! transition region where es is a combination of ice & water values
!
       real(r8) estbl(plenest)
       real(r8) tmin
       real(r8) tmax
       real(r8) ttrice
       real(r8) pcf(6)
       real(r8) epsqs
       real(r8) rgasv
       real(r8) hlatf
       real(r8) hlatv
       real(r8) cp
       real(r8) tmelt
       logical icephs

       integer, parameter :: iulog=6

       contains

       real(r8) function estblf( td )
!
! Saturation vapor pressure table lookup
!
       real(r8), intent(in) :: td
!
       real(r8) :: e
       real(r8) :: ai
       integer :: i
!
       e = max(min(td,tmax),tmin)
       i = int(e-tmin)+1
       ai = aint(e-tmin)
       estblf = (tmin+ai-e+1._r8)* estbl(i)-(tmin+ai-e)* estbl(i+1)
       end function estblf

      subroutine gestbl(tmn ,tmx ,trice ,ip ,epsil , latvap ,latice ,   &
     &                  rh2o ,cpair ,tmeltx )
!-----------------------------------------------------------------------
!
! Purpose:
! Builds saturation vapor pressure table for later lookup procedure.
!
! Method:
! Uses Goff & Gratch (1946) relationships to generate the table
! according to a set of free parameters defined below.  Auxiliary
! routines are also included for making rapid estimates (well with 1%)
! of both es and d(es)/dt for the particular table configuration.
!
! Author: J. Hack
!
!-----------------------------------------------------------------------


!------------------------------Arguments--------------------------------
!
! Input arguments
!
       real(r8), intent(in) :: tmn
       real(r8), intent(in) :: tmx
       real(r8), intent(in) :: epsil
       real(r8), intent(in) :: trice
       real(r8), intent(in) :: latvap
       real(r8), intent(in) :: latice
       real(r8), intent(in) :: rh2o
       real(r8), intent(in) :: cpair
       real(r8), intent(in) :: tmeltx
!
!---------------------------Local variables-----------------------------
!
       real(r8) t
       integer n
       integer lentbl
       integer itype
!            1 -> ice phase, no transition
!           -x -> ice phase, x degree transition
       logical ip
!
!-----------------------------------------------------------------------
!
! Set es table parameters
!
       tmin   = tmn
       tmax   = tmx
       ttrice = trice
       icephs = ip
!
! Set physical constants required for es calculation
!
       epsqs = epsil
       hlatv = latvap
       hlatf = latice
       rgasv = rh2o
       cp    = cpair
       tmelt = tmeltx
!
       lentbl = INT(tmax-tmin+2.000001_r8)
       if (lentbl .gt. plenest) then



       write(*,*) "AHHH wv_sat"
       STOP
       end if
!
! Begin building es table.
! Check whether ice phase requested.
! If so, set appropriate transition range for temperature
!
       if (icephs) then
         if (ttrice /= 0.0_r8) then
           itype = -ttrice
         else
           itype = 1
         end if
       else
         itype = 0
       end if
!
       t = tmin - 1.0_r8
       do n=1,lentbl
         t = t + 1.0_r8
         call gffgch(t,estbl(n),tmelt,itype)
       end do
!
       do n=lentbl+1,plenest
         estbl(n) = -99999.0_r8
       end do
!
! Table complete -- Set coefficients for polynomial approximation of
! difference between saturation vapor press over water and saturation
! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
! is valid in the range -40 < t < 0 (degrees C).
!
!                  --- Degree 5 approximation ---
!
       pcf(1) = 5.04469588506e-01_r8
       pcf(2) = -5.47288442819e+00_r8
       pcf(3) = -3.67471858735e-01_r8
       pcf(4) = -8.95963532403e-03_r8
       pcf(5) = -7.78053686625e-05_r8
!
!                  --- Degree 6 approximation ---
!
!-----pcf(1) =  7.63285250063e-02
!-----pcf(2) = -5.86048427932e+00
!-----pcf(3) = -4.38660831780e-01
!-----pcf(4) = -1.37898276415e-02
!-----pcf(5) = -2.14444472424e-04
!-----pcf(6) = -1.36639103771e-06
!

!++jtb (comment out)
!   if (masterproc) then
!      !!write(iulog,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
!   end if
!--jtb

       return
!
9000   format('GESTBL: FATAL ERROR *********************************    &
     &',/, ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH',   &
     & ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/,      &
     & ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
!
      end subroutine gestbl

      subroutine aqsat(t ,p ,es ,qs ,ii , ilen ,kk ,kstart ,kend )
!-----------------------------------------------------------------------
!
! Purpose:
! Utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
! This routine is useful for evaluating only a selected region in the
! vertical.
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: J. Hack
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: ii
       integer,  intent(in) :: kk
       integer,  intent(in) :: ilen
       integer,  intent(in) :: kstart
       integer,  intent(in) :: kend
       real(r8), intent(in) :: t(ii,kk)
       real(r8), intent(in) :: p(ii,kk)
!
! Output arguments
!
       real(r8), intent(out) :: es(ii,kk)
       real(r8), intent(out) :: qs(ii,kk)
!
!---------------------------Local workspace-----------------------------
!
       real(r8) omeps
       integer i, k
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
       do k=kstart,kend
         do i=1,ilen
           es(i,k) = min(estblf(t(i,k)),p(i,k))
!
! Saturation specific humidity
!
           qs(i,k) = min(1.0_r8, epsqs*es(i,k)/(p(i,k)-omeps*es(i,k)))
!
! The following check is to avoid the generation of negative values
! that can occur in the upper stratosphere and mesosphere
!
!      if (qs(i,k) < 0.0_r8) then
!      qs(i,k) = 1.0_r8
!      es(i,k) = p(i,k)
!      end if

         end do
       end do
!
       return
      end subroutine aqsat

!++xl
      subroutine aqsat_water(t, p, es, qs, ii, ilen, kk, kstart,kend)
!-----------------------------------------------------------------------
!
! Purpose:
! Utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
! This routine is useful for evaluating only a selected region in the
! vertical.
!
! Method:
! <Describe the algorithm(s) used in the routine.>
! <Also include any applicable external references.>
!
! Author: J. Hack
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: ii
       integer,  intent(in) :: kk
       integer,  intent(in) :: ilen
       integer,  intent(in) :: kstart
       integer,  intent(in) :: kend
       real(r8), intent(in) :: t(ii,kk)
       real(r8), intent(in) :: p(ii,kk)
!
! Output arguments
!
       real(r8), intent(out) :: es(ii,kk)
       real(r8), intent(out) :: qs(ii,kk)
!
!---------------------------Local workspace-----------------------------
!
       real(r8) omeps
       integer i, k
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
       do k=kstart,kend
         do i=1,ilen
!          es(i,k) = estblf(t(i,k))
#ifdef GEOS5
           es(i,k) = min(polysvp(t(i,k),0), p(i,k))
#endif
#ifdef NEMS_GSM
           es(i,k) = min(fpvsl(t(i,k)), p(i,k))
#endif
!
! Saturation specific humidity
!
           qs(i,k) = min(1.0_r8, epsqs*es(i,k)/(p(i,k)-omeps*es(i,k)))
!
! The following check is to avoid the generation of negative values
! that can occur in the upper stratosphere and mesosphere
!
!          if (qs(i,k) < 0.0_r8) then
!            qs(i,k) = 1.0_r8
!            es(i,k) = p(i,k)
!          end if
         end do
       end do
!
       return
      end subroutine aqsat_water
!--xl


      subroutine aqsatd(t, p, es, qs, gam, ii, ilen, kk, kstart, kend)
!-----------------------------------------------------------------------
!
! Purpose:
! Utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g).
!
! Method:
! Differs from aqsat by also calculating and returning
! gamma (l/cp)*(d(qsat)/dT)
! Input arrays temperature and pressure (dimensioned ii,kk).
!
! Author: J. Hack
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: ii
       integer,  intent(in) :: ilen
       integer,  intent(in) :: kk
       integer,  intent(in) :: kstart
       integer,  intent(in) :: kend

       real(r8), intent(in) :: t(ii,kk)
       real(r8), intent(in) :: p(ii,kk)

!
! Output arguments
!
       real(r8), intent(out) :: es(ii,kk)
       real(r8), intent(out) :: qs(ii,kk)
       real(r8), intent(out) :: gam(ii,kk)
!
!---------------------------Local workspace-----------------------------
!
       logical lflg
       integer i
       integer k
       real(r8) omeps
       real(r8) trinv
       real(r8) tc
       real(r8) weight
       real(r8) hltalt
       real(r8) hlatsb
       real(r8) hlatvp
       real(r8) tterm
       real(r8) desdt
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
       do k=kstart,kend
         do i=1,ilen
           es(i,k) = min(p(i,k), estblf(t(i,k)))
!
! Saturation specific humidity
!
           qs(i,k) = min(1.0_r8, epsqs*es(i,k)/(p(i,k)-omeps*es(i,k)))
!
! The following check is to avoid the generation of negative qs
! values which can occur in the upper stratosphere and mesosphere
!
!
!          if (qs(i,k) < 0.0_r8) then
!            qs(i,k) = 1.0_r8
!            es(i,k) = p(i,k)
!          end if
         end do
       end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
       trinv = 0.0_r8
       if ((.not. icephs) .or. (ttrice == 0.0_r8)) go to 10
       trinv = 1.0_r8/ttrice
!
       do k=kstart,kend
         do i=1,ilen
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where constant slope is given by -2369 j/(kg c) =cpv - cw
!
           tc     = t(i,k) - tmelt
           lflg   = (tc >= -ttrice .and. tc < 0.0_r8)
           weight = min(-tc*trinv,1.0_r8)
           hlatsb = hlatv + weight*hlatf
           hlatvp = hlatv - 2369.0_r8*tc
           if (t(i,k) < tmelt) then
             hltalt = hlatsb
           else
             hltalt = hlatvp
           end if
           if (lflg) then
             tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4)      &
     &                      + tc*pcf(5))))
           else
             tterm = 0.0_r8
           end if
             desdt = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k)) + tterm*trinv
             gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/(cp*es(i,k)*(p(i,k) &
     &                  - omeps*es(i,k)))
           if (qs(i,k) == 1.0_r8) gam(i,k) = 0.0_r8
         end do
       end do
!
       go to 20
!
! No icephs or water to ice transition
!
10     do k=kstart,kend
         do i=1,ilen
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
           hlatvp = hlatv - 2369.0_r8*(t(i,k)-tmelt)
           if (icephs) then
             hlatsb = hlatv + hlatf
           else
             hlatsb = hlatv
           end if
           if (t(i,k) < tmelt) then
             hltalt = hlatsb
           else
             hltalt = hlatvp
           end if
           desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k))
           gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/(cp*es(i,k)*(p(i,k)   &
     &                  - omeps*es(i,k)))
           if (qs(i,k) == 1.0_r8) gam(i,k) = 0.0_r8
         end do
       end do
!
20     return
      end subroutine aqsatd

      subroutine vqsatd(t ,p ,es ,qs ,gam , len )
!-----------------------------------------------------------------------
!
! Purpose:
! Utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
! function as qsatd, but operates on vectors of temperature and pressure
!
! Method:
!
! Author: J. Hack
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: len
       real(r8), intent(in) :: t(len)
       real(r8), intent(in) :: p(len)
!
! Output arguments
!
       real(r8), intent(out) :: es(len)
       real(r8), intent(out) :: qs(len)
       real(r8), intent(out) :: gam(len)
!
!--------------------------Local Variables------------------------------
!
       logical lflg
!
       integer i
!
       real(r8) omeps
       real(r8) trinv
       real(r8) tc
       real(r8) weight
       real(r8) hltalt
!
       real(r8) hlatsb
       real(r8) hlatvp
       real(r8) tterm
       real(r8) desdt
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
       do i=1,len
         es(i) = min(estblf(t(i)), p(i))
!
! Saturation specific humidity
!
         qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
         qs(i) = min(1.0_r8,qs(i))
!
!        if (qs(i) < 0.0_r8) then
!          qs(i) = 1.0_r8
!          es(i) = p(i)
!        end if

       end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
       trinv = 0.0_r8
       if ((.not. icephs) .or. (ttrice.eq.0.0_r8)) go to 10
       trinv = 1.0_r8/ttrice
       do i=1,len
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
       tc = t(i) - tmelt
       lflg = (tc >= -ttrice .and. tc < 0.0_r8)
       weight = min(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       if (t(i) < tmelt) then
       hltalt = hlatsb
       else
       hltalt = hlatvp
       end if
       if (lflg) then
       tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4)            &
     &                + tc*pcf(5))))
       else
       tterm = 0.0_r8
       end if
       desdt = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_r8) gam(i) = 0.0_r8
       end do
       return
!
! No icephs or water to ice transition
!
10      do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       if (icephs) then
       hlatsb = hlatv + hlatf
       else
       hlatsb = hlatv
       end if
       if (t(i) < tmelt) then
       hltalt = hlatsb
       else
       hltalt = hlatvp
       end if
       desdt = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_r8) gam(i) = 0.0_r8
       end do
!
       return
!
      end subroutine vqsatd

!++xl
      subroutine vqsatd_water(t, p, es, qs, gam, len)

!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: len
       real(r8), intent(in) :: t(len)
       real(r8), intent(in) :: p(len)

!
! Output arguments
!
       real(r8), intent(out) :: es(len)
       real(r8), intent(out) :: qs(len)
       real(r8), intent(out) :: gam(len)

!
!--------------------------Local Variables------------------------------
!
!
       integer i
!
       real(r8) omeps
       real(r8) hltalt
!
       real(r8) hlatsb
       real(r8) hlatvp
       real(r8) desdt
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
       do i=1,len
#ifdef NEMS_GSM
         es(i) = min(fpvsl(t(i)), p(i))
#else
         es(i) = min(polysvp(t(i),0), p(i))
#endif
!
! Saturation specific humidity
!
         qs(i) = min(1.0_r8, epsqs*es(i) / (p(i)-omeps*es(i)))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
!        qs(i) = min(1.0_r8,qs(i))
!
!        if (qs(i) < 0.0_r8) then
!          qs(i) = 1.0_r8
!          es(i) = p(i)
!        end if

       end do
!
! No icephs or water to ice transition
!
       do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       hlatsb = hlatv
       if (t(i) < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_r8) gam(i) = 0.0_r8
       end do
!
       return
!
      end subroutine vqsatd_water

      function polysvp (T,typ)
!  Compute saturation vapor pressure by using
! function from Goff and Gatch (1946)
!  Polysvp returned in units of pa.
!  T is input in units of K.
!  type refers to saturation with respect to liquid (0) or ice (1)

!!DONIFF Changed to Murphy and Koop (2005) (03/04/14)


       real(r8) dum

       real(r8) t,polysvp

       integer typ


      if (.true.) then
!ice
       if (typ == 1) then
         polysvp = MurphyKoop_svp_ice(t)
       end if
       if (typ == 0) then
         polysvp = MurphyKoop_svp_water(t)
       end if

      else

! ice
       if (typ.eq.1) then



       polysvp = 10._r8**(-9.09718_r8*(273.16_r8/t-1._r8)-3.56654_r8*   &
     & log10(273.16_r8/t)+0.876793_r8*(1._r8-t/273.16_r8)+              &
     & log10(6.1071_r8))*100._r8

       end if



       if (typ.eq.0) then
       polysvp = 10._r8**(-7.90298_r8*(373.16_r8/t-1._r8)+ 5.02808_r8*   &
     &log10(373.16_r8/t)- 1.3816e-7_r8*(10._r8**(11.344_r8*(1._r8-t/     &
     &373.16_r8))-1._r8)+ 8.1328e-3_r8*(10._r8**(-3.49149_r8*(373.16_r8/ &
     &t-1._r8))-1._r8)+ log10(1013.246_r8))*100._r8
       end if

       end if

       end function polysvp
!--xl



      integer function fqsatd(t ,p ,es ,qs ,gam , len )





       integer, intent(in) :: len
       real(r8), intent(in) :: t(len)
       real(r8), intent(in) :: p(len)

       real(r8), intent(out) :: es(len)
       real(r8), intent(out) :: qs(len)
       real(r8), intent(out) :: gam(len)

       call vqsatd(t ,p ,es ,qs ,gam , len )
       fqsatd = 1
       return
      end function fqsatd

      real(r8) function qsat_water(t,p)

       real(r8) t
       real(r8) p
       real(r8) es
       real(r8) ps, ts, e1, e2, f1, f2, f3, f4, f5, f





       ps = 1013.246_r8
       ts = 373.16_r8
       e1 = 11.344_r8*(1.0_r8 - t/ts)
       e2 = -3.49149_r8*(ts/t - 1.0_r8)
       f1 = -7.90298_r8*(ts/t - 1.0_r8)
       f2 = 5.02808_r8*log10(ts/t)
       f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
       f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
       f5 = log10(ps)
       f = f1 + f2 + f3 + f4 + f5
       es = (10.0_r8**f)*100.0_r8

       qsat_water = epsqs*es/(p-(1.-epsqs)*es)
       if(qsat_water < 0.) qsat_water = 1.

      end function qsat_water

      subroutine vqsat_water(t,p,qsat_water,len)

       integer, intent(in) :: len
       real(r8) t(len)
       real(r8) p(len)
       real(r8) qsat_water(len)
       real(r8) es
       real(r8), parameter :: t0inv = 1._r8/273._r8
       real(r8) coef
       integer :: i

       coef = hlatv/rgasv
       do i=1,len
       es = 611._r8*exp(coef*(t0inv-1./t(i)))
       qsat_water(i) = epsqs*es/(p(i)-(1.-epsqs)*es)
       if(qsat_water(i) < 0.) qsat_water(i) = 1.
       enddo

       return

      end subroutine vqsat_water

      real(r8) function qsat_ice(t,p)

       real(r8) t
       real(r8) p
       real(r8) es
       real(r8), parameter :: t0inv = 1._r8/273._r8
       es = 611.*exp((hlatv+hlatf)/rgasv*(t0inv-1./t))
       qsat_ice = epsqs*es/(p-(1.-epsqs)*es)
       if(qsat_ice < 0.) qsat_ice = 1.

      end function qsat_ice

      subroutine vqsat_ice(t,p,qsat_ice,len)

       integer,intent(in) :: len
       real(r8) t(len)
       real(r8) p(len)
       real(r8) qsat_ice(len)
       real(r8) es
       real(r8), parameter :: t0inv = 1._r8/273._r8
       real(r8) coef
       integer :: i

       coef = (hlatv+hlatf)/rgasv
       do i=1,len
       es = 611.*exp(coef*(t0inv-1./t(i)))
       qsat_ice(i) = epsqs*es/(p(i)-(1.-epsqs)*es)
       if(qsat_ice(i) < 0.) qsat_ice(i) = 1.
       enddo

       return

      end subroutine vqsat_ice

! Sungsu
! Below two subroutines (vqsatd2_water,vqsatd2_water_single) are by Sungsu
! Replace 'gam -> dqsdt'
! Sungsu

      subroutine vqsatd2_water(t ,p ,es ,qs ,dqsdt , len )

!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer, intent(in) :: len
       real(r8), intent(in) :: t(len)
       real(r8), intent(in) :: p(len)

!
! Output arguments
!
       real(r8), intent(out) :: es(len)
       real(r8), intent(out) :: qs(len)


       real(r8), intent(out) :: dqsdt(len)


!
!--------------------------Local Variables------------------------------
!
!
       integer i
!
       real(r8) omeps
       real(r8) hltalt
!
       real(r8) hlatsb
       real(r8) hlatvp
       real(r8) desdt


       real(r8) gam(len)


!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
       do i=1,len
#ifdef GEOS5
       es(i) = min(polysvp(t(i),0), p(i))
#endif
#ifdef NEMS_GSM
       es(i) = min(fpvsl(t(i)), p(i))
#endif
!
! Saturation specific humidity
!
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
       qs(i) = min(1.0_r8,qs(i))
!
!      if (qs(i) < 0.0_r8) then
!        qs(i) = 1.0_r8
!        es(i) = p(i)
!      end if

       end do
!
! No icephs or water to ice transition
!
       do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       hlatsb = hlatv
       if (t(i) < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_r8) gam(i) = 0.0_r8

       dqsdt(i) = (cp/hltalt)*gam(i)

       end do
!
       return
!
      end subroutine vqsatd2_water

      subroutine vqsatd2_water_single(t ,p ,es ,qs ,dqsdt)

!------------------------------Arguments--------------------------------
!
! Input arguments
!

       real(r8), intent(in) :: t, p

!
! Output arguments
!
       real(r8), intent(out) :: es, qs, dqsdt
!
!--------------------------Local Variables------------------------------
!
!      integer i
!
       real(r8) omeps, hltalt, hlatsb, hlatvp, desdt, gam
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
!  do i=1,len
#ifdef GEOS5
       es = min(p, polysvp(t,0))
#endif
#ifdef NEMS_GSM
       es = min(p, fpvsl(t))
#endif
!
! Saturation specific humidity
!
       qs = min(1.0_r8, epsqs*es/(p-omeps*es))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
!      if (qs < 0.0_r8) then
!      qs = 1.0_r8
!      es = p
!      end if
!  end do
!
! No icephs or water to ice transition
!
!  do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_r8*(t-tmelt)
       hlatsb = hlatv
       if (t < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt = hltalt*es/(rgasv*t*t)
       gam   = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
       if (qs >= 1.0_r8) gam = 0.0_r8

       dqsdt = (cp/hltalt)*gam

!  end do
!
       return
!
      end subroutine vqsatd2_water_single


      subroutine vqsatd2(t ,p ,es ,qs ,dqsdt , len )
!-----------------------------------------------------------------------
! Sungsu : This is directly copied from 'vqsatd' but 'dqsdt' is output
!          instead of gam for use in Sungsu's equilibrium stratiform
!          macrophysics scheme.
!
! Purpose:
! Utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
! function as qsatd, but operates on vectors of temperature and pressure
!
! Method:
!
! Author: J. Hack
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer, intent(in) :: len
       real(r8), intent(in) :: t(len)
       real(r8), intent(in) :: p(len)
!
! Output arguments
!
       real(r8), intent(out) :: es(len)
       real(r8), intent(out) :: qs(len)


       real(r8), intent(out) :: dqsdt(len)


!
!--------------------------Local Variables------------------------------
!
       logical lflg
!
       integer i
!
       real(r8) omeps
       real(r8) trinv
       real(r8) tc
       real(r8) weight
       real(r8) hltalt
!
       real(r8) hlatsb
       real(r8) hlatvp
       real(r8) tterm
       real(r8) desdt


       real(r8) gam(len)

!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
       do i=1,len
#ifdef GEOS5
         es(i) = min(p(i), estblf(t(i)))
#endif
#ifdef NEMS_GSM
         es(i) = min(p(i), fpvsi(t(i)))
#endif
!
! Saturation specific humidity
!
         qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
         qs(i) = min(1.0_r8,qs(i))
!
!        if (qs(i) < 0.0_r8) then
!        qs(i) = 1.0_r8
!        es(i) = p(i)
!        end if
       end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
       trinv = 0.0_r8
       if ((.not. icephs) .or. (ttrice == 0.0_r8)) go to 10
       trinv = 1.0_r8/ttrice
       do i=1,len
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
       tc = t(i) - tmelt
       lflg = (tc >= -ttrice .and. tc < 0.0_r8)
       weight = min(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       if (t(i) < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       if (lflg) then
         tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4)          &
     &                  + tc*pcf(5))))
       else
         tterm = 0.0_r8
       end if
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_r8) gam(i) = 0.0_r8

       dqsdt(i) = (cp/hltalt)*gam(i)

       end do
       return
!
! No icephs or water to ice transition
!
10     do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       if (icephs) then
       hlatsb = hlatv + hlatf
       else
       hlatsb = hlatv
       end if
       if (t(i) < tmelt) then
       hltalt = hlatsb
       else
       hltalt = hlatvp
       end if
       desdt = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_r8) gam(i) = 0.0_r8

       dqsdt(i) = (cp/hltalt)*gam(i)

       end do
!
       return
!
      end subroutine vqsatd2


! Below routine is by Sungsu

      subroutine vqsatd2_single(t ,p ,es ,qs ,dqsdt)
!-----------------------------------------------------------------------
! Sungsu : This is directly copied from 'vqsatd' but 'dqsdt' is output
!          instead of gam for use in Sungsu's equilibrium stratiform
!          macrophysics scheme.
!
! Purpose:
! Utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
! function as qsatd, but operates on vectors of temperature and pressure
!
! Method:
!
! Author: J. Hack
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       real(r8), intent(in) :: t, p
!
! Output arguments
!
       real(r8), intent(out) :: es, qs, dqsdt
!
!--------------------------Local Variables------------------------------
!
       logical lflg
!
!  integer i      ! index for vector calculations
!
       real(r8) omeps
       real(r8) trinv
       real(r8) tc
       real(r8) weight
       real(r8) hltalt
!
       real(r8) hlatsb
       real(r8) hlatvp
       real(r8) tterm
       real(r8) desdt

       real(r8) gam

!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs

!  do i=1,len

#ifdef GEOS5
       es = estblf(t)
#endif
#ifdef NEMS_GSM
       es = min(fpvs(t), p)
#endif
!
! Saturation specific humidity
!
       qs = epsqs*es/(p - omeps*es)
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
       qs = min(1.0_r8,qs)
!
!      if (qs < 0.0_r8) then
!        qs = 1.0_r8
!        es = p
!      end if

!  end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
       trinv = 0.0_r8
       if ((.not. icephs) .or. (ttrice == 0.0_r8)) go to 10
       trinv = 1.0_r8/ttrice

!  do i=1,len
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
       tc = t - tmelt
       lflg = (tc >= -ttrice .and. tc < 0.0_r8)
       weight = min(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       if (t < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       if (lflg) then
         tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4)          &
     &                  + tc*pcf(5))))
       else
         tterm = 0.0_r8
       end if
       desdt = hltalt*es/(rgasv*t*t) + tterm*trinv
       gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
       if (qs == 1.0_r8) gam = 0.0_r8

       dqsdt = (cp/hltalt)*gam

!  end do
       return
!
! No icephs or water to ice transition
!

10     continue

!10 do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_r8*(t-tmelt)
       if (icephs) then
         hlatsb = hlatv + hlatf
       else
         hlatsb = hlatv
       end if
       if (t < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt = hltalt*es/(rgasv*t*t)
       gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
       if (qs == 1.0_r8) gam = 0.0_r8

       dqsdt = (cp/hltalt)*gam


!  end do
!
       return
!
      end subroutine vqsatd2_single

!----------------------------------------------------------------------

!----------------------------------------------------------------------

      subroutine gffgch(t ,es ,tmelt ,itype )
!-----------------------------------------------------------------------
!
! Purpose:
! Computes saturation vapor pressure over water and/or over ice using
! Goff & Gratch (1946) relationships.
! <Say what the routine does>
!
! Method:
! T (temperature), and itype are input parameters, while es (saturation
! vapor pressure) is an output parameter.  The input parameter itype
! serves two purposes: a value of zero indicates that saturation vapor
! pressures over water are to be returned (regardless of temperature),
! while a value of one indicates that saturation vapor pressures over
! ice should be returned when t is less than freezing degrees.  If itype
! is negative, its absolute value is interpreted to define a temperature
! transition region below freezing in which the returned
! saturation vapor pressure is a weighted average of the respective ice
! and water value.  That is, in the temperature range 0 => -itype
! degrees c, the saturation vapor pressures are assumed to be a weighted
! average of the vapor pressure over supercooled water and ice (all
! water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
!
! Author: J. Hack
!
!-----------------------------------------------------------------------





       implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       real(r8), intent(in) :: t ,tmelt
!
! Output arguments
!
       integer, intent(inout) :: itype

       real(r8), intent(out) :: es
!
!---------------------------Local variables-----------------------------
!
       real(r8) e1
       real(r8) e2
       real(r8) eswtr
       real(r8) f
       real(r8) f1
       real(r8) f2
       real(r8) f3
       real(r8) f4
       real(r8) f5
       real(r8) ps
       real(r8) t0
       real(r8) term1
       real(r8) term2
       real(r8) term3
       real(r8) tr
       real(r8) ts
       real(r8) weight
       integer itypo
!
!-----------------------------------------------------------------------
!
! Check on whether there is to be a transition region for es
!
       if (itype < 0) then
       tr = abs(real(itype,r8))
       itypo = itype
       itype = 1
       else
       tr = 0.0_r8
       itypo = itype
       end if
       if (tr > 40.0_r8) then
       write(iulog,900) tr

       end if
!
       if(t < (tmelt - tr) .and. itype == 1) go to 10
!
! Water
!
       ps = 1013.246_r8
       ts = 373.16_r8
       e1 = 11.344_r8*(1.0_r8 - t/ts)
       e2 = -3.49149_r8*(ts/t - 1.0_r8)
       f1 = -7.90298_r8*(ts/t - 1.0_r8)
       f2 = 5.02808_r8*log10(ts/t)
       f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
       f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
       f5 = log10(ps)
       f = f1 + f2 + f3 + f4 + f5
       es = (10.0_r8**f)*100.0_r8
       eswtr = es
!
       if(t >= tmelt .or. itype == 0) go to 20
!
! Ice
!
10     continue
       t0 = tmelt
       term1 = 2.01889049_r8/(t0/t)
       term2 = 3.56654_r8*log(t0/t)
       term3 = 20.947031_r8*(t0/t)
       es = 575.185606e10_r8*exp(-(term1 + term2 + term3))
!
       if (t < (tmelt - tr)) go to 20
!
! Weighted transition between water and ice
!
       weight = min((tmelt - t)/tr,1.0_r8)
       es = weight*es + (1.0_r8 - weight)*eswtr
!
20     continue
       itype = itypo
       return
!
900    format('GFFGCH: FATAL ERROR ******************************',/,    &
     & 'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', ' PRESSURE, &
     & TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', ' 40.0 DEGREES C',/,     &
     & ' TR = ',f7.2)
!
      end subroutine gffgch


!!DONIF USe Murphy and Koop (2005) (Written by Andrew Gettelman)

       function MurphyKoop_svp_water(tx) result(es)
       real(r8), intent(in) :: tx
       real(r8) :: es
       real(r8):: t

       t=min(tx, 332.0_r8)
       t=max(123.0_r8, tx)

       es = exp(54.842763_r8 - (6763.22_r8 / t) - (4.210_r8 * log(t)) + &
     & (0.000367_r8 * t) + (tanh(0.0415_r8 * (t - 218.8_r8)) *          &
     & (53.878_r8 - (1331.22_r8 / t) - (9.44523_r8 * log(t)) +          &
     & 0.014025_r8 * t)))

      end function MurphyKoop_svp_water

       function MurphyKoop_svp_ice(tx) result(es)
       real(r8), intent(in) :: tx
       real(r8) :: t
       real(r8) :: es

       t=max(100.0_r8, tx)
       t=min(274.0_r8, tx)


       es = exp(9.550426_r8 - (5723.265_r8 / t) + (3.53068_r8 *         &
     & log(t)) - (0.00728332_r8 * t))

      end function MurphyKoop_svp_ice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine vqsatd2_ice_single(t ,p ,es ,qs ,dqsdt)

!------------------------------Arguments--------------------------------
!
! Input arguments
!
       real(r8), intent(in) :: t, p
!
! Output arguments
!
       real(r8), intent(out) :: es, qs, dqsdt

!
!--------------------------Local Variables------------------------------
!
!      integer i
!
       real(r8) omeps, hltalt, hlatsb, hlatvp, desdt, gam
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_r8 - epsqs
!  do i=1,len
#ifdef GEOS5
       es = min(polysvp(t,1),p)
#endif
#ifdef NEMS_GSM
       es = min(fpvsi(t),p)
#endif
!
! Saturation specific humidity
!
       qs = min(1.0_r8, epsqs*es/(p-omeps*es))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
!      if (qs < 0.0_r8) then
!        qs = 1.0_r8
!        es = p
!      end if
!  end do
!
! No icephs or water to ice transition
!
!  do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hltalt = hlatv + hlatf
       desdt  = hltalt*es/(rgasv*t*t)
       if (qs < 1.0_r8) then
         gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
       else
         gam = 0.0_r8
       endif

       dqsdt = (cp/hltalt)*gam

!  end do
!
       return
!
      end subroutine vqsatd2_ice_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module wv_saturation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
