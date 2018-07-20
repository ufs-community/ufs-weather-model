!>  \file mfpbl.f
!!  This file contains the subroutine that calculates the updraft properties and mass flux for use in the Hybrid EDMF PBL scheme.

!>  \ingroup GFS_edmf_main
!!  \brief This subroutine is used for calculating the mass flux and updraft properties.
!!
!!  The mfpbl routines works as follows: if the PBL is convective, first, the ascending parcel entrainment rate is calculated as a
!!  function of height. Next, a surface parcel is initiated according to surface layer properties and the updraft buoyancy is calculated
!!  as a function of height. Next, using the buoyancy and entrainment values, the parcel vertical velocity is calculated using a well
!!  known steady-state budget equation. With the profile of updraft vertical velocity, the PBL height is recalculated as the height
!!  where the updraft vertical velocity returns to 0, and the entrainment profile is updated with the new PBL height. Finally, the mass
!!  flux profile is calculated using the updraft vertical velocity and assumed updraft fraction and the updraft properties are calculated
!!  using the updated entrainment profile, surface values, and environmental profiles.
!!  \param[in] im      integer, number of used points
!!  \param[in] ix      integer, horizontal dimension
!!  \param[in] km      integer, vertical layer dimension
!!  \param[in] ntrac   integer, number of tracers
!!  \param[in] delt    real, physics time step
!!  \param[in] cnvflg  logical, im, flag to denote a strongly unstable (convective) PBL
!!  \param[in] zl      real, (im, km), height of grid centers
!!  \param[in] zm      real, (im, km+1), height of grid interfaces
!!  \param[in] thvx    real, (im, km), virtual potential temperature at grid centers (\f$ K \f$)
!!  \param[in] q1      real, (ix, km, ntrac), layer mean tracer concentration (units?)
!!  \param[in] t1      real, (ix, km), layer mean temperature (\f$ K \f$)
!!  \param[in] u1      real, (ix, km), u component of layer wind (\f$ m s^{-1} \f$)
!!  \param[in] v1      real, (ix, km), v component of layer wind (\f$ m s^{-1} \f$)
!!  \param[in,out] hpbl  real, im, PBL top height (m)
!!  \param[in,out] kpbl  integer, im, PBL top index
!!  \param[in] sflx      real, im, total surface heat flux (units?)
!!  \param[in] ustar   real, im, surface friction velocity
!!  \param[in] wstar   real, im, convective velocity scale
!!  \param[out] xmf    real, (im, km), updraft mass flux
!!  \param[in,out] tcko  real, (im, km), updraft temperature (\f$ K \f$)
!!  \param[in,out] qcko  real, (im, km, ntrac), updraft tracer concentration (units?)
!!  \param[in,out] ucko  real, (im, km), updraft u component of horizontal momentum (\f$ m s^{-1} \f$)
!!  \param[in,out] vcko  real, (im, km), updraft v component of horizontal momentum (\f$ m s^{-1} \f$)
!!
!!  \section general_mfpbl mfpbl General Algorithm
!!  -# Determine an updraft parcel's entrainment rate, buoyancy, and vertical velocity.
!!  -# Recalculate the PBL height (previously calculated in moninedmf) and the parcel's entrainment rate.
!!  -# Calculate the mass flux profile and updraft properties.
!!  \section detailed_mfpbl mfpbl Detailed Algorithm
!!  @{
      subroutine mfpbl(im,ix,km,ntrac,delt,cnvflg,                      &
     &   zl,zm,thvx,q1,t1,u1,v1,hpbl,kpbl,                              &
     &   sflx,ustar,wstar,xmf,tcko,qcko,ucko,vcko)
!
      use machine , only : kind_phys
      use physcons, grav => con_g, cp => con_cp
!
      implicit none
!
      integer              im, ix, km, ntrac
!    &,                    me
      integer              kpbl(im)
      logical              cnvflg(im)
      real(kind=kind_phys) delt
      real(kind=kind_phys) q1(ix,km,ntrac), t1(ix,km),                  &
     &                     u1(ix,km),  v1(ix,km),                       &
     &                     thvx(im,km),                                 &
     &                     zl(im,km),  zm(im,km+1),                     &
     &                     hpbl(im),   sflx(im),    ustar(im),          &
     &                     wstar(im),  xmf(im,km),                      &
     &                     tcko(im,km),qcko(im,km,ntrac),               &
     &                     ucko(im,km),vcko(im,km)
!
c  local variables and arrays
!
      integer   i, j, k, n, kmpbl
!
      real(kind=kind_phys) dt2,     dz,      ce0,
     &                     h1,      factor,  gocp,
     &                     g,       c1,      d1,
     &                     b1,      f1,      bb1,     bb2,
     &                     alp,     a1,      qmin,    zfmin,
     &                     xmmx,    rbint,   tau,
!    &                     rbint,   tau,
     &                     tem,     tem1,    tem2,
     &                     ptem,    ptem1,   ptem2,
     &                     pgcon
!
      real(kind=kind_phys) sigw1(im),   usws3(im),  xlamax(im),
     &                     rbdn(im),    rbup(im),   delz(im)
!
      real(kind=kind_phys) wu2(im,km),     xlamue(im,km),
     &                     thvu(im,km),    zi(im,km),
     &                     buo(im,km)
!
      logical totflg, flg(im)
!
c  physical parameters
      parameter(g=grav)
      parameter(gocp=g/cp)
!     parameter(ce0=0.37,qmin=1.e-8,alp=1.0,pgcon=0.55)
      parameter(ce0=0.38,qmin=1.e-8,alp=1.0,pgcon=0.55)
      parameter(a1=0.08,b1=0.5,f1=0.15,c1=0.3,d1=2.58,tau=500.)
      parameter(zfmin=1.e-8,h1=0.33333333)
!
c-----------------------------------------------------------------------
!
!************************************************************************
!
      kmpbl = km/2 + 1
      dt2 = delt
!> Since the mfpbl subroutine is called regardless of whether the PBL is convective, a check of the convective PBL flag is performed and the subroutine returns back to moninedmf (with the output variables set to the initialized values) if the PBL is not convective.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
      do k = 1, km
        do i=1,im
          if (cnvflg(i)) then
            zi(i,k) = zm(i,k+1)
          endif
        enddo
      enddo
!>  ## Determine an updraft parcel's entrainment rate, buoyancy, and vertical velocity.
!!  Calculate the entrainment rate according to equation 16 in \cite siebesma_et_al_2007 for all levels (xlamue) and a default entrainment rate (xlamax) for use above the PBL top.
      do i=1,im
        if(cnvflg(i)) then
          k = kpbl(i) / 2
          k = max(k, 1)
          delz(i) = zl(i,k+1) - zl(i,k)
          xlamax(i) = ce0 / delz(i)
        endif
      enddo
      do k = 1, kmpbl
        do i=1,im
          if(cnvflg(i)) then
            if(k < kpbl(i)) then
              ptem = 1./(zi(i,k)+delz(i))
              tem = max((hpbl(i)-zi(i,k)+delz(i)) ,delz(i))
              ptem1 = 1./tem
              xlamue(i,k) = ce0 * (ptem+ptem1)
            else
              xlamue(i,k) = xlamax(i)
            endif
          endif
        enddo
      enddo
c
c  compute thermal excess
c
!>  Using equations 17 and 7 from \cite siebesma_et_al_2007 along with \f$u_*\f$, \f$w_*\f$, and the previously diagnosed PBL height, the initial \f$\theta_v\f$ of the updraft (and its surface buoyancy) is calculated.
      do i=1,im
        if(cnvflg(i)) then
          tem = zl(i,1)/hpbl(i)
          usws3(i) = (ustar(i)/wstar(i))**3.
          tem1 = usws3(i) + 0.6*tem
          tem2 = max((1.-tem), zfmin)
          ptem = (tem1**h1) * sqrt(tem2)
          sigw1(i) = 1.3 * ptem * wstar(i)
          ptem1 = alp * sflx(i) / sigw1(i)
          thvu(i,1) = thvx(i,1) + ptem1
          buo(i,1) = g * (thvu(i,1)/thvx(i,1)-1.)
        endif
      enddo
c
c  compute potential temperature and buoyancy for updraft air parcel
c
!>  From the second level to the middle of the vertical domain, the updraft virtual potential temperature is calculated using the entraining updraft equation as in equation 10 of \cite siebesma_et_al_2007, discretized as
!!  \f[
!!  \frac{\theta_{v,u}^k - \theta_{v,u}^{k-1}}{\Delta z}=-\epsilon^{k-1}\left[\frac{1}{2}\left(\theta_{v,u}^k + \theta_{v,u}^{k-1}\right)-\frac{1}{2}\left(\overline{\theta_{v}}^k + \overline{\theta_v}^{k-1}\right)\right]
!!  \f]
!!  where the superscript \f$k\f$ denotes model level, and subscript \f$u\f$ denotes an updraft property, and the overbar denotes the grid-scale mean value.
      do k = 2, kmpbl
        do i=1,im
          if(cnvflg(i)) then
            dz = zl(i,k) - zl(i,k-1)
            tem = xlamue(i,k-1) * dz
            ptem = 2. + tem
            ptem1 = (2. - tem) / ptem
            tem1 = tem  * (thvx(i,k)+thvx(i,k-1)) / ptem
            thvu(i,k) = ptem1 * thvu(i,k-1) + tem1
            buo(i,k) = g * (thvu(i,k)/thvx(i,k)-1.)
          endif
        enddo
      enddo
c
c  compute updraft velocity square(wu2)
c
!>  Rather than use the vertical velocity equation given as equation 15 in \cite siebesma_et_al_2007 (which parameterizes the pressure term in terms of the updraft vertical velocity itself), this scheme uses the more widely used form of the steady state vertical velocity equation given as equation 6 in \cite soares_et_al_2004 discretized as
!!  \f[
!!  \frac{w_{u,k}^2 - w_{u,k-1}^2}{\Delta z} = -2b_1\frac{1}{2}\left(\epsilon_k + \epsilon_{k-1}\right)\frac{1}{2}\left(w_{u,k}^2 + w_{u,k-1}^2\right) + 2b_2B
!!  \f]
!! The constants used in the scheme are labeled \f$bb1 = 2b_1\f$ and \f$bb2 = 2b_2\f$ and are tuned to be equal to 1.8 and 3.5, respectively, close to the values proposed by \cite soares_et_al_2004 .
!     tem = 1.-2.*f1
!     bb1 = 2. * b1 / tem
!     bb2 = 2. / tem
!  from soares et al. (2004,qjrms)
!     bb1 = 2.
!     bb2 = 4.
!
!  from bretherton et al. (2004, mwr)
!     bb1 = 4.
!     bb2 = 2.
!
!  from our tuning
      bb1 = 1.8
      bb2 = 3.5
!
      do i = 1, im
        if(cnvflg(i)) then
!
!         tem = zi(i,1)/hpbl(i)
!         tem1 = usws3(i) + 0.6*tem
!         tem2 = max((1.-tem), zfmin)
!         ptem = (tem1**h1) * sqrt(tem2)
!         ptem1 = 1.3 * ptem * wstar(i)
!         wu2(i,1) = d1*d1*ptem1*ptem1
!
          dz   = zi(i,1)
          tem  = 0.5*bb1*xlamue(i,1)*dz
          tem1 = bb2 * buo(i,1) * dz
          ptem1 = 1. + tem
          wu2(i,1) = tem1 / ptem1
!
        endif
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if(cnvflg(i)) then
            dz    = zi(i,k) - zi(i,k-1)
            tem  = 0.25*bb1*(xlamue(i,k)+xlamue(i,k-1))*dz
            tem1 = bb2 * buo(i,k) * dz
            ptem = (1. - tem) * wu2(i,k-1)
            ptem1 = 1. + tem
            wu2(i,k) = (ptem + tem1) / ptem1
          endif
        enddo
      enddo
c
c  update pbl height as the height where updraft velocity vanishes
c
!>  ## Recalculate the PBL height and the parcel's entrainment rate.
!!  Find the level where the updraft vertical velocity is less than zero and linearly interpolate to find the height where it would be exactly zero. Set the PBL height to this determined height.
      do i=1,im
         flg(i)  = .true.
         if(cnvflg(i)) then
           flg(i)  = .false.
           rbup(i) = wu2(i,1)
         endif
      enddo
      do k = 2, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          rbup(i) = wu2(i,k)
          kpbl(i) = k
          flg(i)  = rbup(i).le.0.
        endif
      enddo
      enddo
      do i = 1,im
        if(cnvflg(i)) then
           k = kpbl(i)
           if(rbdn(i) <= 0.) then
              rbint = 0.
           elseif(rbup(i) >= 0.) then
              rbint = 1.
           else
              rbint = rbdn(i)/(rbdn(i)-rbup(i))
           endif
           hpbl(i) = zi(i,k-1) + rbint*(zi(i,k)-zi(i,k-1))
        endif
      enddo
c
!>  Recalculate the entrainment rate as before except use the updated value of the PBL height.
      do i=1,im
        if(cnvflg(i)) then
          k = kpbl(i) / 2
          k = max(k, 1)
          delz(i) = zl(i,k+1) - zl(i,k)
          xlamax(i) = ce0 / delz(i)
        endif
      enddo
c
c  update entrainment rate
c
!     do k = 1, kmpbl
!       do i=1,im
!         if(cnvflg(i)) then
!           if(k < kpbl(i)) then
!             tem = tau * sqrt(wu2(i,k))
!             tem1 = 1. / tem
!             ptem = ce0 / zi(i,k)
!             xlamue(i,k) = max(tem1, ptem)
!           else
!             xlamue(i,k) = xlamax(i)
!           endif
!         endif
!       enddo
!     enddo
!
      do k = 1, kmpbl
        do i=1,im
          if(cnvflg(i)) then
            if(k < kpbl(i)) then
              ptem = 1./(zi(i,k)+delz(i))
              tem = max((hpbl(i)-zi(i,k)+delz(i)) ,delz(i))
              ptem1 = 1./tem
              xlamue(i,k) = ce0 * (ptem+ptem1)
            else
              xlamue(i,k) = xlamax(i)
            endif
          endif
        enddo
      enddo
c
c  updraft mass flux as a function of sigmaw
c   (0.3*sigmaw[square root of vertical turbulence variance])
c
!>  ## Calculate the mass flux profile and updraft properties.
!     do k = 1, kmpbl
!       do i=1,im
!         if(cnvflg(i) .and. k < kpbl(i)) then
!           tem = zi(i,k)/hpbl(i)
!           tem1 = usws3(i) + 0.6*tem
!           tem2 = max((1.-tem), zfmin)
!           ptem = (tem1**h1) * sqrt(tem2)
!           ptem1 = 1.3 * ptem * wstar(i)
!           xmf(i,k) = c1 * ptem1
!         endif
!       enddo
!     enddo
c
c  updraft mass flux as a function of updraft velocity profile
c
!>  Calculate the mass flux:
!!  \f[
!!  M = a_uw_u
!!  \f]
!!  where \f$a_u\f$ is the tunable parameter that represents the fractional area of updrafts (currently set to 0.08). Limit the computed mass flux to be less than \f$\frac{\Delta z}{\Delta t}\f$. This is different than what is done in \cite siebesma_et_al_2007 where the mass flux is the product of a tunable constant and the diagnosed standard deviation of \f$w\f$.
      do k = 1, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k < kpbl(i)) then
             xmf(i,k) = a1 * sqrt(wu2(i,k))
             dz   = zl(i,k+1) - zl(i,k)
             xmmx = dz / dt2
             xmf(i,k) = min(xmf(i,k),xmmx)
          endif
        enddo
      enddo
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  compute updraft property
c
!>  The updraft properties are calculated according to the entraining updraft equation
!!  \f[
!!  \frac{\partial \phi}{\partial z}=-\epsilon\left(\phi_u - \overline{\phi}\right)
!!  \f]
!!  where \f$\phi\f$ is \f$T\f$ or \f$q\f$. The equation is discretized according to
!!  \f[
!!  \frac{\phi_{u,k} - \phi_{u,k-1}}{\Delta z}=-\epsilon_{k-1}\left[\frac{1}{2}\left(\phi_{u,k} + \phi_{u,k-1}\right)-\frac{1}{2}\left(\overline{\phi}_k + \overline{\phi}_{k-1}\right)\right]
!!  \f]
!!  The exception is for the horizontal momentum components, which have been modified to account for the updraft-induced pressure gradient force, and use the following equation, following \cite han_and_pan_2006
!!  \f[
!!  \frac{\partial v}{\partial z} = -\epsilon\left(v_u - \overline{v}\right)+d_1\frac{\partial \overline{v}}{\partial z}
!!  \f]
!!  where \f$d_1=0.55\f$ is a tunable constant.
      do k = 2, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k <= kpbl(i)) then
             dz   = zl(i,k) - zl(i,k-1)
             tem  = 0.5 * xlamue(i,k-1) * dz
             factor = 1. + tem
             ptem = tem + pgcon
             ptem1= tem - pgcon
!
             tcko(i,k) = ((1.-tem)*tcko(i,k-1)+tem*
     &                    (t1(i,k)+t1(i,k-1))-gocp*dz)/factor
             ucko(i,k) = ((1.-tem)*ucko(i,k-1)+ptem*u1(i,k)
     &                    +ptem1*u1(i,k-1))/factor
             vcko(i,k) = ((1.-tem)*vcko(i,k-1)+ptem*v1(i,k)
     &                    +ptem1*v1(i,k-1))/factor
          endif
        enddo
      enddo
      do n = 1, ntrac
      do k = 2, kmpbl
        do i = 1, im
          if (cnvflg(i) .and. k <= kpbl(i)) then
             dz   = zl(i,k) - zl(i,k-1)
             tem  = 0.5 * xlamue(i,k-1) * dz
             factor = 1. + tem

             qcko(i,k,n) = ((1.-tem)*qcko(i,k-1,n)+tem*
     &                    (q1(i,k,n)+q1(i,k-1,n)))/factor
          endif
        enddo
      enddo
      enddo
!
      return
      end
!> @}
