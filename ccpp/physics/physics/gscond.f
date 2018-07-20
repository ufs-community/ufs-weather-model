!> \file gscond.f
!! This file contains the subroutine that calculates grid-scale
!! condensation and evaporation for use in 
!! \cite zhao_and_carr_1997 scheme.

      module zhaocarr_gscond
      contains


! \brief Brief description of the subroutine
!
!> \section arg_table_gscond_init  Argument Table
!!
       subroutine zhaocarr_gscond_init
       end subroutine zhaocarr_gscond_init

! \brief Brief description of the subroutine
!
!> \section arg_table_gscond_finalize  Argument Table
!!
       subroutine zhaocarr_gscond_finalize
       end subroutine zhaocarr_gscond_finalize

!> \defgroup condense GFS gscond Main
!! \brief This subroutine computes grid-scale condensation and evaporation of
!! cloud condensate.
!!
!! \section arg_table_zhaocarr_gscond_run Argument Table
!! | local_name     | standard_name                                              | long_name                                                | units   | rank |  type     |   kind    | intent | optional |
!! |----------------|------------------------------------------------------------|----------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                     | horizontal loop extent                                   | count   |    0 | integer   |           | in     | F        |
!! | ix             | horizontal_dimension                                       | horizontal dimension                                     | count   |    0 | integer   |           | in     | F        |
!! | km             | vertical_dimension                                         | vertical layer dimension                                 | count   |    0 | integer   |           | in     | F        |
!! | dt             | time_step_for_physics                                      | physics time step                                        | s       |    0 | real      | kind_phys | in     | F        |
!! | dtf            | time_step_for_dynamics                                     | dynamics time step                                       | s       |    0 | real      | kind_phys | in     | F        |
!! | prsl           | air_pressure                                               | layer mean air pressure                                  | Pa      |    2 | real      | kind_phys | in     | F        |
!! | ps             | surface_air_pressure                                       | surface pressure                                         | Pa      |    1 | real      | kind_phys | in     | F        |
!! | q              | water_vapor_specific_humidity_updated_by_physics           | water vapor specific humidity                            | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | clw1           | cloud_ice_specific_humidity                                | cloud ice specific humidity                              | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | clw2           | cloud_liquid_water_specific_humidity                       | cloud water specific humidity                            | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | cwm            | cloud_condensed_water_specific_humidity_updated_by_physics | cloud condensed water specific humidity                  | kg kg-1 |    2 | real      | kind_phys | out    | F        |
!! | t              | air_temperature_updated_by_physics                         | layer mean air temperature                               | K       |    2 | real      | kind_phys | inout  | F        |
!! | tp             | air_temperature_two_time_steps_back                        | air temperature two time steps back                      | K       |    2 | real      | kind_phys | inout  | F        |
!! | qp             | water_vapor_specific_humidity_two_time_steps_back          | water vapor specific humidity two time steps back        | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | psp            | surface_air_pressure_two_time_steps_back                   | surface air pressure two time steps back                 | Pa      |    1 | real      | kind_phys | inout  | F        |
!! | tp1            | air_temperature_at_previous_time_step                      | air temperature at previous time step                    | K       |    2 | real      | kind_phys | inout  | F        |
!! | qp1            | water_vapor_specific_humidity_at_previous_time_step        | water vapor specific humidity at previous time step      | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | psp1           | surface_air_pressure_at_previous_time_step                 | surface air surface pressure at previous time step       | Pa      |    1 | real      | kind_phys | inout  | F        |
!! | u              | critical_relative_humidity                                 | critical relative humidity                               | frac    |    2 | real      | kind_phys | in     | F        |
!! | lprnt          | flag_print                                                 | flag for printing diagnostics to output                  | flag    |    0 | logical   |           | in     | F        |
!! | ipr            | horizontal_index_of_printed_column                         | horizontal index of printed column                       | index   |    0 | integer   |           | in     | F        |
!! | errmsg         | error_message                                              | error message for error handling in CCPP                 | none    |    0 | character | len=*     | out    | F        |
!! | errflg         | error_flag                                                 | error flag for error handling in CCPP                    | flag    |    0 | integer   |           | out    | F        |
!!
!> \section general_gscond GFS gscond Scheme General Algorithm
!! -# Calculate ice-water identification number \f$IW\f$ in order to make a distinction betwee
!! cloud water and cloud ice (table2 of \cite zhao_and_carr_1997).
!! -# Calculate the changes in \f$t\f$, \f$q\f$ and \f$p\f$ due to all the processes except microphysics.
!! -# Calculate cloud evaporation rate (\f$E_c\f$, eq. 19 of \cite zhao_and_carr_1997)
!! -# Calculate cloud condensation rate (\f$C_g\f$, eq.8 of \cite zhao_and_carr_1997) 
!! -# update t,q,cwm due to cloud evaporation and condensation process
!> \section Zhao-Carr_cond_detailed GFS gscond Scheme Detailed Algorithm
!> @{
        subroutine zhaocarr_gscond_run (im,ix,km,dt,dtf,prsl,ps,q,clw1  &
     &,                  clw2, cwm, t, tp, qp, psp                      &
     &,                  tp1, qp1, psp1, u, lprnt, ipr, errmsg, errflg)

!
!     ******************************************************************
!     *                                                                *
!     *  subroutine for grid-scale condensation & evaporation          *
!     *  for the mrf model at ncep.                                    *
!     *                                                                *
!     ******************************************************************
!     *                                                                *
!     *  created by:   q.  zhao         jan. 1995                      *
!     *  modified by:  h.-l. pan        sep. 1998                      *
!     *  modified by:  s. moorthi       aug. 1998, 1999, 2000          *
!     *                                                                *
!     *  references:                                                   *
!     *                                                                *
!     ******************************************************************
!
      use machine , only : kind_phys
      use funcphys , only : fpvs
      use physcons, psat => con_psat, hvap => con_hvap, grav => con_g
     &,             hfus => con_hfus, ttp => con_ttp, rd => con_rd
     &,             cp => con_cp, eps => con_eps, epsm1 => con_epsm1
     &,             rv => con_rv
!      use namelist_def, only: nsdfi,fhdfi
      implicit none
!
! Interface variables
      integer,              intent(in)    :: im, ix, km, ipr
      real(kind=kind_phys), intent(in)    :: dt, dtf
      real(kind=kind_phys), intent(in)    :: prsl(ix,km), ps(im)
      real(kind=kind_phys), intent(inout) :: q(ix,km)
      real(kind=kind_phys), intent(in)    :: clw1(ix,km), clw2(ix,km)
      real(kind=kind_phys), intent(out)   :: cwm(ix,km)
      real(kind=kind_phys), intent(inout) :: t(ix,km)                   &
     &,                     tp(ix,km),   qp(ix,km),   psp(im)           &
     &,                     tp1(ix,km),  qp1(ix,km),  psp1(im)
      real(kind=kind_phys), intent(in)    :: u(im,km)
      logical,              intent(in)    :: lprnt
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
! Local variables
      real (kind=kind_phys) h1
     &,                     d00,  elwv, eliv
     &,                     epsq
     &,                     r,     cpr,  rcp
      parameter (h1=1.e0,       d00=0.e0
     &,          elwv=hvap,     eliv=hvap+hfus
     &,          epsq=2.e-12,   r=rd
     &,          cpr=cp*r,      rcp=h1/cp)
!
      real(kind=kind_phys), parameter :: cons_0=0.0, cons_m15=-15.0
!
      real (kind=kind_phys)  qi(im), qint(im), ccrik, e0
     &,                      cond,   rdt, us, cclimit, climit
     &,                      tmt0, tmt15, qik, cwmik
     &,                      ai, qw, u00ik, tik, pres, pp0, fi
     &,                      at, aq, ap, fiw, elv, qc, rqik
     &,                      rqikk, tx1, tx2, tx3, es, qs
     &,                      tsq, delq, condi, cone0, us00, ccrik1
     &,                      aa, ab, ac, ad, ae, af, ag
     &,                      el2orc, albycp
!     real (kind=kind_phys) vprs(im)
      integer iw(im,km), i, k, iwik
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!-----------------GFS interstitial in driver ----------------------------
       do i = 1,im
         do k= 1,km
              cwm(i,k) = clw1(i,k)+clw2(i,k)
         enddo
       enddo
!-----------------prepare constants for later uses-----------------
!
      el2orc = hvap*hvap / (rv*cp)
      albycp = hvap / cp
!     write(0,*)' in gscond im=',im,' ix=',ix
!
      rdt     = h1/dt
      us      = h1
      cclimit = 1.0e-3
      climit  = 1.0e-20
!
      do  i = 1, im
        iw(i,km) = d00
      enddo
!
!  check for first time step
!
!      if (tp(1,1) < 1.) then
!        do k = 1, km
!          do i = 1, im
!            tp(i,k) = t(i,k)
!            qp(i,k) = max(q(i,k),epsq)
!            tp1(i,k) = t(i,k)
!            qp1(i,k) = max(q(i,k),epsq)
!          enddo
!        enddo
!        do i = 1, im
!          psp(i)  = ps(i)
!          psp1(i) = ps(i)
!        enddo
!      endif
!
!*************************************************************
!> -# Begining of  grid-scale condensation/evaporation loop (start of
!! k-loop, i-loop)
!*************************************************************
!
!     do k = km-1,2,-1
      do k = km,1,-1
!       vprs(:) = 0.001 * fpvs(t(:,k))       ! fpvs in pa
!-----------------------------------------------------------------------
!------------------qw, qi and qint--------------------------------------
        do i = 1, im
          tmt0  = t(i,k)-273.16
          tmt15 = min(tmt0,cons_m15)
          qik   = max(q(i,k),epsq)
          cwmik = max(cwm(i,k),climit)
!
!         ai    = 0.008855
!         bi    = 1.0
!         if (tmt0 .lt. -20.0) then
!           ai = 0.007225
!           bi = 0.9674
!         end if
!
!  the global qsat computation is done in pa
          pres    = prsl(i,k)
!
!         qw      = vprs(i)
          qw      = min(pres, fpvs(t(i,k)))
!
          qw      = eps * qw / (pres + epsm1 * qw)
          qw      = max(qw,epsq)
!         qi(i)   = qw *(bi+ai*min(tmt0,cons_0))
!         qint(i) = qw *(1.-0.00032*tmt15*(tmt15+15.))
          qi(i)   = qw
          qint(i) = qw
!         if (tmt0 .le. -40.) qint(i) = qi(i)

!> -# Compute ice-water identification number IW.
!!\n  The distinction between cloud water and cloud ice is made by the
!! cloud identification number IW, which is zero for cloud water and
!! unity for cloud ice (Table 2 in 
!! \cite zhao_and_carr_1997):
!!  - All clouds are defined to consist of liquid water below the
!! freezing level (\f$T\geq 0^oC\f$) and of ice particles above the
!! \f$T=-15^oC\f$ level.
!!  - In the temperature region between \f$-15^oC\f$ and \f$0^oC\f$,
!! clouds may be composed of liquid water or ice. If there are cloud
!! ice particles above this point at the previous or current time step,
!! or if the cloud at this point at the previous time step consists of
!! ice particles, then the cloud substance at this point is considered
!! to be ice particles because of the cloud seeding effect and the
!! memory of its content. Otherwise, all clouds in this region are
!! considered to contain supercooled cloud water.

!-------------------ice-water id number iw------------------------------
          if(tmt0.lt.-15.0) then
            u00ik = u(i,k)
            fi    = qik - u00ik*qi(i)
            if(fi > d00.or.cwmik > climit) then
               iw(i,k) = 1
            else
              iw(i,k) = 0
            end if
          end if
!
          if(tmt0.ge.0.0) then
            iw(i,k) = 0
          end if
!
          if (tmt0 < 0.0 .and. tmt0 >= -15.0) then
            iw(i,k) = 0
            if (k < km) then
            if (iw(i,k+1)  == 1 .and. cwmik > climit) iw(i,k) = 1
            endif
          end if
        enddo
!> -# Condensation and evaporation of cloud
!--------------condensation and evaporation of cloud--------------------
        do i = 1, im
!>  - Compute the changes in t, q and p (\f$A_{t}\f$,\f$A_{q}\f$ and
!! \f$A_{p}\f$) caused by all the processes except grid-scale
!! condensation and evaporation.
!!\f[
!!   A_{t}=(t-tp)/dt
!!\f]
!!\f[
!!   A_{q}=(q-qp)/dt
!!\f]
!!\f[
!!   A_{p}=(prsl-\frac{prsl}{ps} \times psp)/dt
!!\f]
!------------------------at, aq and dp/dt-------------------------------
          qik   = max(q(i,k),epsq)
          cwmik = max(cwm(i,k),climit)
          iwik  = iw(i,k)
          u00ik = u(i,k)
          tik   = t(i,k)
          pres  = prsl(i,k)
          pp0   = (pres / ps(i)) * psp(i)
          at    = (tik-tp(i,k)) * rdt
          aq    = (qik-qp(i,k)) * rdt
          ap    = (pres-pp0)    * rdt
!>  - Calculate the saturation specific humidity \f$q_{s}\f$ and the
!! relative humidity \f$f\f$ using IW.
!----------------the satuation specific humidity------------------------
          fiw   = float(iwik)
          elv   = (h1-fiw)*elwv    + fiw*eliv
          qc    = (h1-fiw)*qint(i) + fiw*qi(i)
!     if (lprnt) print *,' qc=',qc,' qint=',qint(i),' qi=',qi(i)
!----------------the relative humidity----------------------------------
          if(qc.le.1.0e-10) then
            rqik=d00
          else
            rqik = qik/qc
          endif

!>  - According to \cite sundqvist_et_al_1989,
!! estimate cloud fraction \f$b\f$ at a grid point from relative
!! humidity \f$f\f$ using the equation
!!\f[
!!       b=1-\left ( \frac{f_{s}-f}{f_{s}-u} \right )^{1/2}
!!\f]
!! for \f$f>u\f$; and \f$b=0\f$ for \f$f<u\f$. where \f$f_{s}=1.0\f$ is
!! the relative humidity in a cloud region and \f$u\f$ ,which is an
!! input parameter accounts for the effects of subgrid-scale variations
!! in moisture on large-scale condensation. Since both temperature and
!! moisture may vary at scales smaller than the model grid scale, it is
!! possible for condensation to occur before the grid-average relative
!! humidity reaches 100%. Therefore \f$u\f$ needs to be less than 1.0
!! to account for the subgrid-scale variation of temperature and
!! moisture fields and allow subgrid-scale condensation.
!!  - If cloud fraction \f$b\leq 1.0\times10^{-3}\f$, then evaporate
!! any existing cloud condensate using evaporation rate \f$E_{c}\f$ as
!! computed below.
!!\n If \f$q_{0}\f$ represents the specific humidity at relative
!! humidity \f$u\f$, then
!!\f[
!!           q_{0}=uq_{s}
!!\f]
!!\n if the cloud water/ice at this point is enough to be evaporated
!! until \f$u\f$ is reached, then the evaporation rate \f$E_{c}\f$,
!! assuming that the evaporation process occurs in one time step, is
!! determined by
!!\f[
!!           E_{c}=\frac{q_{0}-q}{dt}
!!\f]
!!\n  Using \f$q_{0}=uq_{s}\f$ and the equation \f$q=fq_{s}\f$,
!! \f$E_{c}\f$ then becomes
!!\f[
!!  E_{c}=\frac{q_{s}}{dt}(u-f)
!!\f]
!! where \f$dt\f$ is the time step for precipitation calculation in the
!! model. It is a simplified version of a higher-order cloud
!! evaporation algorithm (
!! \cite rutledge_and_hobbs_1983). In the case where all clouds will
!! evaporate before \f$u\f$ is reached, the following equation is used:
!! \f[
!!  E_{c}=\frac{cwm}{dt}
!! \f]
!!  - If cloud fraction \f$b>1.0\times10^{-3}\f$, condense water vapor
!! into cloud condensate (\f$C_{g}\f$).
!!\n Using \f$q=fq_{s}\f$, \f$q_{s}=\epsilon e_{s}/p\f$, and the
!! Clausius-Clapeyron equation \f$de_{s}/dT=\epsilon Le_{s}/RT^{2}\f$,
!! where \f$q_{s}\f$ is the saturation specific humidity,\f$e_{s}\f$
!! is the saturation vapor pressure, \f$R\f$ is the specific gas
!! constant for dry air, \f$f\f$ is the relative humidity, and
!! \f$\epsilon=0.622\f$, the expression for \f$C_{g}\f$ has the form
!!\f[
!!  C_{g}=\frac{M-q_{s}f_{t}}{1+(f\epsilon L^{2}q_{s}/RC_{p}T^{2})}+E_{c}
!!\f]
!! where
!!\f[
!!   M=A_{q}-\frac{f\epsilon Lq_{s}}{RT^{2}}A_{t}+\frac{fq_{s}}{p}A_{p}
!!\f]
!! To close the system, an equation for the relative humidity tendency
!! \f$f_{t}\f$ was derived by 
!! \cite sundqvist_et_al_1989 using the hypothesis that the quantity
!! \f$M+E_{c}\f$ is divided into one part,\f$bM\f$,which condenses
!! in the already cloudy portion of a grid square, and another part,
!! \f$(1-b)M+E_{c}\f$,which is used to increase the relative humidity
!! of the cloud-free portion and the cloudiness in the square. The
!! equation is written as
!!\f[
!!  f_{t}=\frac{2(1-b)(f_{s}-u)[(1-b)M+E_{c}]}{2q_{s}(1-b)(f_{s}-u)+cwm/b}
!!\f]
!!  - Check and correct if over condensation occurs.
!!  - Update  t, q and cwm (according to Eqs(6) and (7) in \cite zhao_and_carr_1997)
!!\f[
!!   cwm=cwm+(C_{g}-E_{c})\times dt
!!\f]
!!\f[
!!   q=q-(C_{g}-E_{c})\times dt
!!\f]
!!\f[
!!   t=t+\frac{L}{C_{p}}(C_{g}-E_{c})\times dt
!!\f]
!!\n where \f$L\f$ is the latent heat of condensation/deposition, and
!! \f$C_{p}\f$ is the specific heat of air at constant pressure.

!----------------cloud cover ratio ccrik--------------------------------
          if (rqik .lt. u00ik) then
             ccrik = d00
          elseif(rqik.ge.us) then
             ccrik = us
          else
             rqikk  = min(us,rqik)
             ccrik = h1-sqrt((us-rqikk)/(us-u00ik))
          endif
!-----------correct ccr if it is too small in large cwm regions--------
!         if(ccrik.ge.0.01.and.ccrik.le.0.2.and
!    &          .cwmik.ge.0.2e-3) then
!          ccrik=min(1.0,cwmik*1.0e3)
!         end if
!----------------------------------------------------------------------
!   if no cloud exists then evaporate any existing cloud condensate
!----------------evaporation of cloud water-----------------------------
          e0 = d00
          if (ccrik <= cclimit.and. cwmik > climit)  then
!
!   first iteration - increment halved
!
            tx1 = tik
            tx3 = qik
!
            es   = min(pres, fpvs(tx1))
            qs   = u00ik * eps * es / (pres + epsm1*es)
            tsq  = tx1 * tx1
            delq = 0.5 * (qs - tx3) * tsq / (tsq + el2orc * qs)
!
            tx2   = delq
            tx1   = tx1 - delq * albycp
            tx3   = tx3 + delq
!
!   second iteration
!
            es   = min(pres, fpvs(tx1))
            qs   = u00ik * eps * es / (pres + epsm1*es)
            tsq  = tx1 * tx1
            delq = (qs - tx3) * tsq / (tsq + el2orc * qs)
!
            tx2  = tx2 + delq
            tx1  = tx1 - delq * albycp
            tx3  = tx3 + delq
!
!   third iteration
!
            es   = min(pres, fpvs(tx1))
            qs   = u00ik * eps * es / (pres + epsm1*es)
            tsq  = tx1 * tx1
            delq = (qs - tx3) * tsq / (tsq + el2orc * qs)
            tx2  = tx2 + delq
!
            e0   = max(tx2*rdt, cons_0)
!     if (lprnt .and. i .eq. ipr .and. k .eq. 34)
!    & print *,' tx2=',tx2,' qc=',qc,' u00ik=',u00ik,' rqik=',rqik
!    &,' cwmik=',cwmik,' e0',e0

!           e0 = max(qc*(u00ik-rqik)*rdt, cons_0)
            e0 = min(cwmik*rdt,   e0)
            e0 = max(cons_0,e0)
          end if
!   if cloud cover > 0.2 condense water vapor in to cloud condensate
!-----------the eqs. for cond. has been reorganized to reduce cpu------
          cond = d00
!         if (ccrik .gt. 0.20 .and. qc .gt. epsq) then
          if (ccrik .gt. cclimit .and. qc .gt. epsq) then
             us00   = us  - u00ik
             ccrik1 = 1.0 - ccrik
             aa     = eps*elv*pres*qik
             ab     = ccrik*ccrik1*qc*us00
             ac     = ab + 0.5*cwmik
             ad     = ab * ccrik1
             ae     = cpr*tik*tik
             af     = ae * pres
             ag     = aa * elv
             ai     = cp * aa
             cond   = (ac-ad)*(af*aq-ai*at+ae*qik*ap)/(ac*(af+ag))
!-----------check & correct if over condensation occurs-----------------
             condi  = (qik   -u00ik   *qc*1.0)*rdt
             cond   = min(cond, condi)
!----------check & correct if supersatuation is too high----------------
!             qtemp=qik-max(0.,(cond-e0))*dt
!             if(qc.le.1.0e-10) then
!               rqtmp=0.0
!             else
!               rqtmp=qtemp/qc
!             end if
!             if(rqtmp.ge.1.10) then
!               cond=(qik-1.10*qc)*rdt
!             end if
!-----------------------------------------------------------------------
             cond = max(cond, d00)
!-------------------update of t, q and cwm------------------------------
          end if
          cone0    = (cond-e0) * dt
          cwm(i,k) = cwm(i,k) + cone0
!     if (lprnt .and. i .eq. ipr) print *,' t=',t(i,k),' cone0',cone0
!    &,' cond=',cond,' e0=',e0,' elv=',elv,' rcp=',rcp,' k=',k
!    &,' cwm=',cwm(i,k)
          t(i,k)   = t(i,k)   + elv*rcp*cone0
          q(i,k)   = q(i,k)   - cone0
        enddo                                  ! end of i-loop!
      enddo                                    ! end of k-loop!
!
!*********************************************************************
!> -# End of the condensation/evaporation loop (end of i-loop,k-loop).
!*********************************************************************
!

!> -# Store \f$t\f$, \f$q\f$, \f$ps\f$ for next time step.

      if (dt > dtf+0.001) then     ! three time level
        do k = 1, km
          do i = 1, im
            tp(i,k)  = tp1(i,k)
            qp(i,k)  = qp1(i,k)
!
            tp1(i,k) = t(i,k)
            qp1(i,k) = max(q(i,k),epsq)
          enddo
        enddo
        do i = 1, im
          psp(i)  = psp1(i)
          psp1(i) = ps(i)
        enddo
      else                   ! two time level scheme - tp1, qp1, psp1 not used
        do k = 1, km
!     write(0,*)' in gscond k=',k,' im=',im,' km=',km
          do i = 1, im
!     write(0,*)' in gscond i=',i
            tp(i,k)  = t(i,k)
            qp(i,k)  = max(q(i,k),epsq)
!           qp(i,k)  = q(i,k)
            tp1(i,k) = tp(i,k)
            qp1(i,k) = qp(i,k)
          enddo
        enddo
        do i = 1, im
          psp(i)  = ps(i)
          psp1(i) = ps(i)
        enddo
      endif
!-----------------------------------------------------------------------
      return
      end subroutine zhaocarr_gscond_run
!> @}
! @}


! @}

      end module zhaocarr_gscond
