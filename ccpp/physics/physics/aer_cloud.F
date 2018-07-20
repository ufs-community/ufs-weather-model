       MODULE aer_cloud

#ifdef GEOS5
       use MAPL_ConstantsMod, r8 => MAPL_R8
#endif
#ifdef NEMS_GSM
       use physcons, only: MAPL_PI=>con_pi
       use machine, only : r8 => kind_phys
#endif

! according to the models of Nenes & Seinfeld (2003), Fountoukis and Nenes (2005) and Barahona and Nenes (2008, 2009).
! *** Code Developer: Donifan Barahona donifan.o.barahona@nasa.gov
!
!=======================================================================
!

       implicit none
       private

       public :: aerosol_activate
       public :: AerConversion
       public :: AerConversion1
       public :: AerProps
       public :: getINsubset
       public :: init_Aer
       public :: aer_cloud_init
       public :: vertical_vel_variance


       integer, parameter :: nsmx_par=20, npgauss=10



       type :: AerProps
!      sequence
       real, dimension(nsmx_par) :: num, dpg, sig, den, kap
     &,                             fdust, fsoot, forg
       integer :: nmods
       end type AerProps

       interface assignment (=)
       module procedure copy_aer
       end interface


!==================================================================


!
!=================================================================

!
       real*8, dimension(npgauss) :: xgs_par, wgs_par
!

!Global aux variables

       type(AerProps) :: AerPr_base_clean, AerPr_base_polluted

       real*8         :: base_mass_so4_polluted, base_mass_so4_clean,
     &                   base_mass_ss, frac_dust(5), frac_bc, frac_org,
     &                   ahet_dust(5), ahet_bc

!==================================================================

!
!=================================================================

       real*8 ::  sh_ice, doin_ice,vmin_ice, acorr_dust, acorr_bc
     &,           denw_par, cpair_par, dhv_par

       integer :: typeofspec_ice
       logical :: purehet_ice, purehom_ice,  is_gocart

!==================================================================
!
!==========================

       integer, parameter :: maxit_par=100

       real, parameter :: amw_par=18d-3,      ama_par=29d-3
     &,                   grav_par=9.81d0,    rgas_par=8.31d0
     &,                   accom_par=1.0d0,    eps_par=1d-6
     &,                   zero_par=1.0e-20,   great_par=1d20
     &,                   pi_par=3.1415927d0, sq2pi_par=sqrt(pi_par)
     &,                   sq2_par=1.41421356237d0
!
     &,                   wmw_ice=018d0,      amw_ice=0.029d0
     &,                   rgas_ice=8.314d0,   grav_ice=9.81d0
     &,                   cpa_ice=1005.1d0,   pi_ice=pi_par
     &,                   depcoef_ice=0.1d0,  thaccom_ice=0.7d0
!
     &,                   To_ice=272.15d0,    Tmin_ice=185.d0
     &,                   Pmin_ice=100.0d0,   Thom=236.0d0
     &,                   rv_ice=rgas_ice/wmw_ice


       CONTAINS

       subroutine aer_cloud_init()

       real*8  :: daux, sigaux
       integer ::ix

       call AerConversion_base

       acorr_dust = 2.7e7
       acorr_bc   = 8.0e7

       do ix = 1, 5
         daux          = AerPr_base_polluted%dpg(ix)
         sigaux        = AerPr_base_polluted%sig(ix)
         frac_dust(ix) = 0.5d0*(1d0
     &                      - erfapp(log(0.1e-6/daux)/(sigaux*sq2_par)))

         ahet_dust(ix) = daux*daux*daux*0.52*acorr_dust
     &                 * exp(4.5*sigaux*sigaux)

       end do


       daux     = AerPr_base_polluted%dpg(12)
       sigaux   = AerPr_base_polluted%sig(12)
       frac_bc  = 0.5d0*(1d0-erfapp(log(0.1e-6/daux)/(sigaux*sq2_par)))
       ahet_bc  = daux*daux*daux*0.52*acorr_bc* exp(4.5*sigaux*sigaux)

       daux     = AerPr_base_polluted%dpg(13)
       sigaux   = AerPr_base_polluted%sig(13)
       frac_org = 0.5d0*(1d0-erfapp(log(0.1e-6/daux)/(sigaux*sq2_par)))

       end subroutine aer_cloud_init



! SUBROUTNE AEROSOL ACTIVATE: sets the variables needed for
!the activation  subroutines and  return the activated droplet and ice number concentration

! ===============Input=============:
! tparc_in = T (K)
! pparc_in = P (pa)
! sigwparc_in = variance of the distribution of updraft velocity (m s-1)
! wparc_ls = mean of the distribution of updraft velocity (m s-1)
! Aer_Props = AerProps  structure containing the aerosol properties Aerosol number concentration (Kg-1)
! npre_in =  number concentration of prexisting ice crystals (#/Kg)
! dpre_in =  mass-weighted diameter of prexisting ice crystals (m)
! ccn_diagr8 = array of supersaturations for CCN diagnostics (in-out)

! Ndropr8 = Current droplet number concentration (Kg -1)
! qc = Liquid water mixing ratio (Kg/Kg)
! use_average_v =  .false. integrate over the updraft distribution. True: use the mean vertical velocity
! CCN_param = CCN activation parameterization. 1- Fountoukis and Nenes (2005), 2-Abdul_Razzak and Ghan (2002) (def = 2)
! IN_param =  IN activation spectrum (default is 5)
! ===============Output=============:

! cdncr8 = Activated cloud droplet number concentration (Kg-1)
! smaxliqr8 = Maximum supersaturation w.r.t liquid during droplet activation
! incr8  = Nucleated ice crystal concentration (Kg-1)
! smaxicer8 = Maximum supersaturation w.r.t. ice during ice nucleation
! nheticer8 = Nucleated ice crystal concentration by het freezing (Kg-1)
! INimmr8 =  Nucleated nc by droplet immersion freezing in mixed-phase clouds (Kg-1)
! dINimmr8 = Ice crystal number tendency by immersion freezing (Kg-1 s-1)
! Ncdepr8 = Nucleated nc by deposition ice nucleation (Kg-1)
! Ncdhfr8 = Nucleated nc by immersion in aerosol (Kg -1)
! sc_icer8 = Critical saturation ratio in cirrus
! fdust_depr8 = Fraction of deposition ice nuclei that are dust
! fdust_immr8 = Fraction of immersion mixed-phase ice nuclei that are dust
! fdust_dhfr8 = Fraction of immersion ice nuclei that are dust (not mixed-phase)
! nlimr8 =  Limiting ice nuclei concentration (m-3)

!===================================================================================


       subroutine aerosol_activate(tparc_in, pparc_in, sigwparc_in,
     & wparc_ls, Aer_Props, npre_in, dpre_in, ccn_diagr8, Ndropr8, 
     & cdncr8, smaxliqr8, incr8, smaxicer8, nheticer8, INimmr8,
     & dINimmr8, Ncdepr8, Ncdhfr8, sc_icer8, fdust_immr8, fdust_depr8,
     & fdust_dhfr8, nlimr8, use_average_v, CCN_param, IN_param, fd_dust,
     & fd_soot, pfrz_inc_r8, sigma_nuc, rhi_cell,nccn)




       type(AerProps), intent(in) :: Aer_Props

       logical :: use_average_v

       real(r8), intent(in) :: tparc_in, pparc_in, sigwparc_in,
     & wparc_ls, npre_in, dpre_in, Ndropr8,  fd_soot, fd_dust,
     & sigma_nuc, rhi_cell
       integer, intent(in) :: CCN_param, IN_param, nccn

       real(r8), dimension(:), intent(inout) :: ccn_diagr8

       real(r8), intent(out) :: cdncr8, smaxliqr8, incr8, smaxicer8,
     & nheticer8, INimmr8, dINimmr8, Ncdepr8, Ncdhfr8, sc_icer8,
     & fdust_immr8, fdust_depr8, fdust_dhfr8, nlimr8, pfrz_inc_r8

       type(AerProps) :: Aeraux


       integer :: k, n, I, J, naux


       real*8 :: nact, wparc, tparc,pparc, accom,sigw, smax, antot,
     &           ccn_at_s, sigwparc
!      real*8, allocatable, dimension(:) :: smax_diag
       real*8, dimension(nccn) :: smax_diag


       real*8 :: nhet, nice, smaxice, nlim, air_den, frac, norg, nbc,
     &           nhom, dorg, dbc, kappa, INimm, dINimm, aux

!Anning Cheng move allocable array here for thread safety, 5/4/2016
!      real*8, allocatable, dimension(:) :: sg_par, tp_par, dpg_par,
       real*8, dimension(nsmx_par)  :: sg_par,    tp_par,    dpg_par
     &,                                sig_par,   vhf_par,   ams_par
     &,                                dens_par,  deni_par,  amfs_par
     &,                                kappa_par, ndust_ice, sigdust_ice
     &,                                ddust_ice
!      real*8, allocatable, dimension(:) :: ndust_ice, sigdust_ice,
!    & ddust_ice
       real*8 :: temp_par,  pres_par
       real*8 :: akoh_par,  alfa_par,  bet2_par
       real*8 aka_par, dv_par, psat_par, dair_par,surt_par,ddry_ice, 
     & np_ice,nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,
     & g1_ice, g2_ice,gdoin_ice, z_ice, norg_ice, sigorg_ice, 
     & dorg_ice, dbc_ice,sigbc_ice,lambda_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,sc_ice,Nhet_dep,fdust_imm,fdust_dep,fdust_dhf,
     & waux_ice,fdrop_dust,fdrop_bc,D_preex, N_preex,
     & one_over_tao_preex,P_ice, T_ice,miuv_ice,Nhet_dhf,denice_ice,
     & vpresw_ice,vpresi_ice,denair_ice
       real*8  :: ntot
       integer :: nmodes,act_param,nbindust_ice
       logical :: use_av_v

!=============inputs================
       tparc     = tparc_in
       pparc     = pparc_in
       sigwparc  = sigwparc_in
       miuv_ice  = wparc_ls
       air_den   = pparc*28.8d-3/rgas_par/tparc
       N_preex   = max(npre_in*air_den, zero_par)
       D_preex   = max(dpre_in, 1.0e-9)
       use_av_v  = use_average_v
       act_param = 2
       typeofspec_ice = 5



       smaxicer8   = zero_par
       smaxice     = zero_par
       cdncr8      = zero_par
       smaxliqr8   = zero_par
       incr8       = zero_par
       smaxice     = max(2.349d0-(tparc/259d0) -1.0 , 0.0)
       nheticer8   = zero_par
       nlimr8      = zero_par
       sc_ice      = max(2.349d0-(tparc/259d0), 1.0)
       If (tparc > Thom) sc_ice =1.0

       INimmr8     = zero_par
       dINimmr8    = zero_par
       Ncdepr8     = zero_par
       Ncdhfr8     = zero_par
       fdust_immr8 = zero_par
       fdust_dhfr8 = zero_par
       fdust_depr8 = zero_par
       fdust_imm   = zero_par
       fdust_dhf   = zero_par
       fdust_dep   = zero_par
       pfrz_inc_r8 = zero_par

       nact     = zero_par
       smax     = zero_par
       sc_icer8 = sc_ice

       is_gocart = .true.

       if (sum(Aer_Props%num) <= 1.0e2) then
         return
       end if

       nmodes = max(Aer_Props%nmods, 1)

!      allocate (dpg_par(nmodes))
!      allocate (vhf_par(nmodes))
!      allocate (ams_par(nmodes))
!      allocate (dens_par(nmodes))
!      allocate (sig_par(nmodes))
!      allocate (tp_par(nmodes))
!      allocate (amfs_par(nmodes))
!      allocate (deni_par(nmodes))
!      allocate (sg_par(nmodes))
!      allocate (smax_diag(size(ccn_diagr8)))

!      allocate (kappa_par(nmodes))


       smax_diag = 0.01
       sigw      = zero_par
       do n=1,nmodes
         dpg_par(n)   = zero_par
         vhf_par(n)   = zero_par
         ams_par(n)   = zero_par
         dens_par(n)  = zero_par
         sig_par(n)   = 1d0
         tp_par(n)    = zero_par
         amfs_par(n)  = zero_par
         deni_par(n)  = zero_par
         kappa_par(n) = zero_par
       enddo

       call init_Aer(Aeraux)

       do n=1,nmodes
         tp_par(n)    = DBLE(Aer_Props%num(n))*air_den
         dpg_par(n)   = max(DBLE(Aer_Props%dpg(n)), 1.0e-10)
         sig_par(n)   = DBLE(Aer_Props%sig(n))
         kappa_par(n) = max(DBLE(Aer_Props%kap(n)), 0.001)
         dens_par(n)  = DBLE(Aer_Props%den(n))
         vhf_par(n)   = 3.0
         if (kappa_par(n) > 0.01) then
           ams_par(n) = 18.0e-3*1.7*3.0/kappa_par(n)
         else
           ams_par(n) = 900.0e-3
           tp_par(n)  = 0.0
         endif
         amfs_par(n)  = 1.0
         deni_par(n) = dens_par(n)
       enddo

       kappa_par = max(kappa_par, 0.001)
       dpg_par   = max(dpg_par, 1.0e-10)
       temp_par  = max(tparc, 245.0)
       pres_par  = max(pparc, 34000.0)

       antot = sum(tp_par)
       ntot  = antot

       wparc = max(max(0.8d0*sigwparc, 0.01)+ wparc_ls, 0.01)


       act_param = CCN_param

!============== Calculate cloud droplet number concentration===================

       if (tparc > 245.0) then
         if (antot > 1.0) then

           call ccnspec (tparc,pparc,nmodes
     &,      amfs_par,dens_par,deni_par,vhf_par,ams_par
     &,      sg_par,tp_par,akoh_par,surt_par,temp_par,pres_par
     &,      dv_par,act_param,aka_par, psat_par,dair_par
     &,      ntot,dpg_par)

           if (wparc >= 0.005) then
             if (act_param > 1) then

               call arg_activ (wparc,0.d0,nact,smax,nmodes,tp_par
     &,             dpg_par,kappa_par,sig_par,temp_par, pres_par)

             else

               call pdfactiv (wparc,0.d0,nact,smax,nmodes
     &,             alfa_par,bet2_par,akoh_par,sg_par,tp_par
     &,             temp_par, pres_par,aka_par, dv_par, psat_par
     &,             dair_par,ntot,sig_par)

             endif
           endif

           cdncr8    = max(nact/air_den, zero_par)
           smaxliqr8 = max(smax, zero_par)

!============ Calculate diagnostic CCN number concentration==================

           smax_diag = ccn_diagr8

!          do k =1, size (smax_diag)
           do k =1, nccn
             call ccn_at_super (smax_diag(k), ccn_at_s,nmodes,
     &                          sig_par,sg_par,tp_par)
             ccn_diagr8 (k) = ccn_at_s
           end do

         end if
       end if


! ==========================================================================================
! ==========================================================================================
!==========================  Ice crystal nucleation parameterization  ======================
! ==========================================================================================
       dbc_ice      = 1.0e-9
       nbc_ice      = zero_par
       norg_ice     = zero_par
       dorg_ice     = 1.0e-9
       sigbc_ice    = zero_par
       sigorg_ice   = zero_par
       ddry_ice     = 1.0e-9
       np_ice       = zero_par
       nin_ice      = 0.
       kdust_ice    = 0.
       kbc_ice      = 0.
       shdust_ice   = 0.
       shbc_ice     = 0.
       effdust_ice  = 0.
       effbc_ice    = 0.  
       del1dust_ice = 0.
       si0dust_ice  = 0.
       del1bc_ice   = 0.
       si0bc_ice    = 0.

       naux = 0
       do k = 1, nmodes
         if (kappa_par(k) > 0.1) then
           np_ice   = np_ice   + tp_par(k)
           ddry_ice = ddry_ice + dpg_par(k)
           naux     = naux + 1
         end if
       end do
       ddry_ice = ddry_ice/max(naux , 1)
       frac     = 1.0
       np_ice   = frac*np_ice

!get dust from input structure

       call getINsubset(1, Aer_Props, Aeraux)
       nbindust_ice = max(Aeraux%nmods, 1)

!      allocate(ndust_ice(nbindust_ice))
!      allocate(sigdust_ice(nbindust_ice))
!      allocate(ddust_ice(nbindust_ice))


       do n=1,nbindust_ice
         ddust_ice(n)   = DBLE(Aeraux%dpg(n))
         ndust_ice(n)   = DBLE(Aeraux%num(n))*air_den
         sigdust_ice(n) = DBLE(Aeraux%sig(n))
       enddo


!Black carbon. Only a single mode considered. Use average size and sigma

       call getINsubset(2, Aer_Props, Aeraux)
       naux      = max(Aeraux%nmods, 1)
       dbc_ice   = DBLE(sum(Aeraux%dpg(1:naux)))/naux
       nbc_ice   = DBLE(sum(Aeraux%num(1:naux)))*air_den
       sigbc_ice = DBLE(sum(Aeraux%sig(1:naux)))/naux



       call getINsubset(3, Aer_Props, Aeraux)
       naux       = max(Aeraux%nmods, 1)
       dorg_ice   = DBLE(sum(Aeraux%dpg(1:naux)))/naux
       norg_ice   = DBLE(sum(Aeraux%num(1:naux)))*air_den
       sigorg_ice = DBLE(sum(Aeraux%sig(1:naux)))/naux


       nhet     = zero_par
       nice     = zero_par
       nlim     = zero_par
       INimm    = zero_par
       dINimm   = zero_par
       Nhet_dep = zero_par
       Nhet_dhf = zero_par
       antot=sum(ndust_ice)+ norg_ice+ nbc_ice+ np_ice

       sigwparc=max(0.01, sigwparc)
       sigwparc=min(5.0, sigwparc)
       waux_ice=max(wparc_ls + sigwparc*0.8, 0.01)


!===========Calculate nucleated crystal number. Follows Barahona and Nenes (2008, 2009)==============

       purehet_ice= .FALSE.
       purehom_ice= .FALSE.


      if (antot > 1.0e2) then
        if (tparc < To_ice) then

          CALL prop_ice(tparc, pparc, denice_ice,ddry_ice
     &,    nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice
     &,    g1_ice, g2_ice,gdoin_ice, z_ice,lambda_ice
     &,    kdust_ice, kbc_ice, shdust_ice, shbc_ice
     &,    effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice
     &,    si0bc_ice,nbc_ice,D_preex, N_preex
     &,    one_over_tao_preex,P_ice, T_ice,act_param
     &,    ndust_ice,vpresw_ice,vpresi_ice,denair_ice)


          if (tparc > Thom) then

            fdrop_dust = fd_dust
            fdrop_bc   = fd_soot

            if (sum(ndust_ice)+ norg_ice+ nbc_ice .gt. 1.e3) then

              call INimmersion(INimm, dINimm, waux_ice,dbc_ice,sigbc_ice
     &,            nbc_ice,fdust_imm,fdrop_dust,fdrop_bc,ndust_ice
     &,            sigdust_ice,ddust_ice,nbindust_ice
     &,            vpresw_ice,vpresi_ice,T_ice)


              ndust_ice = max(ndust_ice*(1.0-fdrop_dust), 0.0)
              nbc_ice   = max(nbc_ice*(1.0-fdrop_bc), 0.0)


              call IceParam (sigwparc, denice_ice,ddry_ice,np_ice
     &,         nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice
     &,         g1_ice, g2_ice,gdoin_ice, z_ice,lambda_ice,sc_ice
     &,         norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice
     &,         kdust_ice, kbc_ice, shdust_ice, shbc_ice
     &,         effdust_ice, effbc_ice, del1dust_ice, si0dust_ice
     &,         del1bc_ice,  si0bc_ice,nbc_ice,one_over_tao_preex
     &,         nhet, nice, smaxice, nlim,Nhet_dep,Nhet_dhf,fdust_dep
     &,         P_ice,T_ice,ndust_ice,sigdust_ice,ddust_ice,nbindust_ice
     &,         use_av_v,miuv_ice,vpresw_ice,vpresi_ice,denair_ice)
            end if

            sc_ice = 1.0

          else

            call IceParam (sigwparc, denice_ice,ddry_ice,np_ice
     &,      nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice
     &,      g1_ice, g2_ice,gdoin_ice, z_ice,lambda_ice,sc_ice
     &,      norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice
     &,      kdust_ice, kbc_ice, shdust_ice, shbc_ice
     &,      effdust_ice, effbc_ice, del1dust_ice, si0dust_ice
     &,      del1bc_ice,  si0bc_ice,nbc_ice,one_over_tao_preex
     &,      nhet, nice, smaxice, nlim,Nhet_dep,Nhet_dhf,fdust_dep
     &,      P_ice, T_ice,ndust_ice, sigdust_ice,ddust_ice,nbindust_ice
     &,      use_av_v,miuv_ice,vpresw_ice,vpresi_ice,denair_ice)

          end if

          aux = (sc_ice - rhi_cell)/(sq2_par*sigma_nuc)

          pfrz_inc_r8 = 1.0d0- 0.5d0*(1.0d0+erf(aux))
          pfrz_inc_r8 = min(max(pfrz_inc_r8, 0.0), 0.999)

        end if
      end if


!======================== use sc_ice only for cirrus
       If (tparc > Thom) sc_ice =1.0
!==========================

!All # m-3 except those passed to MG later
       smaxicer8   = min(max(smaxice, zero_par), 2.0)
       nheticer8   = min(max(nhet, zero_par), 1e10)
       incr8       = min(max(nice/air_den, zero_par), 1e10)
       nlimr8      = min(max(nlim, zero_par), 1e10)
       sc_icer8    = min(max(sc_ice, 1.0), 2.0)
       INimmr8     = min(max(INimm, zero_par), 1e10)
       dINimmr8    = min(max(dINimm/air_den, zero_par), 1e10)
       Ncdepr8     = min(max(Nhet_dep, zero_par), 1e10)
       Ncdhfr8     = min(max(Nhet_dhf, zero_par), 1e10)
       fdust_immr8 = min(max(fdust_imm, zero_par), 1e10)
       fdust_depr8 = min(max(fdust_dep, zero_par), 1e10)
       fdust_dhfr8 = min(max(fdust_dhf, zero_par), 1e10)

!      deallocate (ndust_ice)
!      deallocate (sigdust_ice)
!      deallocate (ddust_ice)
!      deallocate (dpg_par)
!      deallocate (vhf_par)
!      deallocate (ams_par)
!      deallocate (dens_par)
!      deallocate (sig_par)
!      deallocate (tp_par)
!      deallocate (amfs_par)
!      deallocate (deni_par)
!      deallocate (sg_par)
!      deallocate (smax_diag)

!      deallocate (kappa_par)



2033    return

       END subroutine aerosol_activate
!


!=======================================================================
!
! *** SUBROUTINE AerConversion_base
! *** This subrotine sets basic properties of the aerosol size distributions when using GOCART aerosol
!****Mass-number conversion based on Barahona at al. GMD, 2014.
!=======================================================================

!Output:



!
       SUBROUTINE AerConversion_base ()

       integer, parameter    :: NMDM = 20
       real, dimension(NMDM) :: TPI, DPGI, SIGI, DENSI, KAPPAS, FDUST,
     &                          FSOOT, FORG, TPI_aux, DPGI_aux, SIGI_aux

       real:: aux
       integer:: nmod, K

       do k=1,NMDM
         TPI(k)    = 0.0
         DPGI(k)   = 1.0e-9
         SIGI(k)   = 2.0
         DENSI(k)  = 2200.0
         KAPPAS(k) = 0.01
         FDUST(k)  = 0.0
         FSOOT(k)  = 0.0
         FORG(k)   = 0.0
       enddo
       nmod = 13

! Gocart aerosol size distributions for dust

!!!!!!!!!!!!!!======================================
!!!!!!!!!   Dust
!!!!!!!!!!!!!!======================================


       do k=1,5
         SIGI(k)   = log(1.8)
         DENSI(k)  = 1700.0
         KAPPAS(k) = 0.0001
         FDUST(k)  = 1.0
       enddo



       DPGI (1) = 1.46e-6

       DPGI (2) = 2.80e-6

       DPGI (3) = 4.80e-6
!! Dust 4: 3-6
       DPGI (4) = 9.0e-6
!!  Dust 5: 6-10
       DPGI (5) = 16.0e-6

       DO K =1 , 5
         TPI(K) = 6.0/(DENSI(K)*pi_par*exp(4.5*SIGI(K)*SIGI(K))*DPGI(K)*
     &                 DPGI(K)*DPGI(K))
       END DO

!!!!!!!!!!!!!!======================================
!!!!!!!!!   Sea Salt (Using 3 modes based on Barahona et al. GMD. 2014.
!!!!!!!!!!!!!!======================================



       DENSI(6:8)  = 2200.0
       KAPPAS(6:8) = 1.28



       TPI  (6) = 100e6
       DPGI (6) = .01e-6
       SIGI (6) = log(1.6)

       TPI  (7) = 60.0e6
       DPGI (7) = 0.071e-6
       SIGI (7) = log(2.0)

       TPI  (8) = 3.1e6
       DPGI (8) = 0.62e-6
       SIGI (8) = log(2.7)

       aux = 0.
       DO K =6 , 8
         aux = (TPI(K)*DENSI(K)*pi_par*exp(4.5*SIGI(K)*SIGI(K))*DPGI(K)*
     &       DPGI(K)*DPGI(K))/6.0 + aux
       END DO
       base_mass_ss = aux




       KAPPAS(9:11) = 0.65
       DENSI(9:11)  = 1650.0

! Different size distributions for polluted and clean environments


       TPI  (9) = 1.06e11
       DPGI (9) = .014e-6
       SIGI (9) = log(1.8d0)

       TPI  (10) = 3.2e10
       DPGI (10) = 0.054e-6
       SIGI (10) = log(2.16)

       TPI  (11) = 5.4e6
       DPGI (11) = 0.86e-6
       SIGI (11) = log(2.21)

       aux = 0.
       DO K =9, 11
         aux = (TPI(K)*DENSI(K)*pi_par*exp(4.5*SIGI(K)*SIGI(K))*DPGI(K)*
     &          DPGI(K)*DPGI(K))/6.0 + aux
       END DO
       base_mass_so4_polluted = aux




       TPI_aux  (9) = 1.0e9
       DPGI_aux (9) = .016e-6
       SIGI_aux (9) = log(1.6d0)

       TPI_aux  (10) = 8.0e8
       DPGI_aux (10) = 0.067e-6
       SIGI_aux (10) = log(2.1)

       TPI_aux  (11) = 2.0e6
       DPGI_aux (11) = 0.93e-6
       SIGI_aux (11) = log(2.2)


       aux = 0.
       DO K =9, 11
         aux =(TPI_aux(K)*DENSI(K)*pi_par*exp(4.5*SIGI_aux(K)*
     &    SIGI_aux(K))*DPGI_aux(K)*DPGI_aux(K)*DPGI_aux(K))/6.0 + aux
       END DO
       base_mass_so4_clean = aux




       DPGI (12)  = 0.0118*2.e-6
       SIGI (12)  = log(2.00)
       DENSI(12)  = 1600.0
       KAPPAS(12) = 0.0001
       FSOOT(12)  = 1.0
       TPI(12)    = 6.0/(DENSI(12)*pi_par*exp(4.5*SIGI(12)*SIGI(12))*
     &                   DPGI(12)*DPGI(12)*DPGI(12))



       DPGI (13)  = 0.0212*2.e-6
       SIGI (13)  = log(2.20)
       DENSI(13)  = 900.0
       KAPPAS(13) = 0.0001
       FORG(13)   = 1.0
       TPI(13)    = 6.0/(DENSI(13)*pi_par*exp(4.5*SIGI(13)*SIGI(13))*
     &                   DPGI(13)*DPGI(13)*DPGI(13))

       call init_Aer(AerPr_base_polluted)
       call init_Aer(AerPr_base_clean)


       do k=1,nmod
         AerPr_base_polluted%num(k)   = TPI(k)
         AerPr_base_polluted%dpg(k)   = DPGI(k)
         AerPr_base_polluted%sig(k)   = SIGI(k)
         AerPr_base_polluted%kap(k)   = KAPPAS(k)
         AerPr_base_polluted%den(k)   = DENSI(k)
         AerPr_base_polluted%fdust(k) = FDUST(k)
         AerPr_base_polluted%fsoot(k) = FSOOT(k)
         AerPr_base_polluted%forg(k)  = FORG(k)
       enddo
       AerPr_base_polluted%nmods = nmod

       AerPr_base_clean           = AerPr_base_polluted
       AerPr_base_clean%num(9:11) = TPI_aux(9:11)
       AerPr_base_clean%dpg(9:11) = DPGI_aux(9:11)
       AerPr_base_clean%sig(9:11) = SIGI_aux(9:11)

       RETURN
!
       END SUBROUTINE AerConversion_base


!=======================================================================
!
! *** SUBROUTINE AerConversion
! *** This subrotine sets the properties of the aerosol distributions
!****Mass-number conversion based on Barahona at al. GMD, 2014.
!=======================================================================
!Input.

!Output:



!
       SUBROUTINE AerConversion (aer_mass, AerPr, kappa, SULFATE, ORG,
     &                           BCARBON, DUST, SEASALT)


       type(AerProps), intent (out), dimension(:,:,:) :: AerPr


       real, dimension(:,:,:,:), intent(in) :: aer_mass
       real, intent (out), dimension(:,:,:) :: kappa, DUST, SULFATE,
     &                                         BCARBON, ORG, SEASALT
       real:: aux, densSO4, densORG, k_SO4, k_ORG, k_SS, tot_mass, dens,
     &        kappa_aux, tx1

       integer :: i,j,k,l
       integer :: im, jm, lm
       type(AerProps) :: AeroAux
       real, dimension(size(aer_mass,4)) :: aer_mass_tmp

       im = size(aer_mass,1)
       jm = size(aer_mass,2)
       lm = size(aer_mass,3)

       call init_Aer(AeroAux)

       do k = 1, lm
         do j = 1, jm
           do i = 1, im
             aer_mass_tmp(:) = aer_mass(i,j,k,:)


             tot_mass  = aer_mass_tmp(11) + aer_mass_tmp (15)
     &                 + aer_mass_tmp (14)
             densSO4   = 1700.0
             densORG   = 1600.0
             k_SO4     = 0.65
             k_ORG     = 0.2
             kappa_aux = 0.65



             if (tot_mass > 2.0e-20) then
               tx1 = 1.0 / tot_mass
               dens = (aer_mass_tmp(11)*densSO4
     &              + sum(aer_mass_tmp(14:15))*densORG) * tx1
               kappa_aux = (aer_mass_tmp(11)*k_SO4
     &                   + sum(aer_mass_tmp(14:15))*k_ORG) * tx1
             else
               dens = 1750.0
               kappa_aux = 0.65
             end if


             if (tot_mass > 5.0e-8) then
               AeroAux = AerPr_base_polluted
               AeroAux%num(9:11) = AeroAux%num(9:11)*tot_mass
     &                           / base_mass_so4_polluted
             else
               AeroAux = AerPr_base_clean
               AeroAux%num(9:11) = AeroAux%num(9:11)*tot_mass
     &                           / base_mass_so4_clean
             end if

             AeroAux%kap(9:11) = max(kappa_aux, 0.001)
             AeroAux%den(9:11) = dens
             SULFATE(i,j,k)    = SUM(AeroAux%num(9:11))
             kappa_aux         = kappa_aux*tot_mass




             AeroAux%num(1:5) = AeroAux%num(1:5) *aer_mass_tmp(1:5)
             kappa_aux        = kappa_aux
     &                        + AeroAux%kap(1)*sum(aer_mass_tmp(1:5))
             DUST(i,j,k)      = sum(AeroAux%num(1:5))


             tot_mass         = sum(aer_mass_tmp(6:10))
             AeroAux%num(6:8) = AeroAux%num(6:8)*tot_mass/base_mass_ss
             kappa_aux        = kappa_aux + AeroAux%kap(6)*tot_mass
             SEASALT(i,j,k )  = sum(AeroAux%num(6:8))


             AeroAux%num(12) = AeroAux%num(12) *aer_mass_tmp(13)
             kappa_aux       = kappa_aux
     &                       + AeroAux%kap(12)*aer_mass_tmp(13)
             BCARBON(i,j,k)  = AeroAux%num(12)


             AeroAux%num(13) = AeroAux%num(13) *aer_mass_tmp(15)
             ORG(i, j, k)    = AeroAux%num(13)
             tot_mass        = sum(aer_mass_tmp)

             if (tot_mass > 0.0) then
               kappa(i,j,k) = kappa_aux/tot_mass
             end if

             AerPr(i,j,k) = AeroAux

           end do
         end do
       end do

       RETURN
!
       END SUBROUTINE AerConversion

!=======================================================================
!=======================================================================
!
       SUBROUTINE AerConversion1 (aer_mass, AerPr)


       type(AerProps),               dimension(:,:) :: AerPr
!      type(AerProps), intent (out), dimension(:,:) :: AerPr


       real, dimension(:,:,:)             :: aer_mass
!      real, dimension(:,:,:), intent(in) :: aer_mass
       real:: aux, densSO4, densORG, k_SO4, k_ORG, k_SS, tot_mass, dens,
     &        kappa_aux, tx1, tx2

       integer :: i,k,l,im,lm
       type(AerProps) :: AeroAux
       real, dimension(size(aer_mass,3)) :: aer_mass_tmp

       im = size(aer_mass,1)
       lm = size(aer_mass,2)

       call init_Aer(AeroAux)

       do k = 1, lm
         do i = 1, im
           aer_mass_tmp(:) = aer_mass(i,k,:)


           tot_mass  = aer_mass_tmp(11) + aer_mass_tmp (15)
     &               + aer_mass_tmp (14)
           densSO4   = 1700.0
           densORG   = 1600.0
           k_SO4     = 0.65
           k_ORG     = 0.2
           kappa_aux = 0.65



           if (tot_mass > 2.0e-20) then
             tx1       = 1.0 / tot_mass
             tx2       = aer_mass_tmp(14) + aer_mass_tmp(15)
             dens      = (aer_mass_tmp(11)*densSO4 + tx2*densORG) * tx1
             kappa_aux = (aer_mass_tmp(11)*k_SO4   + tx2*k_ORG) * tx1
           else
             dens      = 1750.0
             kappa_aux = 0.65
           end if


           if (tot_mass > 5.0e-8) then
             AeroAux = AerPr_base_polluted
             AeroAux%num(9:11) = AeroAux%num(9:11)*tot_mass
     &                         / base_mass_so4_polluted
           else
             AeroAux = AerPr_base_clean
             AeroAux%num(9:11) = AeroAux%num(9:11)*tot_mass
     &                         / base_mass_so4_clean
           end if

           AeroAux%kap(9:11) = max(kappa_aux, 0.001)
           AeroAux%den(9:11) = dens
           kappa_aux         = kappa_aux*tot_mass




           AeroAux%num(1:5) = AeroAux%num(1:5) *aer_mass_tmp(1:5)
           kappa_aux        = kappa_aux
     &                      + AeroAux%kap(1)*sum(aer_mass_tmp(1:5))


           tot_mass         = sum(aer_mass_tmp(6:10))
           AeroAux%num(6:8) = AeroAux%num(6:8)*tot_mass/base_mass_ss
           kappa_aux        = kappa_aux + AeroAux%kap(6)*tot_mass


           AeroAux%num(12) = AeroAux%num(12)*aer_mass_tmp(13)
           kappa_aux       = kappa_aux+AeroAux%kap(12)*aer_mass_tmp(13)


           AeroAux%num(13) = AeroAux%num(13) *aer_mass_tmp(15)
           tot_mass        = sum(aer_mass_tmp)


           AerPr(i,k) = AeroAux

         end do
       end do

       RETURN
!
       END SUBROUTINE AerConversion1

!=======================================================================
!=======================================================================
!=======================================================================

!*********** Calculate subgrid scale distribution of vertical velocity****
! ====================================================================

      subroutine vertical_vel_variance(omeg, lc_turb, tm_gw, pm_gw,
     & rad_cool, uwind_gw, tausurf_gw, nm_gw, LCCIRRUS, Nct, Wct, ksa1,
     & fcn, KH, FRLAND, ZPBL, Z, maxkhpbl, wparc_ls, wparc_gw,
     & wparc_cgw, wparc_turb, LTS)


       real(r8), intent(in) :: omeg, tm_gw, lc_turb, rad_cool, uwind_gw,
     & pm_gw
       real , intent(in) :: LCCIRRUS, KH, ZPBL, Z, FRLAND, nm_gw,
     & tausurf_gw, ksa1, fcn, maxkhpbl, Nct, Wct, LTS

       real(r8), intent(out) :: wparc_ls, wparc_gw, wparc_cgw,
     & wparc_turb

       real(r8) :: rho_gw, k_gw, h_gw, c2_gw, dummyW, maxkh, Wbreak



!!!:========= mean V Large scale and radiative cooling
       rho_gw = pm_gw*28.8d-3/rgas_par/tm_gw


       wparc_ls =-omeg/rho_gw/grav_ice + cpa_ice*rad_cool/grav_ice

!!!======== Orographic Gravity gwave (and brackground) initiated (According to Barahona et al. 2013 GMD)

       wparc_gw = 0.0
       k_gw = 2d0*pi_par/LCCIRRUS

       h_gw= k_gw*rho_gw*uwind_gw*nm_gw

       if (h_gw .gt. 0.0) then
       h_gw=sqrt(2.0*tausurf_gw/h_gw)
       else
       h_gw = 0.0
       end if

       Wbreak = 0.133*k_gw*uwind_gw/nm_gw

       wparc_gw = k_gw*uwind_gw*h_gw*0.133
       wparc_gw = min(wparc_gw, Wbreak)
       wparc_gw = wparc_gw*wparc_gw


!!!======== Subgrid variability from Convective Sources According to Barahona et al. 2014 in prep

       wparc_cgw = 0.0
       c2_gw     = (nm_gw+Nct)/Nct
       wparc_cgw = sqrt(ksa1)*fcn*c2_gw*Wct*0.6334
       wparc_cgw = min(wparc_cgw, Wbreak)
       wparc_cgw = wparc_cgw*wparc_cgw

!!!:=========Subgrid scale variance from turbulence


       wparc_turb=KH/lc_turb


       if (LTS > 20.0) then
         wparc_turb=maxkhpbl/lc_turb
       end if


       wparc_turb= wparc_turb*wparc_turb



      end subroutine vertical_vel_variance
!=======================================================================
!=======================================================================
!=======================================================================

!=======================================================================
!Extracts aerosol props with INactive = typ
      subroutine getINsubset(typ, aerin, aerout)

! typ  : type of aerosol needed: 1 dust, 2 soot, 3 organics
! nbins: number of modes with num>0

       type(AerProps), intent (in)    :: aerin
       type(AerProps), intent (inout) :: aerout
       integer, intent(in)            :: typ

       integer:: k, bin

       call init_Aer(aerout)
       bin = 0

       do k=1, aerin%nmods

         if (typ == 1) then
           if (aerin%fdust(k) > 0.9) then
             bin = bin + 1
             call copy_mode(aerout,aerin, k,bin)
           end if
         elseif (typ == 2) then
           if (aerin%fsoot(k) > 0.9) then
             bin = bin + 1
             call copy_mode(aerout,aerin, k,bin)
           end if
         elseif (typ == 3) then
           if (aerin%forg(k) > 0.9) then
             bin = bin + 1
             call copy_mode(aerout,aerin, k,bin)
           end if
         end if

       end do

       aerout%nmods = max(bin, 1)

      end subroutine getINsubset

!========================subroutines to handle aer strucuture=====================================


       subroutine copy_Aer(a,b)

       type (AerProps), intent(out) :: a
       type (AerProps), intent(in)  :: b

       a%num   = b%num
       a%sig   = b%sig
       a%dpg   = b%dpg
       a%kap   = b%kap
       a%den   = b%den
       a%fdust = b%fdust
       a%fsoot = b%fsoot
       a%forg  = b%forg
       a%nmods = b%nmods

       end subroutine copy_Aer

       subroutine copy_mode(a_out,a_in, mode_in, mode_out)
       type (AerProps), intent(out) :: a_out
       type (AerProps), intent(in) :: a_in
       integer, intent (in) :: mode_in, mode_out

       a_out%num(mode_out)= a_in%num(mode_in)
       a_out%sig(mode_out) = a_in%sig(mode_in)
       a_out%dpg(mode_out) = a_in%dpg(mode_in)
       a_out%kap(mode_out) = a_in%kap(mode_in)
       a_out%den(mode_out) = a_in%den(mode_in)
       a_out%fdust(mode_out) = a_in%fdust(mode_in)
       a_out%fsoot(mode_out) = a_in%fsoot(mode_in)
       a_out%forg(mode_out) = a_in%forg(mode_in)

       end subroutine copy_mode

       subroutine init_Aer(aerout)

       type (AerProps), intent(inout) :: aerout
       integer n

       do n=1,nsmx_par
         aerout%num(n)   =  0.
         aerout%dpg(n)   = 1.0e-9
         aerout%sig(n)   = 2.0
         aerout%kap(n)   = 0.2
         aerout%den(n)   = 2200.0
         aerout%fdust(n) = 0.0
         aerout%fsoot(n) = 0.0
         aerout%forg(n)  = 0.0
       enddo
       aerout%nmods = 1

       end subroutine init_Aer



!!!!!!!!!!!!!!======================================
!!!!!!!!!   Subroutine ARG_act: finds the activated droplet number following Abdul_Razzak and Ghan 2000.
!Written by Donifan Barahona
!!donifan.o.barahona@nasa.gov
!!!!!!!!!!!!!!====================================

       subroutine arg_activ (wparc,sigw,nact,smax,nmodes,tp_par,
     &  dpg_par,kappa_par,sig_par,temp_par, pres_par)


       real*8, intent(in) :: wparc, sigw
       integer, intent(in):: nmodes
       real*8, intent(out) :: nact, smax
       real*8 :: SMI(nmodes),tp_par(nmodes)
       real*8, dimension(nmodes)::dpg_par,kappa_par,sig_par

       real*8 :: alfa, beta, Akoh, G, T, fi, gi, nui, citai, ui, aux1,
     & PACT, auxx, aux,temp_par, pres_par
       integer :: I, INDEX

       T = min(max(temp_par, 243.0), 323.0)
       alfa=2.8915E-08*(T*T) - 2.1328E-05*T + 4.2523E-03
       beta=exp(3.49996E-04*T*T - 2.27938E-01*T + 4.20901E+01)
       G=exp(-2.94362E-06*T*T*T + 2.77941E-03*T*T - 8.92889E-01*T +
     & 1.18787E+02)
       Akoh= 0.66e-6/T


       DO I = 1, nmodes
       aux =0.667*Akoh/dpg_par(I)
       SMI (I) = (aux*sqrt(aux))/SQRT(2.0*kappa_par(I))
       END DO



       SMI=MAX(SMI, 1.0e-5)
       aux =alfa*wparc*G
       aux = aux*sqrt(aux)/(2.d0*pi_par*980.d0*beta)
       citai = 0.667*Akoh*SQRT(alfa*wparc*G)

       auxx=0.0
       DO INDEX =1, nmodes
       if (tp_par(INDEX) .gt. 0.0) then
       fi=0.5*exp(2.5*sig_par(INDEX))
       gi=1.0+0.25*sig_par(INDEX)
       nui=aux/tp_par(INDEX)
       aux1=fi*((citai/nui)*sqrt(citai/nui)) + gi*(SMI(INDEX)*
     &SMI(INDEX) /(nui+(3.0*citai)))**0.75
       aux1=aux1/(SMI(INDEX)*SMI(INDEX))
       auxx=auxx+aux1

       end if
       end do


       if (auxx .gt. 0.0) then
         smax = 1/sqrt(auxx)
       else
         nact = 0.0
         return
       end if


       auxx = 0.0


       DO INDEX = 1, nmodes
         ui   = sq2_par*log(SMI(INDEX)/smax)/3.0
         aux1 = 0.5*tp_par(INDEX)*(1.0-ERFAPP(ui))
         auxx = auxx + aux1
       END DO

       nact = auxx


       End Subroutine arg_activ

!===================================================================================================

!===================================================================================================




!=====================================================================================================
       subroutine ccn_at_super (super,ccn_at_s,nmodes,
     &                          sig_par,sg_par,tp_par)

       integer :: j, I,nmodes
       real*8  :: dlgsg, dlgsp, orism5, ndl, nds, super, ccn_at_s
       real*8, dimension(nmodes) :: sig_par,sg_par,tp_par

       ndl = 0d0

       do j=1, nmodes

         dlgsg = sig_par(j)
 
         if (sg_par(j) .gt. 0.0) then
           if (super .gt. 0.0) then
             dlgsp = dlog(sg_par(j)/super)
           else
             dlgsp = dlog(sg_par(j)/0.01)
           end if
         else
           dlgsp = 0.0

         end if
         orism5 = 2.d0*dlgsp/(3d0*sq2_par*dlgsg)
         ndl    = (tp_par(j)/2.0)*(1.0-erf(orism5))+ndl

       end do

       ccn_at_s = ndl

       end subroutine ccn_at_super

!=======================================================================

!    subroutine ccnspec
!------------------------------
!     DESCRIPTION
!
! *** subroutine ccnspec
! *** this subroutine calculates the ccn spectrum of the aerosol using
!     the appropriate form of kohler theory and assuming a lognormal aerosol size dist
!
! *** written by athanasios nenes
!!
!      Code Developer
!      Donifan Barahona
!      donifanb@umbc.edu
!=======================================================================
!
       subroutine ccnspec (tparc,pparc,nmodes,
     & amfs_par,dens_par,deni_par,vhf_par,ams_par,
     & sg_par,tp_par,akoh_par,surt_par,temp_par,pres_par
     & ,dv_par,act_param,aka_par,  psat_par,dair_par,
     & ntot,dpg_par)
!

       integer :: nmodes,i,j,k
       real*8  :: tparc, pparc, amfi,denp,vlfs,par1, par2
       real*8, dimension(nmodes)::amfs_par,dens_par,deni_par,
     &                   vhf_par,ams_par,sg_par,tp_par,dpg_par
       real*8 akoh_par,surt_par,temp_par,pres_par,
     &        aka_par, dv_par, psat_par,dair_par,ntot
       integer act_param


       ntot=zero_par
       temp_par = max(tparc, 245.0)
       pres_par = max(pparc, 34000.0)
!
       call props(pres_par,temp_par,surt_par,dv_par,act_param,
     & aka_par, psat_par,dair_par)


!
       akoh_par = 4d0*amw_par*surt_par/rgas_par/temp_par/denw_par
!

       do k=1,nmodes


       amfi = max(1.0d0-amfs_par(k),zero_par)
       denp = amfs_par(k)*dens_par(k) + amfi*deni_par(k)
       vlfs = amfs_par(k)/dens_par(k)/(amfs_par(k)/dens_par(k)+ amfi/
     &deni_par(k))
       par1 = 4d0*denw_par*ams_par(k)/27d0/vhf_par(k)/denp/amw_par/
     & dpg_par(k)**3


       par1 = par1/vlfs
       par2 = sqrt(max(par1*akoh_par*akoh_par*akoh_par, zero_par))
       sg_par(k)= max(exp(par2) - 1d0, zero_par)
       ntot=ntot+tp_par(k)
       enddo
!


! *** end of subroutine ccnspec ****************************************
!
       return
       end subroutine ccnspec


!=======================================================================
!    -------------------------------
!     DESCRIPTION
!
! *** subroutine pdfactiv
! *** this subroutine calculates the ccn activation fraction according
!     to the nenes and seinfeld (2003) parameterization, with
!     modification for non-contunuum effects as proposed by fountoukis
!     and nenes (2005). this routine calculates for a pdf of
!     updraft velocities.
!
! *** written by athanasios nenes
!
!      Code Developer
!      Donifan Barahona
!      donifanb@umbc.edu

!=======================================================================
!
       subroutine pdfactiv (wparc,sigw, nact,smax,nmodes,
     &  alfa_par,bet2_par,akoh_par,sg_par,tp_par,
     &  temp_par, pres_par,aka_par, dv_par, psat_par,
     &  dair_par,ntot,sig_par)
!
!
       integer :: i, isec,nmodes
       real*8 :: tpart, nact, nacti,wparc,smax,ntot, tx1
       real*8 :: pdf, dpnmx,sigw,plimt,probi,whi,wlo, scal, wpi,smaxi,
     & alfa_par,bet2_par,akoh_par,sg_par(nmodes),tp_par(nmodes),
     & tparc,pparc,temp_par, pres_par,aka_par, dv_par, psat_par,
     & dair_par,sig_par(nmodes)
       real*8 :: wpdbg(npgauss), pddbg(npgauss), nadbg(npgauss),
     &           smdbg(npgauss)
!
! *** case where updraft is very small
!
       if (wparc <= 1d-6) then
         smax = 0d0
         nact = 0d0
         isec = 1
         dpnmx = great_par
         return
       endif
!
! *** single updraft case
!
       if (sigw < 1e-10) then
         call activate (wparc,nact,smax,nmodes,
     &         alfa_par,bet2_par,akoh_par,sg_par,tp_par,
     &         temp_par, pres_par,aka_par, dv_par, psat_par,
     &         dair_par,ntot,sig_par)
         wpdbg(1) = wparc
         pddbg(1) = 1.0
         nadbg(1) = nact
         smdbg(1) = smax
!
! *** pdf of updrafts
!
       else
         nact  = zero_par
         smax  = zero_par
         plimt = 1e-3
         probi = sqrt(-2.0*log(plimt*sigw*sq2pi_par))
         whi   = wparc + sigw*probi
         wlo   = 0.05
         scal  = 0.5*(whi-wlo)
         do i=1,npgauss
           wpi = wlo + scal*(1.0-xgs_par(i))
           call activate (wpi,nacti,smaxi,nmodes,
     &          alfa_par,bet2_par,akoh_par,sg_par,tp_par,
     &          temp_par, pres_par,aka_par, dv_par, psat_par,
     &          dair_par,ntot,sig_par)
           tx1  = (wpi-wparc)/ sigw
           pdf  = (1.0/sq2pi_par/sigw) * exp(-0.5*tx1*tx1)
           nact = nact + wgs_par(i)*(pdf*nacti)
           smax = smax + wgs_par(i)*(pdf*smaxi)
           wpdbg(i) = wpi
           pddbg(i) = pdf
           nadbg(i) = nacti
           smdbg(i) = smaxi
           if (pdf < plimt) exit
         enddo
         nact = nact*scal
         smax = smax*scal
       endif
!
       return
!
! *** end of subroutine pdfactiv ****************************************
!
       end subroutine pdfactiv



!=======================================================================

!     DESCRIPTION
!
! *** subroutine activate
! *** this subroutine calculates the ccn activation fraction according
!     to the nenes and seinfeld (2003) parameterization, with
!     modification for non-contunuum effects as proposed by fountoukis
!     and nenes (2005).
!

!=======================================================================
!
       subroutine activate (wparc,ndroplet,smax,nmodes,
     &  alfa_par,bet2_par,akoh_par,sg_par,tp_par,temp_par, pres_par,
     &  aka_par, dv_par, psat_par,dair_par,ntot,sig_par)

!
       integer :: i,niter,nmodes

       real*8 :: ndrpl, wparc, beta,cf1, cf2,x1,sinteg1,sinteg2,y1, x2,
     &y2,x3,y3, sign,smax, ent_par, ndroplet, nrdpl,bet1_par,
     &  alfa_par,bet2_par,akoh_par,sg_par(nmodes),tp_par(nmodes),
     &  temp_par, pres_par,aka_par, dv_par, psat_par,dair_par,ntot,
     &  wparcel,sig_par(nmodes)
!
! *** setup common block variables
!
       wparcel = wparc

       sinteg1=zero_par
       sinteg2=zero_par
       nrdpl = zero_par

! *** setup constants
!
       alfa_par = grav_par*amw_par*dhv_par/cpair_par/rgas_par/temp_par /
     &temp_par -grav_par*ama_par/rgas_par/temp_par

       bet1_par = pres_par*ama_par/psat_par/amw_par + amw_par*dhv_par*
     & dhv_par/cpair_par/rgas_par/temp_par/temp_par
       bet2_par = rgas_par*temp_par*denw_par/psat_par/dv_par/amw_par/
     & 4d0 +dhv_par*denw_par/4d0/aka_par/temp_par*(dhv_par* amw_par/
     &rgas_par/temp_par - 1d0)
       beta = 0.5d0*pi_par*bet1_par*denw_par/bet2_par/alfa_par/ wparc/
     &dair_par


       cf1 = 0.5d0*(((1.d0/bet2_par)/(alfa_par*wparc))**0.5d0)
       cf2 = akoh_par/3d0

!
! *** INITIAL VALUES FOR BISECTION **************************************

! *** initial values for bisection **************************************
!
       x1 = 1.0d-5

       if (ntot .gt. zero_par) then
       call sintegral (x1,ndrpl,sinteg1,sinteg2,wparcel,nmodes,
     &  alfa_par,bet2_par,akoh_par,sg_par,tp_par,sig_par)
       end if

       y1 = (sinteg1*cf1+sinteg2*cf2)*beta*x1 - 1d0
!
       x2 = 1d0
       if (ntot .gt. zero_par) then
       call sintegral (x2,ndrpl,sinteg1,sinteg2,wparcel,nmodes,
     &  alfa_par,bet2_par,akoh_par,sg_par,tp_par,sig_par)
       end if

       y2 = (sinteg1*cf1+sinteg2*cf2)*beta*x2 - 1d0
!
! *** perform bisection *************************************************
!
20     do 30 i=1,maxit_par
       x3 = 0.5*(x1+x2)
!
       if (ntot .gt. zero_par) then
       call sintegral (x3,ndrpl,sinteg1,sinteg2,wparcel,nmodes,
     &  alfa_par,bet2_par,akoh_par,sg_par,tp_par,sig_par)
       end if

       y3 = (sinteg1*cf1+sinteg2*cf2)*beta*x3 - 1d0

       if (sign(1.d0,y1)*sign(1.d0,y3) .le. zero_par) then
!                                          ! (y1*y3 .le. zero)
       y2 = y3
       x2 = x3
       else
       y1 = y3
       x1 = x3
       endif
!
       if (abs(x2-x1) .le. eps_par*x1) goto 40
       niter = i


30     continue
!
! *** converged ; return ************************************************
!
40      x3 = 0.5*(x1+x2)
!
       if (ntot .gt. zero_par) then
       call sintegral (x2,ndrpl,sinteg1,sinteg2,wparcel,nmodes,
     &  alfa_par,bet2_par,akoh_par,sg_par,tp_par,sig_par)
       end if


       y3 = (sinteg1*cf1+sinteg2*cf2)*beta*x3 - 1d0
       smax = x3

       ndroplet=ndrpl
       return
!
! *** end of subroutine activate ****************************************
!
       end subroutine activate


!=======================================================================
!
! *** subroutine sintegral
! *** this subroutine calculates the condensation integrals, according
!     to the population splitting algorithm of nenes and seinfeld (2003)
!     modal formulation according to fountoukis and nenes (2004)
!
! *** written by athanasios nenes
!
!!     Code Developer
!      Donifan Barahona
!      donifanb@umbc.edu
!=======================================================================
!
       subroutine sintegral (spar, summa, sum, summat,wparcel,nmodes,
     &  alfa_par,bet2_par,akoh_par,sg_par,tp_par,sig_par)

!
       integer ::j,i,k,nmodes
       real*8 :: sum, summat, summa, nd(nsmx_par), integ1(nsmx_par),
     & integ2(nsmx_par),alfa_par,bet2_par,akoh_par,sg_par(nmodes),
     & tp_par(nmodes),sig_par(nmodes)
       real*8 :: descr,spar,ratio, ssplt2, ssplt1, sqrt, ssplt, sqtwo,
     & dlgsg,dlgsp, ekth,wparcel
       real*8 :: orism1, orism2, orism3, orism4, orism5
!
! *** here is where the criterion with the descriminant is put. when it
!     is < 0, then set crit2 = .true. otherwise, set the two values of
!     ssplt and continue.
!
       descr = 1d0 - (16d0/9d0)*alfa_par*wparcel*bet2_par* (akoh_par/
     &spar**2)**2
!
       if (descr.le.0d0) then
       ratio = (2.0d7/3.0)*akoh_par*spar**(-0.3824)
       if (ratio.gt.1.0) ratio = 1.0
       ssplt2 = spar*ratio
       else
       ssplt1 = 0.5d0*(1d0-sqrt(descr))
       ssplt2 = 0.5d0*(1d0+sqrt(descr))
       ssplt1 = sqrt(ssplt1)*spar
       ssplt2 = sqrt(ssplt2)*spar
       endif
!
       ssplt = ssplt2
!
! *** calculate integrals
!
       sum = 0
       summat = 0
       summa = 0
!
       sqtwo = sqrt(2d0)

!
       do 999 j = 1, nmodes
!
       dlgsg = sig_par(j)
       dlgsp = dlog(sg_par(j)/spar)
       orism1 = 2.d0*dlog(sg_par(j)/ssplt2)/(3.d0*sqtwo*dlgsg )
       orism2 = orism1 - 3.d0*dlgsg/(2.d0*sqtwo)
       orism3 = 2.d0*dlgsp/(3.d0*sqtwo*dlgsg)-3.d0*dlgsg/(2.d0*sqtwo)
       orism4 = orism1 + 3.d0*dlgsg/sqtwo
       orism5 = 2.d0*dlgsp/(3*sqtwo*dlgsg)
       ekth = exp(9d0/2d0*dlgsg*dlgsg)
       integ1(j) = tp_par(j)*spar*((1-erf(orism1)) - 0.5d0*((sg_par(j)/
     &spar)**2)*ekth*(1-erf(orism4)))

       integ2(j) = (exp(9d0/8d0*dlgsg*dlgsg)*tp_par(j)/sg_par(j))*
     & (erf(orism2) - erf(orism3))
!
! *** calculate number of drops
!
       nd(j) = (tp_par(j)/2.0)*(1.0-erf(orism5))
!
       sum = sum + integ1(j)
       summat = summat + integ2(j)
       summa = summa + nd(j)
999     continue
!
       return
       end subroutine sintegral

!=======================================================================

!     DESCRIPTION
!
! *** subroutine props
! *** this subroutine calculates the thermophysical properties for the CCN activ param
!
! *** written by athanasios nenes
!      Code Developer
!      Donifan Barahona
!      donifanb@umbc.edu

!=======================================================================
!
       subroutine props(pres_par,temp_par,surt_par,dv_par,act_param,
     &  aka_par, psat_par,dair_par)
!
!
       real*8 :: presa,dbig,dlow,coef,aka_par, dv_par, psat_par
       real*8 :: pres_par,temp_par,surt_par,dair_par
       integer act_param
!
       denw_par = 1d3
       dhv_par = 2.25d6
       cpair_par = 1.0061d3
       presa = pres_par/1.013d5
       dair_par = pres_par*ama_par/rgas_par/temp_par
       aka_par = (4.39+0.071*temp_par)*1d-3
       surt_par = sft(temp_par)

       if (act_param .le. 1) then
       dv_par = (0.211d0/presa)*(temp_par/273d0)**1.94
       dv_par = dv_par*1d-4
       dbig = 5.0d-6
       dlow = 0.207683*((accom_par)**(-0.33048))
       dlow = dlow*1d-6



       coef = ((2*pi_par*amw_par/(rgas_par*temp_par))**0.5)

       dv_par = (dv_par/(dbig-dlow))*((dbig-dlow)-(2*dv_par/accom_par) *
     &coef*(dlog((dbig+(2*dv_par/accom_par)*coef)/(dlow+ (2*dv_par/
     &accom_par)*coef))))

       psat_par = vpres(temp_par)*(1e5/1.0d3)



       end if
!
       return
!
! *** end of subroutine props *******************************************
!
       end subroutine props



!PHYSICAL PROPERTIES for Nenes CDNC Activation



!=======================================================================
!
! *** function vpres
! *** this function calculates saturated water vapour pressure as a
!     function of temperature. valid for temperatures between -50 and
!     50 c.
!
! ======================== arguments / usage ===========================
!
!  input:
!     [t]
!     real variable.
!     ambient temperature expressed in kelvin.
!
!  output:
!     [vpres]
!     real variable.
!     saturated vapor pressure expressed in mbar.
!
!=======================================================================
!
       real*8 function vpres (t)
!
       integer ::i
       real*8 :: a(0:6), t,ttemp, vp
       data a/6.107799610e+0, 4.436518521e-1, 1.428945805e-2,
     & 2.650648471e-4, 3.031240396e-6, 2.034080948e-8, 6.136820929e-11/
!
! calculate polynomial (without exponentiation).
!
       ttemp = t-273.0d0
       vp = a(6)*ttemp
       do i=5,1,-1
       vp = (vp + a(i))*ttemp
       enddo
       vpres = vp + a(0)
!
! end of function vpres
!
       return
       end function vpres



!=======================================================================
!
! *** function sft
! *** this function calculates water surface tension as a
!     function of temperature. valid for temperatures between -40 and
!     40 c.
!
! ======================== arguments / usage ===========================
!
!  input:
!     [t]
!     real variable.
!     ambient temperature expressed in kelvin.
!
!  output:
!     [sft]
!     real variable.
!     surface tension expressed in j m-2.
!
!=======================================================================
!
       real*8 function sft (t)
!
       implicit none

!
       real*8 :: t,tpars
!
       tpars = t-273.15d0
       sft = 0.0761-1.55e-4*tpars
!
       return
       end function sft


! ***********************************************************************
!
       subroutine gauleg (x,w,n)
!
! calculation of points and weights for n point gauss integration
! ***********************************************************************
!
       integer :: n,m,i,j
       real*8 :: x(n), w(n),xm,xl,z,p1,p2,p3,pp,z1
       real*8, parameter :: eps_par=1.e-6
       real*8, parameter :: x1=-1.0, x2=1.0
!
! calculation
!
       m=(n+1)/2d0
       xm=0.5d0*(x2+x1)
       xl=0.5d0*(x2-x1)
       do 12 i=1,m
       z=cos(pi_par*(i-.25d0)/(n+.5d0))
 1     continue
       p1=1.d0
       p2=0.d0
       do 11 j=1,n
       p3=p2
       p2=p1
       p1=((2.d0*j-1.)*z*p2-(j-1.d0)*p3)/j
 11    continue
       pp=n*(z*p1-p2)/(z*z-1.d0)
       z1=z
       z=z1-p1/pp
       if(abs(z-z1).gt.eps_par)go to 1
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
       w(n+1-i)=w(i)
 12    continue
       return
       end subroutine gauleg

!C=======================================================================
!C
!C *** REAL FUNCTION erf (overwrites previous versions)
!C *** THIS SUBROUTINE CALCULATES THE ERROR FUNCTION USING A
!C *** POLYNOMIAL APPROXIMATION
!C
!C=======================================================================
!C
       REAL*8 FUNCTION erf(x)
       REAL*8 :: x
       REAL*8 :: AA(4), axx, y
       DATA AA /0.278393d0,0.230389d0,0.000972d0,0.078108d0/

       y = dabs(dble(x))
       axx = 1.d0 + y*(AA(1)+y*(AA(2)+y*(AA(3)+y*AA(4))))
       axx = axx*axx
       axx = axx*axx
       axx = 1.d0 - (1.d0/axx)
       if(x.le.0.) then
       erf = -axx
       else
       erf = axx
       endif
       RETURN
       END FUNCTION


!=======================================================================
!
! *** real function erf
! *** this subroutine calculates the error function
!
! *** obtained from numerical recipies
!
!=======================================================================
!

!
!      real*8  :: x

!      if(x.lt.0.)then
!        erf=-gammp(.5d0,x**2)
!      else
!        erf=gammp(.5d0,x**2)
!      endif
!      return

!      end function erf

!
!=======================================================================
!
       real*8 function gammln(xx)
!
!=======================================================================
!
!
       integer :: j
       real*8 :: cof(6),stp,half,one,fpf,x,tmp,ser,xx
!
       data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0, -
     &1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
       data half,one,fpf/0.5d0,1.0d0,5.5d0/
       x=xx-one
       tmp=x+fpf
       tmp=(x+half)*log(tmp)-tmp
       ser=one
       do 11 j=1,6
       x=x+one
       ser=ser+cof(j)/x
 11    continue
       gammln=tmp+log(stp*ser)
       return
       end function gammln


!
!=======================================================================
!
!      real*8 function gammp(a,x)
!
!=======================================================================
!
!     real*8 :: a,x,gln,gamser,gammcf

!     if(x.lt.0.d0.or.a.le.0.d0)pause
!     if(x.lt.a+1.d0)then
!       call gser(gamser,a,x,gln)
!       gammp=gamser
!     else
!       call gcf(gammcf,a,x,gln)
!       gammp=1.d0-gammcf
!     endif
!     return
!     end function gammp


!
!=======================================================================
!
!     subroutine gcf(gammcf,a,x,gln)
!
!=======================================================================
!
!
!     integer            :: n
!     integer, parameter :: itmax=100
!     real*8, parameter  :: eps_par=3.e-7
!       real*8             :: gln, gold,a0,a1,x,b0,b1,fac,an, &
!                          float,ana,a,anf,g,gammcf
!     gln=gammln(a)
!     gold=0.
!     a0=1.
!     a1=x
!     b0=0.
!     b1=1.
!     fac=1.
!     do 11 n=1,itmax
!       an=float(n)
!       ana=an-a
!       a0=(a1+a0*ana)*fac
!       b0=(b1+b0*ana)*fac
!       anf=an*fac
!       a1=x*a0+anf*a1
!       b1=x*b0+anf*b1
!       if(a1.ne.0.)then
!         fac=1./a1
!         g=b1*fac
!         if(abs((g-gold)/g).lt.eps_par)go to 1
!         gold=g
!       endif
!1    continue
!     pause 'a too large, itmax too small'
!     gammcf=exp(-x+a*log(x)-gln)*g
!     return
!     end subroutine gcf


!
!=======================================================================
!sft
!     subroutine gser(gamser,a,x,gln)
!
!=======================================================================
!

!     integer  :: n
!     integer, parameter :: itmax=100
!     real*8, parameter  :: eps_par=3.e-7
!       real*8  :: gln,x,gamser,a,ap,sum,del,abs

!     gln=gammln(a)
!     if(x.le.0.)then
!       if(x.lt.0.)pause
!       gamser=0.
!       return
!     endif
!     ap=a
!     sum=1./a
!     del=sum
!     do 11 n=1,itmax
!       ap=ap+1.
!       del=del*x/ap
!       sum=sum+del
!       if(abs(del).lt.abs(sum)*eps_par)go to 1
!1    continue
!     pause 'a too large, itmax too small'
!     gamser=sum*exp(-x+a*log(x)-gln)
!     return
!     end subroutine gser

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

!

!*************************************************************************
!
! ICE NUCLEATION PARAMETERIZATION FILES START HERE
!
! ************************************************************************
!
!
!======================================================================
!
!       Code Developer
!       Donifan Barahona, GA TECH
!       donifan@gatech.edu
!    -------------------------------
!     DESCRIPTION
!

!***********************************************************
!** Parameterization  of ICE crystal number concentration
!** for large scale models.
!** Donifan Barahona, Athanasios Nenes
!   JGR, 111, D11211,  2008
!   ACP, 9, 369-381,   2009
!   ACP  9, 5933-5948, 2009
!   Homogeneoeus and heterogeneous nucleation considered
!*** SI units unless otherwise specified.
!
!
! *** WRITTEN BY DONIFAN BARAHONA
!
!=======================================================================



       SUBROUTINE IceParam (sigma_w, denice_ice,ddry_ice,np_ice,
     & nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,
     & g1_ice, g2_ice,gdoin_ice, z_ice,lambda_ice,sc_ice,
     & norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,one_over_tao_preex,
     & nhet, nice, smax, nlim,Nhet_dep,Nhet_dhf,fdust_dep
     & ,P_ice, T_ice,ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,
     & use_av_v,miuv_ice,vpresw_ice,vpresi_ice,denair_ice)

       real*8, intent(in) :: sigma_w,denice_ice,ddry_ice,np_ice,
     & nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,
     & g1_ice, g2_ice,gdoin_ice, z_ice, lambda_ice,
     & norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,sc_ice,one_over_tao_preex,P_ice, T_ice,
     & vpresw_ice,vpresi_ice,denair_ice
       real*8, dimension(:) :: ndust_ice, sigdust_ice,ddust_ice


       real*8, intent(out) :: nice, nhet, smax, nlim,Nhet_dep,Nhet_dhf
     & ,fdust_dep

       real*8 :: wpar_ice, sigmav_ice,vmax_ice,miuv_ice
       integer, intent(in)::nbindust_ice
       logical use_av_v




       vmin_ice=0.005d0
       vmax_ice=2.5d0
       sigmav_ice=sigma_w
       vmax_ice= max(min(miuv_ice+(4d0*sigmav_ice), vmax_ice),
     & vmin_ice +0.01)

       if ((sigmav_ice .lt. 0.05) .or. (T_ice .gt. Thom)) then
       use_av_v= .TRUE.
       end if

       if (vmax_ice .gt. vmin_ice + 0.01) then
       if (use_av_v) then
       wpar_ice=min(max(miuv_ice + sigma_w*0.8, vmin_ice), vmax_ice)
       CALL nice_param(wpar_ice, denice_ice,ddry_ice,np_ice,nin_ice,
     & alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,miuv_ice,
     & sigmav_ice,g1_ice, g2_ice,gdoin_ice, z_ice,vmax_ice,
     & norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,lambda_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,sc_ice,one_over_tao_preex,
     & nice, smax, nhet, nlim,Nhet_dep,Nhet_dhf,fdust_dep,
     & P_ice, T_ice,ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,
     & vpresw_ice,vpresi_ice,denair_ice)
       else
       call nice_Vdist(denice_ice,ddry_ice,np_ice,
     & nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,miuv_ice,
     & sigmav_ice,g1_ice, g2_ice,gdoin_ice, z_ice,vmax_ice,
     & norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,lambda_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,sc_ice,one_over_tao_preex,
     & nice, smax, nhet, nlim,Nhet_dep,Nhet_dhf,fdust_dep,
     & P_ice, T_ice,ndust_ice, sigdust_ice,ddust_ice,nbindust_ice
     & ,vpresw_ice,vpresi_ice,denair_ice)
       end IF
       else
       nice = zero_par
       nhet = zero_par
       nlim = 1.0e10
       smax = zero_par
       end if


       return
       END subroutine IceParam


!*************************************************************
!    Subroutine nice_Vdist. Calculates the ice crystal number concentration
!    at the maximum supersaturation using a PDF of updraft using a
!    sixth order Gauss-Legendre quadrature
!     Inputs:  T, and P all SI units)
!    Output NC, smax, nhet, nlim (m-3)
!    Barahona and Nenes, JGR, D11211 (2008) and ACPD, 15665-15698, (2008)
!   Written by Donifan Barahona
!************************************************************


       subroutine nice_Vdist(denice_ice,ddry_ice,np_ice,
     & nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,miuv_ice,
     & sigmav_ice,g1_ice, g2_ice,gdoin_ice, z_ice,vmax_ice,
     & norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,lambda_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,sc_ice,one_over_tao_preex,
     & nice, smax, nhet, nlim,Nhet_dep,Nhet_dhf,fdust_dep,
     & P_ice, T_ice,ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,
     & vpresw_ice,vpresi_ice,denair_ice)

       real*8 :: quadx(6), wpar, sum1, quadw(6), dp, sum2, sum3, sum4,
     & sum5, x1, x2
       real*8, intent(out) :: nice, smax, nhet, nlim,Nhet_dep,Nhet_dhf
     &  ,fdust_dep
       real*8  wpar_icex,denice_ice,ddry_ice,np_ice,
     & nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,miuv_ice,
     & normv_ice,sigmav_ice,g1_ice, g2_ice,gdoin_ice, z_ice,vmax_ice,
     & norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,lambda_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,sc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,one_over_tao_preex,P_ice, T_ice,
     & vpresw_ice,vpresi_ice,denair_ice
       real*8, dimension(:) :: ndust_ice, sigdust_ice,ddust_ice

       INTEGER :: INDEX,nbindust_ice

       DATA quadx/0.23861918d0, -0.23861918d0, 0.66120939d0, -
     &0.66120939d0, 0.93246951d0, -0.93246951d0/

       DATA quadw/0.46791393d0, 0.46791393d0, 0.36076157d0,
     & 0.36076157d0, 0.17132449d0, 0.17132449d0/



       x1=(vmin_ice-miuv_ice)/(sq2_par*sigmav_ice)
       x2=(vmax_ice-miuv_ice)/(sq2_par*sigmav_ice)

       x2=max(x1 +0.01, x2)
       normv_ice=(ERFAPP(x2)-ERFAPP(x1))*0.5d0

       sum1=0d0
       sum2=0d0
       sum3=0d0
       sum4=0d0
       sum5=0d0

       DO INDEX =1, 6
       wpar=max(0.5d0*(((vmax_ice-vmin_ice)*quadx(INDEX)) +(vmax_ice+
     &      vmin_ice)), 0.01)


       CALL nice_param(wpar, denice_ice,ddry_ice,np_ice,nin_ice,
     & alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,miuv_ice,
     & sigmav_ice,g1_ice, g2_ice,gdoin_ice, z_ice,vmax_ice,
     & norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,lambda_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,sc_ice,one_over_tao_preex,
     & nice, smax, nhet, nlim,Nhet_dep,Nhet_dhf,fdust_dep,
     & P_ice, T_ice,ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,
     & vpresw_ice,vpresi_ice,denair_ice)

       CALL gausspdf(wpar, dp, sigmav_ice,miuv_ice,normv_ice)
       sum1=sum1+(nice*dp*quadw(INDEX))
       sum2=sum2+(smax*dp*quadw(INDEX))
       sum3=sum3+(nhet*dp*quadw(INDEX))
       sum4=sum4+(nlim*dp*quadw(INDEX))
       sum5=sum5+(sc_ice*dp*quadw(INDEX))



       END DO
       nice=sum1*(vmax_ice-vmin_ice)*0.5d0
       smax=sum2*(vmax_ice-vmin_ice)*0.5d0
       nhet=sum3*(vmax_ice-vmin_ice)*0.5d0
       nlim=sum4*(vmax_ice-vmin_ice)*0.5d0
       sc_ice=sum5*(vmax_ice-vmin_ice)*0.5d0
       RETURN


       END subroutine nice_Vdist

!*************************************************************

!*************************************************************
       real*8 function ERFAPP(x)

       real*8, intent(in) :: x
       real*8 :: a
       a = x*x*(1.27324d0+(0.147d0*x*x))/(1d0+(0.147d0*x*x))
       ERFAPP = SQRT(1d0-exp(-a))

       if (x .lt. 0.0) then
         ERFAPP = - ERFAPP
       end if

       end function ERFAPP



!*************************************************************
!    Subroutine Het_freezing. Use only for mixed phase clouds . Inputs: Wpar, T, and P all SI units)
!    Output Nc=Nhet (m-3)
! !    Written by Donifan Barahona
!************************************************************


!      SUBROUTINE Het_freezing(nhet, nice, smax)

!      real*8, intent(out) :: nice, nhet, smax
!      real*8, intent(in) :: np_ice

!      real*8 :: I, SX, NHET_, DSH

!      nhet=zero_par
!      smax=zero_par
!      nice=zero_par

!      SX=vpresw_ice/vpresi_ice
!      call INSPEC_ice(SX-1.0, NHET_, DSH,np_ice)
!      sc_ice=1.0
!      nhet=max(NHET_, zero_par)
!      nice=nhet
!      smax=max(SX-1.0, zero_par)


!      END SUBROUTINE Het_freezing


!*************************************************************
!    Subroutine nice_param. This is the implementation of the Barahona and Nenes(2008, 2009)
!    parameterization
!    Written by Donifan Barahona
!************************************************************


       subroutine nice_param(wpar_icex,denice_ice,ddry_ice,np_ice,
     & nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,miuv_ice,
     & sigmav_ice,g1_ice, g2_ice,gdoin_ice, z_ice,vmax_ice,
     & norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,lambda_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,sc_ice,one_over_tao_preex,
     & nice, smax, nhet, nlim_,Nhet_dep,Nhet_dhf,fdust_dep,
     & P_ice, T_ice,ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,
     & vpresw_ice,vpresi_ice,denair_ice)

       real*8, intent(in) :: wpar_icex,denice_ice,np_ice,nin_ice,
     &  alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,miuv_ice,
     &  sigmav_ice,g1_ice, g2_ice,gdoin_ice, z_ice,vmax_ice,
     &  norg_ice,sigorg_ice,dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice, ddry_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,one_over_tao_preex,P_ice, T_ice,
     & vpresw_ice,vpresi_ice,denair_ice
       real*8, intent(out) :: nice, smax, nhet, nlim_,Nhet_dep,Nhet_dhf
     & ,fdust_dep

       real*8 :: AUX1, AUX2, G, DPMAX, MIU, MONVOL, FDS, NLIM, DLIM,
     & DSTAR, DS, NSTAR, NHOM, FC, PHIDO, AUXNC, SIZECORR, DSH, NHET_,
     & F1, F2, SAUX, SUP, SLOW, SHALF, FHALF, DPMIN, GAM, wpar_ice,
     & preex_effect, swsat, sstep,sc_ice,lambda_ice
       real*8, dimension(:) :: ndust_ice, sigdust_ice,ddust_ice

       integer :: INDEX, maxiter_s,nbindust_ice





       preex_effect=1.0-(one_over_tao_preex*(shom_ice)/alfa_ice/
     & ((shom_ice+1.0))/wpar_icex)
       swsat = vpresw_ice/vpresi_ice -1d0


       if (preex_effect .le. 0.0) then
       nhet=0d0
       NHOM=0d0
       smax=shom_ice
       DSH =0.d0
       FDS=1.d0

       sc_ice = shom_ice + 1.d0
       nice = 0.d0
       nlim_=0d0
       return
       else

       wpar_ice = wpar_icex*preex_effect
       end if




       if (np_ice .gt. 1.0) then

       MONVOL=np_ice*1.0d-6*ddry_ice*ddry_ice*ddry_ice
       AUX1=1.6397d-14*T_ice-3.1769d-12
       DPMAX=AUX1*(MONVOL**(-0.373d0))*(wpar_ice**(-0.05))
       IF (DPMAX.GT.1.0d-4) THEN
       DPMAX=1.0d-4
       END IF

       else
       DPMAX=1.0d-4

       endif


       DPMIN=dliq_ice+(0.02/sqrt(alfa_ice*wpar_ice*g1_ice))
       DPMAX=max(DPMAX,DPMIN)

       AUX1=DPMAX-dliq_ice
       AUX2=dlog((g2_ice+(g1_ice*DPMAX))/(g2_ice+(g1_ice*dliq_ice)))
       G=1.3346d0*((g1_ice*AUX1)-(g2_ice*AUX2))/(AUX1*g1_ice*g1_ice)
       lambda_ice=lambda_ice/sqrt(wpar_ice)
       AUX1 = g1_ice*alfa_ice*wpar_ice
       NSTAR=(AUX1*SQRT(AUX1))/beta_ice/z_ice/sq2_par

       GAM=g2_ice/g1_ice




       FDS=1d0
       NHOM=0d0



       if (typeofspec_ice .ge. 0d0) then

       call INSPEC_ice(shom_ice,NHET_,DSH,np_ice,norg_ice,sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     &   kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     & ndust_ice, sigdust_ice,ddust_ice,nbindust_ice
     & ,vpresw_ice,vpresi_ice)
       SIZECORR=EXP(-2d0/lambda_ice/shom_ice)
       DSTAR=((4d0*DSH*DSH/3d0)+(2d0*DSH*(shom_ice-DSH))) /(shom_ice-
     &DSH+1d0)

       if (DSTAR .gt. 0.0) then
       NLIM=min(NSTAR*(shom_ice+1d0)/shom_ice/sqrt(DSTAR)/SIZECORR,
     & 1.e10)
       else
       NLIM=1.e10
       end if

       else
       DSH=shom_ice-sh_ice
       DSTAR=((4d0*DSH*DSH/3d0)+(2d0*DSH*(shom_ice-DSH))) /(shom_ice-
     &DSH+1d0)
       DLIM=-GAM+sqrt((GAM*GAM)+(2d0*DSTAR/g1_ice/alfa_ice/wpar_ice))
       NLIM=alfa_ice*wpar_ice*(shom_ice+1d0)/z_ice/beta_ice/shom_ice
       NLIM=NLIM*((g1_ice*DLIM)+g2_ice)/DLIM/DLIM
       NHET_=nin_ice
       end if


       nlim_=min(NLIM, 1d10)
       nlim_=max(NLIM, 1d-6)

       if (NHET_ .gt. 0.0) then
       AUX1 =NHET_/nlim_
       FDS=1d0-(AUX1*SQRT(AUX1))
       else
       FDS = 1d0
       end if



       if (purehom_ice) then
       FDS=1.0d0
       end if

       if ((purehet_ice) .or. (T_ice .GE. Thom)) then
       FDS = 0d0
       end if


       IF (FDS .GE. 1.e-10) THEN

       MIU=FDS*alfa_ice*(shom_ice+1d0)*wpar_ice*koft_ice/shom_ice

       PHIDO=sqrt(pi_ice*G/MIU/2d0)*(G/MIU)
       AUXNC=2d0*denair_ice/beta_ice/koft_ice/denice_ice/pi_ice/np_ice
       FC=AUXNC/PHIDO



       if (np_ice .gt.0d0) then

       IF (FC .le. 0.6d0) then
       NHOM=np_ice*EXP(-FC)*(1.0d0-EXP(-FC))
       else
       NHOM=np_ice/(1d0+EXP((9d0-2d0*FC)/7d0))
       end if

       else
       NHOM=0d0
       end if

       smax=shom_ice
       nhet=NHET_
       if (purehom_ice) NHET_ = 0.d0


       ELSE


       NHOM = 0d0


       smax=0d0
       nhet=0d0
       SAUX=0.01d0

       if (typeofspec_ice .lt. 0d0) SAUX=sh_ice+0.00000000001d0

       F1= FINDSMAX(SAUX,DSH,
     & NHET_,np_ice,norg_ice,sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     &  ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,lambda_ice,
     &  gdoin_ice,alfa_ice,wpar_ice,GAM,g1_ice,g2_ice,z_ice,
     &  beta_ice,nin_ice,vpresw_ice,vpresi_ice,NSTAR)
       F2=1.0
       sstep = 0.05d0
       maxiter_s = int(1d0/sstep) + 1

       do INDEX =1, maxiter_s

       if (SAUX .ge. swsat) then
       F2=F1
       SAUX = swsat
       exit
       end if


       SAUX=SAUX+sstep
       F2=FINDSMAX(SAUX,DSH,
     & NHET_,np_ice,norg_ice,sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     &  ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,lambda_ice,
     &  gdoin_ice,alfa_ice,wpar_ice,GAM,g1_ice,g2_ice,z_ice,
     &  beta_ice,nin_ice,vpresw_ice,vpresi_ice,NSTAR)
       IF (F2*F1 .lt. 0d0) exit

       end do

       if (F2*F1 .gt.0d0) then
       nhet=0d0
       smax=SAUX
       else

       if (SAUX .lt. swsat) then


       SUP=SAUX
       SLOW=SAUX-(sstep + 0.001d0)

       DO INDEX=1,50
       SHALF=0.5d0*(SUP+SLOW)
       FHALF=FINDSMAX(SHALF,DSH,
     & NHET_,np_ice,norg_ice,sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     &  ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,lambda_ice,
     &  gdoin_ice,alfa_ice,wpar_ice,GAM,g1_ice,g2_ice,z_ice,
     &  beta_ice,nin_ice,vpresw_ice,vpresi_ice,NSTAR)

       IF (SIGN(1.d0,F1)*SIGN(1.d0,FHALF) .LE. 0d0) THEN
       F2 = FHALF
       SUP = SHALF
       ELSE
       F1 = FHALF
       SLOW = SHALF
       ENDIF


       IF (ABS(SLOW-SUP) .LE. 5d-3) exit
       END DO
       else
       SHALF = swsat
       end if

       smax=SHALF

       end if

       if ((smax .gt. shom_ice) .and.(T_ice .LE. Thom)) smax=shom_ice

       if (typeofspec_ice .ge. 0d0) then
       call INSPEC_ice(smax, nhet,DSH,np_ice,norg_ice,sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     & ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,
     & vpresw_ice,vpresi_ice)
       else
       nhet=nin_ice
       end if

       END IF


       nice=NHOM+nhet
       sc_ice=max(smax+1.0-DSH, 1.0)


      if (.true.) then

       if (FDS .gt. 0.0) then
       sc_ice = (shom_ice+1.0)*FDS + sc_ice*(1.0-FDS)
       end if
       end if

       sc_ice=min(shom_ice+1.0, sc_ice)


       return

       END subroutine nice_param
!*************************************************************
       real*8 function FINDSMAX(SX,DSH,
     & NHET_,np_ice,norg_ice,sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     &  ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,lambda_ice,
     &  gdoin_ice,alfa_ice,wpar_ice,GAM,g1_ice,g2_ice,z_ice,
     &  beta_ice,nin_ice,vpresw_ice,vpresi_ice,NSTAR)
       real*8, intent(in) :: SX
       real*8 :: tao
       real*8 ::  DSTAR, SIZECORR, DSH, NSTAR,DLIM
       integer nbindust_ice
       real*8, dimension(:) :: ndust_ice, sigdust_ice,ddust_ice
       real*8 shom_ice,NHET_,np_ice,norg_ice,sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     &  lambda_ice,
     &  gdoin_ice,alfa_ice,wpar_ice,GAM,g1_ice,g2_ice,z_ice,
     &  beta_ice,nin_ice,vpresw_ice,vpresi_ice

       if (typeofspec_ice .ge. 0d0) then

       call INSPEC_ice(SX,NHET_,DSH,np_ice,norg_ice,sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     &  ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,
     & vpresw_ice,vpresi_ice)
       SIZECORR=EXP(-2d0/lambda_ice/SX)
       DSTAR=((4d0*DSH*DSH/3d0)+(2d0*DSH*(SX-DSH)))/(SX-DSH+1d0)
       DSTAR=DSTAR+(gdoin_ice*alfa_ice*wpar_ice)


       tao=NHET_*SIZECORR*SX*sqrt(DSTAR)/(SX+1d0)/NSTAR
      

       else
       DSH=SX-sh_ice
       DSTAR=((4d0*DSH*DSH/3d0)+(2d0*DSH*(SX-DSH))) /(SX-DSH+1d0)
       DLIM=-GAM+sqrt((GAM*GAM)+(2d0*DSTAR/g1_ice/alfa_ice/wpar_ice))
       tao=alfa_ice*wpar_ice*(SX+1d0)/z_ice/beta_ice/SX
       tao=tao*((g1_ice*DLIM)+g2_ice)/DLIM/DLIM/nin_ice

       end if


       FINDSMAX=1d0-tao

       end function FINDSMAX


!*************************************************************
!    Function VPRESWATER. Calculates the saturated vapor pressure
!    of water (Pa) according to Murphy & Koop (2005)
!    T in K (173.15-373.15)
!************************************************************

       real*8 function VPRESWATER_ice(T)

       real*8, intent(in) :: T
       real*8 :: A(0:9)

       DATA A/54.842763d0, -6763.22d0, -4.21d0, 0.000367d0, 0.0415d0,
     & 218.8d0, 53.878d0, -1331.22d0, -9.44523d0, 0.014025d0/


       VPRESWATER_ice = A(0)+(A(1)/T)+(A(2)*LOG(T))+(A(3)*T)+
     & (TANH(A(4)*(T-A(5)))*((A(6)+(A(7)/T))+ (A(8)*LOG(T))+ (A(9)*T)))

       VPRESWATER_ice=EXP(VPRESWATER_ice)

       return
       END function VPRESWATER_ice

!*************************************************************
!    Function VPRESICE. Calculates the saturated vapor pressure
!    of ice (pa) according to Murphy & Koop (2005)
!    T in K (>110)
!************************************************************

       real*8 function VPRESICE(T)

       real*8, intent(in) :: T
       real*8 :: A(0:3)

       DATA A/9.550426d0, -5723.265d0, 3.53068d0, -0.00728332d0/


       VPRESICE = A(0)+(A(1)/T)+(A(2)*LOG(T))+(A(3)*T)
       VPRESICE=EXP(VPRESICE)

       return
       END function VPRESICE

!*************************************************************
!    Function DHSUB. Calculates the latent heat of sublimation
!    of ice (J/Kg) according to Murphy & Koop (2005)
!    T in K (>30)
!************************************************************

       real*8 function DHSUB_ice(T)

       real*8, intent(in) :: T
       real*8 :: A(0:4)


       DATA A/46782.5d0, 35.8925d0, -0.07414d0, 541.5d0, 123.75d0/

       DHSUB_ice = A(0) + (A(1) * T) + (A(2)*T*T) + (A(3) * EXP(-((T/
     & A(4))**2)))

       DHSUB_ice=1000d0*DHSUB_ice/18d0
       return
       END function DHSUB_ice

!*************************************************************
!    Function ICEDENSITY. Calculates the DENSITY OF ICE
!    of ice (Kg/m3) according to PK97
!    T in K (>30)
!************************************************************

       real*8 function DENSITYICE(T)

       real*8, intent(in) :: T
       real*8 :: A(0:2), TTEMP

       DATA A/0.9167d0, -1.75d-4, -5.0d-7/

       TTEMP=T-273d0

       DENSITYICE= 1000d0*(A(0)+(A(1)*TTEMP)+(A(2)*TTEMP*TTEMP))
       return
       END function DENSITYICE

!*************************************************************
!    Function WATDENSITY. Calculates the DENSITY OF ICE
!    of liquid water (Kg/m3) according to PK97
!    T in K (>240)
!************************************************************

       real*8 function WATDENSITY_ice(T)

       real*8, intent(in) :: T
       real*8 :: A(0:6), TTEMP, WATDENSITY
       INTEGER :: I


       DATA A/0.99986d0, 6.690d-5, -8.486d-6, 1.518d-7, -6.9984d-9, -
     &3.6449d-10, -7.497d-12 /

       TTEMP=T-273D0

       IF (TTEMP .le. -40d0) THEN
       TTEMP=-40d0
       END IF

       WATDENSITY=A(6)*TTEMP

       IF (T .GE. 240.0) THEN
       DO I=5,1, -1
       WATDENSITY= (WATDENSITY+A(I))*(TTEMP)
       ENDDO
       WATDENSITY=WATDENSITY + A(0)
       ELSE
       WATDENSITY=0.979d0
       END IF

       WATDENSITY=WATDENSITY*1000d0
       WATDENSITY_ice=WATDENSITY
       return
       END function WATDENSITY_ice


!*************************************************************
!    Subroutine PROPERTIES. Set physical an thermodynamic
!    properties at T and P for ice param
!************************************************************


       SUBROUTINE prop_ice(T, P, denice_ice,ddry_ice,
     & nin_ice,alfa_ice,beta_ice,shom_ice, koft_ice, dliq_ice,
     & g1_ice, g2_ice,gdoin_ice, z_ice,lambda_ice,
     & kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,D_preex, N_preex,one_over_tao_preex,
     & P_ice, T_ice,act_param,ndust_ice,vpresw_ice,vpresi_ice,
     & denair_ice)

       real*8, intent(in)  :: T, P,ddry_ice,nbc_ice,D_preex, N_preex
       real*8, intent(out) :: alfa_ice,beta_ice,shom_ice,
     &  koft_ice, dliq_ice,g1_ice, g2_ice,gdoin_ice, z_ice,
     &  lambda_ice,kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,one_over_tao_preex,P_ice, T_ice

       real*8 :: AUX, AUX1, AUX2, SW, fice, mice, Tc, hdust, hbc, b0,
     &           b1, b2, b3, x, T0bc, T0dust, gam, gamma,NIN_ICE,
     &           Tr, vw, den_m, Eact, Toact, Dact, acc, n1, Siw, rgo,
     &           ngo, mw_molec,dhs_ice,denwat_ice,denice_ice,vpresw_ice,
     &           vpresi_ice,denair_ice,diff_ice,aircond_ice

       real*8,dimension(:) :: ndust_ice
       integer, intent(in) :: act_param


       T_ice    = min(max(T, Tmin_ice), To_ice)
       P_ice    = max(P, Pmin_ice)
       dliq_ice = 1.0d-6


!      rv_ice     = rgas_ice/wmw_ice
       dhs_ice    = DHSUB_ice(T_ice)
       vpresw_ice = VPRESWATER_ice(T_ice)
       vpresi_ice = VPRESICE(T_ice)
       denice_ice = DENSITYICE(T_ice)
       denwat_ice = WATDENSITY_ice(T_ice)
       denair_ice = P_ice*amw_ice/rgas_ice/T_ice


       diff_ice=(0.211d0*101325d0/P_ice)*((T_ice/273d0)**1.94d0)*1.0d-4
       AUX1=1.0e-3*(4.39d0+0.071d0*T_ice)


       AUX2=(2d0*AUX1/(thaccom_ice*1.0d-6*denair_ice*cpa_ice)) *((58.0d-
     &3*pi_ice/(rgas_ice*T_ice))**0.5d0)

       aircond_ice=AUX1/(1.d0+AUX2)




       AUX1=grav_ice*dhs_ice/rv_ice/T_ice/T_ice/cpa_ice
       AUX2=grav_ice*amw_ice/rgas_ice/T_ice
       alfa_ice=AUX1-AUX2
       beta_ice=amw_ice*P_ice/wmw_ice/vpresi_ice
       gamma=1.5d0*dhs_ice*dhs_ice/rv_ice/T_ice/T_ice/cpa_ice

       beta_ice=beta_ice+gamma



       shom_ice=2.349d0-(T_ice/259d0)
       SW=shom_ice*vpresi_ice/vpresw_ice
       shom_ice=shom_ice-1d0
       koft_ice=(0.0240d0*T_ice*T_ice)-(8.035d0*T_ice)+934.0d0





       if (SW .lt. 0.99) then
       AUX1=(1d0/(1d0-SW))-1.1764d0
       else
       AUX1=(1d0/0.01)-1.1764d0
       end if
       dliq_ice=ddry_ice*0.9344d0*(AUX1**0.333)



       AUX1=denice_ice*rv_ice*T_ice/vpresi_ice/diff_ice
       AUX2=dhs_ice*denice_ice/aircond_ice/T_ice
       AUX2=AUX2*((dhs_ice/rv_ice/T_ice)-1.0d0)
       g1_ice=(AUX1+AUX2)/4.0d0


       one_over_tao_preex = beta_ice*denice_ice*pi_ice*0.5* D_preex*
     &N_preex/g1_ice/denair_ice



       g2_ice=denice_ice*rv_ice*T_ice/2.0d0/vpresi_ice/depcoef_ice
       g2_ice=g2_ice*((2.0d0*pi_ice/rv_ice/T_ice)**0.5d0)

       doin_ice=1.0d-6
       gdoin_ice=(g1_ice*0.5d0*doin_ice*doin_ice)+(g2_ice*doin_ice)
       z_ice=denice_ice*pi_ice/2.0d0/denair_ice

       gam=g2_ice/g1_ice
       lambda_ice=1d0/sqrt(alfa_ice*g1_ice*gam*gam)





       if (typeofspec_ice .lt. 0) then
       sh_ice=0.3d0
       nin_ice=(sum(ndust_ice)+nbc_ice)*0.05d0
       elseif (typeofspec_ice .eq. 3) then

       shdust_ice = 0.2d0
       effdust_ice=0.6d0
       shbc_ice = 0.35d0
       effbc_ice=0.05d0
       mice = 0.96d0
       fice=0.25d0*((mice*mice*mice)-(3d0*mice)+2d0)
       kdust_ice=koft_ice*fice
       mice = 0.76d0
       fice=0.25d0*((mice*mice*mice)-(3d0*mice)+2d0)
       kbc_ice=koft_ice*fice
       elseif (typeofspec_ice .eq. 4) then

       Tc=T_ice-273.15d0
       hdust=0.15d0
       T0dust=-40d0
       b0=-1.0261d0; b1=3.1656d-3; b2=5.3938d-4; b3=8.2584d-6
       x=b0+(b1*Tc)+(b2*Tc*Tc)+(b3*Tc*Tc*Tc)
       si0dust_ice=1d0+(10d0**x)
       del1dust_ice=cubicint_ice(Tc, T0dust, T0dust+5d0, 1d0, hdust)
       hbc=0d0
       T0bc=-50d0
       b0=0.5652d0; b1=1.085d-2; b2=-3.118d-5
       si0bc_ice=b0+(b1*T_ice)+(b2*T_ice*T_ice)-0.1d0
       del1bc_ice=cubicint_ice(Tc, T0bc, T0bc+5d0, 1d0, hbc)
       end if

       RETURN

       END SUBROUTINE prop_ice

!*************************************************************
!   Subroutine gauspdf (normalized width of the updraft distribution).
!************************************************************

       SUBROUTINE gausspdf(x, dp, sigmav_ice,miuv_ice,normv_ice)


       real*8, intent(in) :: x
       real*8 sigmav_ice,miuv_ice,normv_ice
       real*8, intent(out) :: dp

       sigmav_ice =max(sigmav_ice, 0.01)
       normv_ice =max(normv_ice, 0.01)


       dp=EXP(-0.5d0*(x-miuv_ice)*(x-miuv_ice)/sigmav_ice/sigmav_ice) /
     &sigmav_ice/sq2pi_par/(normv_ice + 0.001)


       RETURN


       END SUBROUTINE gausspdf


!*************************************************************
!   Function cubicint_ice (cubic interpolation between y1 and y2 within a and b).
!************************************************************

       real*8 function cubicint_ice(y, y1, y2, a, b)

       real*8, intent(in) :: y, y1, y2, a, b
       real*8 :: A_, B_, a0, a1, a2, a3, d, AUX

       if (y .le. y1) then
       d=a
       goto 5065
       end if

       if (y .ge. y2) then
       d=b
       goto 5065
       end if


       AUX=y2-y1
       A_=6d0*(a-b)/(AUX*AUX*AUX)
       B_=a+(A_*(y1*y1*y1)/6d0)-(A_*(y1*y1)*y2*0.5d0)

       a0=B_
       a1=A_*y1*y2
       a2=-A_*(y1+y2)*0.5d0
       a3=A_/3d0
       d=a0+(a1*y)+(a2*y*y)+(a3*y*y*y)


 5065    cubicint_ice=d


       end function cubicint_ice

!*************************************************************
!   Function dcubicint_ice (used in the PDA08 spectrum).
!************************************************************

       real*8 function dcubicint_ice(y, y1, y2, a, b)

       real*8, intent(in) :: y, y1, y2, a, b
       real*8 :: A_, a0, a1, a2, a3, d, AUX

       if (y .le. y1) then
       d=0
       goto 5065
       end if

       if (y .ge. y2) then
       d=0
       goto 5065
       end if


       AUX=y2-y1
       A_=6d0*(a-b)/(AUX*AUX*AUX)

       a1=A_*y1*y2
       a2=-A_*(y1+y2)*0.5d0
       a3=A_/3d0
       d=(a1)+(2d0*a2*y)+(3d0*a3*y*y)


 5065   dcubicint_ice=d


       end function dcubicint_ice

!*************************************************************
! Function PDG07 (simplified ice nucleation
!                     spectra according to Phillips et. al. 2007).
! si is supersaturation wrt ice and T is in K
!************************************************************

       real*8 function PDG07_ice(si, T)

       real*8, intent(in) :: si, T
       real*8 :: N

       if (T .le. 243d0)then
       N=1000d0*exp(-0.388d0)*(exp(3.88d0*si)-1d0)
       else
       N=60d0*exp(-0.639d0)*(exp(12.96d0*si)-1d0)
       end if

       PDG07_ice=N

       end function PDG07_ice




!*************************************************************
! Subroutine INSPEC_ice
!  Provides the Ice Nuclei concentration (m-3)
! and the chracteristic freezing threeshold, DSh (Barahona & Nenes 2009), at given
! si and T. The variable typeofspec_ice (integer) has the values
! 1 Meyers et. al. 1992
! 2  Phillips et. al. 2007
! 3  Barahona 2011
! 4  Phillips et. al. 2008 (simplifed)
! 5  Phillips et. al. 2013 (simplifed)
! si is supersaturation wrt ice and T is in K

!      Written by Donifan Barahona
!      donifan.o.barahona@nasa.gov

!************************************************************

       subroutine INSPEC_ice(six, N, Dsh,np_ice,norg_ice, sigorg_ice,
     &   dorg_ice, dbc_ice,sigbc_ice,
     &   kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice,
     & ndust_ice, sigdust_ice,ddust_ice,nbindust_ice,
     & vpresw_ice,vpresi_ice)

       real*8, intent(in) :: six,np_ice,norg_ice, sigorg_ice,
     &  dorg_ice, dbc_ice,sigbc_ice,
     &  kdust_ice, kbc_ice, shdust_ice, shbc_ice,
     & effdust_ice, effbc_ice, del1dust_ice, si0dust_ice, del1bc_ice,
     & si0bc_ice,nbc_ice,P_ice, T_ice
       real*8, intent(out) :: N, Dsh,Nhet_dep,Nhet_dhf,fdust_dep
       real*8 :: Nd, Nbc, aux, Si_, SW, del0, ddel0, fc, delw0, ddelw0,
     & SW0, Hdust, Hbc, Nbase, dNd, dNbc, dNbase, dH, dfc, Ndaux,
     & dNdaux, dNorg, Norg, Ndustaux, frac, aux2, Dx2, fdep, Ndep, Ndhf,
     & dNdep, dNdhf, si, dfrac, Ncdep_, Ncdhf_ , fglassy, Nglassy,
     & dNglassy, SIW, D_grid_bio, n_grid_bio,vpresw_ice,vpresi_ice

       real*8, dimension(3) :: sig_array, the_array, frac_array
       real*8, dimension(:) :: ndust_ice, sigdust_ice,ddust_ice

       real :: n_iw, DSh_s , nbc_s, dbc_s, Asolo
       real, dimension (nbindust_ice) :: ndust_s, ddust_s
       integer :: index, kindex, mode,nbindust_ice

       si=six
       Si_=si+1d0
       SW=Si_*vpresi_ice/vpresw_ice
       Nd = zero_par
       Nbc = zero_par
       Norg = zero_par
       Nglassy = zero_par

       dNd = zero_par
       dNbc = zero_par
       dNorg = zero_par
       dNglassy = zero_par

       if ((six .lt. 0.02) .or. (T_ice .gt. 270.0)) then
       N=0d0
       Dsh=si
       return
       end if

       if (SW .ge. 1.0) then
       SW=1.0
       Si_=vpresw_ice/vpresi_ice
       si=Si_-1.0
       end if
       SIW=vpresw_ice/vpresi_ice

       sig_array = 0.0
       the_array = 0.0
       frac_array = 0.0


       select case (typeofspec_ice)

       case(1)
       N=1000d0*exp(-0.639d0)*(exp(12.96d0*si)-1d0)
       Dsh=1d0/12.96d0

       case(2)
       N=PDG07_ice(si, T_ice)
       if (T_ice .le. 243d0)then
       Dsh=1d0/3.88d0
       else
       Dsh=1d0/12.96d0
       end if

       case(3)

       Ndustaux=0.0d0
       DO index=1, nbindust_ice

       Ndustaux=Ndustaux+ndust_ice(index)
       end do


       if (si .le.shdust_ice) then
       Nd=(si/shdust_ice)*Ndustaux*effdust_ice* exp(-kdust_ice*
     &(shdust_ice-si))
       dNd=Nd*((1d0/si)+kdust_ice)

       else
       Nd=Ndustaux*effdust_ice
       dNd=0d0
       end if


       if (si .le.shbc_ice) then
       Nbc=(si/shbc_ice)*nbc_ice*effbc_ice* exp(-kbc_ice*(shbc_ice-si))
       dNbc=Nbc*((1d0/si)+kbc_ice)
       else
       Nbc=nbc_ice*effbc_ice
       dNbc=0d0
       end if

       N=Nd+Nbc
       if (((dNd+dNbc) .gt. 0d0) .and. (N .gt. 0.0)) then
       Dsh=N/(dNd+dNbc)
       else
       Dsh=0.0
       end if

       case(4)



       SW0=0.97d0
       delw0=cubicint_ice(SW, SW0, 1d0, 0d0, 1d0)
       ddelw0=dcubicint_ice(SW, SW0, 1d0, 0d0, 1d0)

       Nbase=PDG07_ice(si, T_ice)


       if (T_ice .le. 243d0)then
       dNbase=3.88d0*Nbase
       else
       dNbase=12.96d0*Nbase
       end if


       del0=cubicint_ice(Si_, si0dust_ice, si0dust_ice+0.1d0, 0d0, 1d0)
       ddel0=dcubicint_ice(Si_, si0dust_ice, si0dust_ice+0.1d0, 0d0,
     & 1d0)

       fc=0.5d0*del1dust_ice*del0
       dfc=0.5d0*del1dust_ice*ddel0

       Hdust=fc+((1d0-fc)*delw0)
       dH=(dfc*(1d0-delw0))+(ddelw0*(1d0-fc))

       if (Hdust .gt. 1d0) then
       Hdust=1d0
       dH=0d0
       end if


       aux=(2d0/3d0)*Hdust*(Nbase/0.76d0)*pi_ice/5.0d-7/4d0
       aux2=(2d0/3d0)*pi_ice/0.76d0/5.0d-7/4d0

       Nd=0d0
       dNd=0d0

       DO index =1, nbindust_ice


       Dx2=ddust_ice(index)*ddust_ice(index)*ddust_ice(index)*0.52*
     &acorr_dust

       frac=0.5d0*(1d0-erfapp(-log(ddust_ice(index)/0.1e-6) /
     &sigdust_ice(index)/sq2_par))

       Ndaux=frac*ndust_ice(index)*(1d0-exp(-aux*Dx2))

       Nd=Nd+Ndaux
       Ndaux=(frac*ndust_ice(index)-Ndaux)
       dNdaux=Ndaux*((dH*Nbase)+(Hdust*dNbase))*aux2*Dx2

       dNd=dNd+dNdaux

       END DO




       del0=cubicint_ice(Si_, si0bc_ice, si0bc_ice+0.1d0, 0d0, 1d0)
       ddel0=dcubicint_ice(Si_, si0bc_ice, si0bc_ice+0.1d0, 0d0, 1d0)

       fc=0.5d0*del1bc_ice*del0
       Hbc=fc+((1d0-fc)*delw0)
       dfc=0.5d0*del1bc_ice*ddel0
       dH=(dfc*(1d0-delw0))+(ddelw0*(1d0-fc))


       if (Hbc .gt. 1d0) then
       Hbc=1d0
       dH=0d0
       end if

       frac=0.5d0*(1d0 -erfapp(-log(dbc_ice/0.1e-6) /sigbc_ice/
     &sq2_par))


       Dx2=dbc_ice*dbc_ice*dbc_ice*0.52*acorr_bc

       aux=((1d0/3d0)-0.06d0)*Hbc*(Nbase/0.76d0)*pi_ice/2.7d-7
       aux2=((1d0/3d0)-0.06d0)*pi_ice/0.76d0/2.7d-7

       Nbc=nbc_ice*frac*(1d0-exp(-aux*Dx2))
       dNbc=(nbc_ice*frac-Nbc)*((dH*Nbase)+(Hbc*dNbase))*aux2*Dx2




       frac=0.5d0*(1d0-erfapp(-log(dorg_ice/0.1e-6) /sigorg_ice/
     &sq2_par))

       Dx2=dorg_ice*dorg_ice

       aux=0.06d0*Hbc*(Nbase/0.76d0)*pi_ice/9.1d-7
       aux2=0.06d0*pi_ice/0.76d0/9.1d-7


       Norg=norg_ice*frac*(1d0-exp(-aux*Dx2))
       dNorg=(norg_ice*frac-Norg)*((dH*Nbase)+(Hbc*dNbase))*aux2*Dx2




       N=Nd+Nbc+Norg

       fdust_dep = Nd
       Nhet_dep = N


       if (.FALSE.) then

       if (T_ice .lt. 210.0) then
       Nglassy= min(0.01 +0.0045*(210.0 -T_ice), 0.1)
       fglassy= 7.7211*1e-3*Si_ - 9.2688*1e-3
       fglassy=min(fglassy, 3.3587e-3)
       fglassy=max(fglassy, 0.0)
       else
       Nglassy = 0.0
       fglassy = 0.0
       end if
       Nglassy = np_ice*Nglassy*fglassy
       dNglassy =Nglassy*7.7211*1e-3
       N=N+Nglassy

       end if



       if ((dNd+dNbc+dNorg+dNglassy) .gt. 0d0) then
       Dsh=N/(dNd+dNbc+dNorg+dNglassy)
       else
       Dsh=0.0
       end if


       case (5)



       D_grid_bio =dorg_ice



       n_grid_bio = 0.0

       if (is_gocart) then
       ndust_s = SNGL(frac_dust*ndust_ice)
       nbc_s = SNGL(frac_bc*nbc_ice)
       Asolo = SNGL(0.25d0*frac_org*norg_ice*pi_par*dorg_ice*dorg_ice)
       else


       DO index =1, nbindust_ice
       frac=0.5d0*(1d0-erfapp(log(0.1e-6/ddust_ice(index)) /
     &sigdust_ice(index)/sq2_par))
       ndust_s(index) = SNGL(frac*ndust_ice(index))
       end do

       frac=0.5d0*(1d0-erfapp(log(0.1e-6/dbc_ice) /sigbc_ice/sq2_par))

       nbc_s = SNGL(frac*nbc_ice)

       frac=0.5d0*(1d0-erfapp(log(0.1e-6/dorg_ice) /sigorg_ice/
     &sq2_par))*0.25d0

       Asolo = SNGL(frac*norg_ice*pi_par*dorg_ice*dorg_ice)
       end if

       ddust_s=SNGL(ddust_ice)
       dbc_s = SNGL(dbc_ice*1.0)

       call EMPIRICAL_PARAM_PHILLIPS(SNGL(Si_), SNGL(SIW), SNGL(SW), (/
     &ddust_s/), (/ndust_s/), 5, (/dbc_s/), (/nbc_s/), 1, (/
     &SNGL(D_grid_bio)/), (/SNGL(n_grid_bio)/), 1, Asolo, n_iw, DSh_s
     & ,Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice)

       N=DBLE(n_iw)
       DSh=DBLE(DSh_s)

       case default
       N=zero_par
       Dsh=0.0
       end select

       if (Dsh .ge. si) then
       Dsh=si
       end if

       if (T_ice .gt. 255.0) then
       N=zero_par
       Dsh=zero_par
       end if

       end subroutine INSPEC_ice


!*************************************************************
!  Subroutine INimmersion
!  Provides the Immersion IN (for activated droplets) concentration at given T(K) according to Barahona et al. GMD, 2014.
!      Written by Donifan Barahona
!      donifan.o.barahona@nasa.gov
!========================================

       subroutine INimmersion(INconc, dINconc, wparcel,dbc_ice,sigbc_ice
     &  ,nbc_ice,fdust_imm,fdrop_dust,fdrop_bc,ndust_ice, sigdust_ice,
     &   ddust_ice,nbindust_ice,vpresw_ice,vpresi_ice,T_ice)

       real*8, intent (in) :: wparcel,dbc_ice,sigbc_ice,nbc_ice
     &  ,fdrop_dust,fdrop_bc
       real*8, intent (out) :: INconc, dINconc,fdust_imm
       real*8 :: Nd, Nd_unc, Nd_coa, Nbc, ahet, frac, dfrac, Naux,dNaux,
     & Tx, nssoot, nsdust, ninbkg, SX , dnsd, dnss, dNd, dNbc, coolr,
     & min_ns_dust, min_ns_soot, dninbkg,vpresw_ice,vpresi_ice,T_ice

       real*8, dimension(3) :: sig_array, the_array, frac_array
       real*8, dimension(:) :: ndust_ice, sigdust_ice,ddust_ice
       integer :: index, kindex, nbindust_ice

       if (T_ice .lt. Thom) then
       INconc = zero_par
       dINconc = zero_par
       return
       end if

       if (T_ice .ge. To_ice) then
       INconc = zero_par
       dINconc = zero_par
       return
       end if

       Tx=T_ice -273.16
       coolr=5.0e-3*wparcel
       min_ns_dust= 3.75e6
       min_ns_soot= 3.75e9

       SX=(vpresw_ice/vpresi_ice)-1.0


       ninbkg = 0.0
       dninbkg=0.0





       nsdust= max(exp(-0.517*Tx + 8.934)-min_ns_dust, 0.0)
       dnsd = max(0.517*nsdust, 0.0)

       nssoot= max(1.0e4*exp(-0.0101*Tx*Tx - 0.8525*Tx + 0.7667)-
     &min_ns_soot, 0.0)
       dnss = max(-(-2.0*0.0101*Tx -0.8525)*nssoot, 0.0)

       Naux=zero_par
       dNaux=zero_par


       if (is_gocart) then

       DO index =1,nbindust_ice

       ahet= ahet_dust(index)
       Naux=(1.0-exp(-nsdust*ahet))*ndust_ice(index)+Naux
       dNaux = exp(-nsdust*ahet)*ndust_ice(index)*dnsd*coolr*ahet+dNaux
       END DO
       ahet = ahet_bc
       else

       DO index =1,nbindust_ice

       ahet= ddust_ice(index)*ddust_ice(index)*ddust_ice(index)*0.52*
     &acorr_dust* exp(4.5*sigdust_ice(index)*sigdust_ice(index))
       Naux=(1.0-exp(-nsdust*ahet))*ndust_ice(index)+Naux
       dNaux = exp(-nsdust*ahet)*ndust_ice(index)*dnsd*coolr*ahet+dNaux
       END DO

       ahet=dbc_ice*dbc_ice*dbc_ice*0.52*acorr_bc* exp(4.5*sigbc_ice*
     &sigbc_ice)
       end if


       Nd=Naux*fdrop_dust
       dNd=dNaux*fdrop_dust

       Nbc=nbc_ice*(1.0-exp(-nssoot*ahet))*fdrop_bc
       dNbc= nbc_ice*exp(-nssoot*ahet)*fdrop_bc*dnss*coolr*ahet




       INconc=Nbc+ Nd + ninbkg
       dINconc=dNbc+dNd + dninbkg
       fdust_imm = Nd


       end subroutine INimmersion

!
!=======================================================================================
!=======================================================================================
!=======================================================================================
!!====================================================================================
!            EMPIRICAL PARAMETERISATION (Phillips et al. 2013, JAS)
!            Code contributed by Vaughan Phillips, University of Leeds
! Implementation:   Donifan Barahona donifan.o.barahona@nasa.gov
!====================================================================================


      SUBROUTINE EMPIRICAL_PARAM_PHILLIPS(SI, SIW, SW, D_grid_dust,
     & n_grid_dust, ijstop_dust, D_grid_soot, n_grid_soot, ijstop_soot,
     & D_grid_bio, n_grid_bio, ijstop_bio, A_solo, n_iw, DSH,
     & Nhet_dep,Nhet_dhf,fdust_dep,P_ice, T_ice)
      implicit none
      real, intent(IN):: SI, SIW, SW, A_solo,P_ice, T_ice
      real, dimension(:), intent(IN):: D_grid_dust, n_grid_dust,
     & D_grid_soot, n_grid_soot, D_grid_bio, n_grid_bio
      integer, intent(IN):: ijstop_dust, ijstop_soot, ijstop_bio
!     integer ijstop_dust, ijstop_soot, ijstop_bio

      real :: nin_a_nuc_dust, nin_a_nuc_soot, nin_a_nuc_bio,
     & nin_a_nuc_solo, num_ic_dust_imm, num_ic_soot_imm, num_ic_bio_imm,
     & num_ic_solo_imm

      real, intent (inout) :: DSH, n_iw
      real, intent (out) :: Nhet_dep,Nhet_dhf,fdust_dep

      real :: dn_in_dust, dn_in_soot, dn_in_bio, dn_in_solo, dNall,
     & dNaux, naux, SS_w, dH_frac_dust, dH_frac_soot, dH_frac_solo, aux,
     & dfdep, temperature_K, P_SAT, ahet_aux



      REAL :: RHO_CFDC, BASE_DUST_OMEGA, BASE_SOOT_PHILIC_OMEGA,
     & BASE_BIO_OMEGA, ALPHA_DUST, ALPHA_SOOT, ALPHA_bio,
     & FRACTION_DEPNUCL_WARM_DUST, PIE, BASE_SOLO_OMEGA,
     & TEMP_MAX_DUST_DEGC, TEMP_MAX_SOOT_DEGC, TEMP_MAX_bio_DEGC,
     & GLASS_FRAC
      PARAMETER( BASE_DUST_OMEGA = 2.0e-6, BASE_SOOT_PHILIC_OMEGA = 1.e-
     &7, BASE_BIO_OMEGA = 0.89e-6, BASE_SOLO_OMEGA = 5.6e-5,
     & GLASS_FRAC = 0.5, ALPHA_DUST = 2./3., ALPHA_SOOT = 1./3. - 0.03,
     & ALPHA_bio = 0.03, RHO_CFDC = 50000./(287.*228.15),
     & FRACTION_DEPNUCL_WARM_DUST = 0.15, PIE = 3.1415926,
     & TEMP_MAX_DUST_DEGC = -10., TEMP_MAX_SOOT_DEGC = -15.,
     & TEMP_MAX_bio_DEGC = -2.)

      real :: FAC_CORRECT_RH = 2., rho_AIDA
      real:: H_frac_dust, n_in, n_in_dust, n_in_ultra, n_in_dust_ultra,
     & CIHENC_dust, ESW, ESI, SS_i, n_in_soot_ultra, H_frac_soot,
     & H_frac_bio, n_in_soot, n_in_bio, n_in_bio_ultra, CIHENC_soot,
     & CIHENC_bio, delta_Si, delta_T, delta_Sw, n_in_max, SS_iw, rho

      real :: H_frac_solO, RHI, n_in_solO, n_in_solO_star, CIHENC_solO,
     & Psi_solO

      real :: mu, S_i_0, RH_crit, S_i_w_warm, S_i_w_cold, S_i_w,
     & tc_HM_degC
      real :: S_w_0, dep_frac, n_in_hat, n_in_tilde
      real*4 :: dH1smooth
      real :: EPS = 0.622
      integer :: ij
!intrinsic :: exp, DEXP, SIZE, DBLE



!print *, SIZE(n_grid_dust(:))
      if(ijstop_dust .ne. SIZE(n_grid_dust)) stop 6366
      if(ijstop_soot .ne. SIZE(n_grid_soot)) stop 6366
      if(ijstop_bio .ne. SIZE(n_grid_bio)) stop 6366
!     ijstop_dust = SIZE(n_grid_dust)
!     ijstop_soot = SIZE(n_grid_soot)
!     ijstop_bio = SIZE(n_grid_bio)



!Naux=12.96 !default
      naux =0.0
!     Nall =dNaux
!     Sh=0.0
      n_iw=0.0
      nin_a_nuc_dust=0.0; nin_a_nuc_soot=0.0; nin_a_nuc_bio=
     &0.0; nin_a_nuc_solo=0.0
      num_ic_dust_imm=0.0; num_ic_soot_imm=0.0; num_ic_bio_imm=
     &0.0; num_ic_solo_imm=0.0
!     AnningC uncom following four lines
      n_in_dust=0.0; dn_in_dust=0.0;n_in_soot=0.0
      dn_in_soot=0.0; dn_in_bio=0.0; n_in_bio=0.0;
      dn_in_solO=0.0;n_in_solO=0.0;n_in_dust_ultra=0.0
      n_in_soot_ultra=0.0;n_in_bio_ultra=0.0;


      H_frac_dust = 0.0
      H_frac_soot = 0.0
      H_frac_solo = 0.0
!     H1smooth=0.0
      aux=0.0

      temperature_K=SNGL(T_ice)
      P_SAT =SNGL(P_ice)

!A_solo = 1e-7 !m2 kg-1


!====================================================================================
!            COMPUTATION BLOCK
!
!====================================================================================
!
       rho_AIDA = 90000./(287.*205.)

       rho = P_SAT/(287.*temperature_K)

       Psi_solO = A_solO/BASE_SOLO_OMEGA
       SS_i = min(max(SI-1.0, 0.0), 1.0)
       SS_w = min(max(SW-1.0, -1.0), 1.0)
       SS_iw = min(max(SIW - 1.0, 0.0), 1.0)


       if(SS_i > 0.0) then
       if(temperature_K < 273.15 .and. temperature_K > 273.15 -
     & 90. ) then


       if(SS_w > 0.) then
       SS_i = SS_iw
       SS_w = 0.0
       end if


!                        S_i_zero = 1.15 !this is taken care of


       delta_Si = H_1_smooth(SS_i + 1, 1.1, 1.2, 0.0, 1.,dH1smooth);
       delta_T = H_1_smooth(-(temperature_K-273.15), 35., 40.,
     & FRACTION_DEPNUCL_WARM_DUST, 1.,dH1smooth);
       delta_Sw = H_1_smooth(SS_w + 1.0, 0.97, 1., 0., 1.,dH1smooth);

       tc_HM_degC = temperature_K - 273.15


       S_i_0 = 1. + 10.**(8.2584e-6*tc_HM_degC*tc_HM_degC*tc_HM_degC +
     & 5.3938E-4*tc_HM_degC*tc_HM_degC + 3.1656E-3*tc_HM_degC - 1.0261)


       S_w_0 = 0.97

       aux =H_1_smooth(-(temperature_K-273.15), 35., 40.,
     & FRACTION_DEPNUCL_WARM_DUST, 1.,dH1smooth)/FAC_CORRECT_RH

       dep_frac = H_1_smooth(SS_i + 1, S_i_0, S_i_0 + 0.1, 0.,1.,
     &  dH1smooth)* aux
       dfdep=dH1smooth*aux

       aux= H_1_smooth(SS_w + 1.0, S_w_0, 1., 0.,1.,dH1smooth)

       H_frac_dust = dep_frac + (1. - dep_frac)*aux

       dH_frac_dust = dfdep + (SIW*(1. - dep_frac)*dH1smooth)- aux*
     &dfdep

       if(H_frac_dust > 1.) H_frac_dust = 1.

       if ((H_frac_dust .gt. 1.0e-6) .and. (H_frac_dust .lt. 1.)) then
       dH_frac_dust = dH_frac_dust/H_frac_dust
       else
       dH_frac_dust =0.0
       end if


       S_i_0 = 1.2

       aux =H_1_smooth(-(temperature_K-273.15), 65., 75., 0.,1.,
     &   dH1smooth)
       dep_frac = H_1_smooth(SS_i + 1, S_i_0, S_i_0+0.1, 0.,1.,
     &    dH1smooth)*aux
       H_frac_solO = dep_frac
       dH_frac_solo=0.0
       if ((H_frac_solo .gt. 1.0e-6) .and. (H_frac_solo .lt. 1.)) then
       dH_frac_solo = dH1smooth/H_frac_solo
       end if
       if(H_frac_solO > 1.) H_frac_solO = 1.


       S_w_0 = 0.97

       S_i_0 = 1.3
!
       aux = H_1_smooth(-(temperature_K-273.15), 40., 50., 0.,1.,
     &   dH1smooth) / FAC_CORRECT_RH
       dep_frac = H_1_smooth(SS_i + 1, S_i_0, S_i_0+0.1, 0.,1.,
     &   dH1smooth)* aux

       dfdep= dH1smooth*aux

       aux = H_1_smooth(SS_w + 1.0, S_w_0, 1., 0.,1.,dH1smooth)
       H_frac_soot = dep_frac + (1. - dep_frac)*aux
       if(H_frac_soot > 1.) H_frac_soot = 1.

       dH_frac_soot = dfdep + (SIW*(1. - dep_frac)*dH1smooth)- aux*
     &dfdep
       if ((H_frac_soot .gt. 1.0e-6) .and. (H_frac_soot .lt. 1.)) then
       dH_frac_soot = dH_frac_soot/H_frac_soot
       else
       dH_frac_soot =0.0
       end if


       H_frac_bio = H_frac_soot

       if(temperature_K < 273.15 .and. temperature_K >= 273.15 -
     & 35.) then
       n_in = 1.E3* (exp(12.96*SS_i - 0.639)/RHO_CFDC) *0.0587*
     &FAC_CORRECT_RH
       if( temperature_K > 273.15 -5. .and. temperature_K < 273.15 -
     & 2. ) then
       n_in = n_in*H_1_smooth(-(temperature_K-273.15), 2., 5., 0., 1.,
     &   dH1smooth)
       endif
       if(temperature_K >= 273.15 - 2. ) n_in = 0.


       if(temperature_K < 273.15 -25. ) then
       n_in_tilde = 1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)*
     &FAC_CORRECT_RH/RHO_CFDC
       n_in_hat = n_in

       if(temperature_K >= 273.15 - 30.) n_in_max = 1.E3* (exp(12.96*
     &SS_iw - 0.639)/RHO_CFDC) *0.0587*FAC_CORRECT_RH
       if(temperature_K < 273.15 - 30.) n_in_max = 1000.*(exp(0.1296*
     &(SS_iw*100.-10.))**0.3)*FAC_CORRECT_RH/RHO_CFDC

       if(n_in_hat > n_in_max) n_in_hat = n_in_max
       if(n_in_tilde > n_in_max) n_in_tilde = n_in_max



       n_in = n_in_hat * ((n_in_tilde/n_in_hat)**(H_1_smooth(-
     &(temperature_K-273.15), 25., 35., 0., 1.,dH1smooth)))


       if(n_in > n_in_max) n_in = n_in_max

       endif
       n_in_dust = 0.
       dn_in_dust = 0.


       if(temperature_K < 273.15 - 30.) then
       dnaux = 3.88
       else
       dnaux = 12.96
       end if


       naux=0.0
       do ij = 1, ijstop_dust
       if (is_gocart) then
       mu = n_in*ALPHA_DUST*H_frac_dust*ahet_dust(ij)/BASE_DUST_OMEGA
       else
       ahet_aux = PIE*D_grid_dust(ij)*D_grid_dust(ij)*D_grid_dust(ij)*
     &4.73*acorr_dust/6.0
       mu = n_in*ALPHA_DUST*H_frac_dust*ahet_aux/BASE_DUST_OMEGA
       end if
       naux = (1. - exp(-mu))*n_grid_dust(ij)
       n_in_dust = n_in_dust + naux
       dn_in_dust = max(mu*(n_grid_dust(ij)-naux)*(dnaux +
     & dH_frac_dust), 0.0) + dn_in_dust

       enddo

       if( temperature_K > 273.15 +TEMP_MAX_DUST_DEGC -
     & 20. .and. temperature_K < 273.15 + TEMP_MAX_DUST_DEGC) then
       n_in_dust = n_in_dust*H_1_smooth(-(temperature_K-273.15),-
     &TEMP_MAX_DUST_DEGC,-TEMP_MAX_DUST_DEGC+20., 0., 1.,dH1smooth)
       endif
       if(temperature_K >= 273.15 + TEMP_MAX_DUST_DEGC) n_in_dust = 0.


       n_in_soot = 0.
       dn_in_soot = 0.
       do ij = 1, ijstop_soot

       if (is_gocart) then
       mu = n_in*ALPHA_SOOT*H_frac_soot*ahet_bc/BASE_SOOT_PHILIC_OMEGA
       else
       ahet_aux = PIE*D_grid_soot(ij)*D_grid_soot(ij)*D_grid_soot(ij)*
     &16.4*acorr_bc/6.0
       mu = n_in*ALPHA_SOOT*H_frac_soot*ahet_aux/BASE_SOOT_PHILIC_OMEGA
       end if
       naux = (1. - exp(-mu))*n_grid_soot(ij)
       n_in_soot = n_in_soot + naux
       dn_in_soot = max(mu*(n_grid_soot(ij)-naux)*(dnaux+dH_frac_soot),
     & 0.0) + dn_in_soot

       enddo

       if( temperature_K > 273.15 + TEMP_MAX_SOOT_DEGC -
     & 10. .and. temperature_K < 273.15 + TEMP_MAX_SOOT_DEGC) then
       n_in_soot = n_in_soot*H_1_smooth(-(temperature_K-273.15),-
     &TEMP_MAX_SOOT_DEGC,-TEMP_MAX_SOOT_DEGC+10., 0., 1.,dH1smooth)

       endif
       if(temperature_K >= 273.15 + TEMP_MAX_SOOT_DEGC) n_in_soot = 0.

       n_in_bio = 0.
       dn_in_bio = 0.


       do ij = 1, ijstop_bio
       mu = n_in*ALPHA_bio*H_frac_bio*PIE*(D_grid_bio(ij)**2.) /
     &BASE_BIO_OMEGA

       mu = n_in*ALPHA_bio*H_frac_bio
       naux = (1. - exp(-mu))*n_grid_bio(ij)




       n_in_bio = n_in_bio + naux
       dn_in_bio = max(mu*(n_grid_bio(ij)-naux)*dnaux, 0.0) + dn_in_bio


       enddo


       if( temperature_K > 273.15 + TEMP_MAX_bio_DEGC -
     & 3. .and. temperature_K < 273.15 + TEMP_MAX_bio_DEGC) then
       n_in_bio = n_in_bio*H_1_smooth(-(temperature_K-273.15),-
     &TEMP_MAX_bio_DEGC,-TEMP_MAX_bio_DEGC+3., 0., 1.,dH1smooth)

       endif
       if(temperature_K >= 273.15 + TEMP_MAX_bio_DEGC ) n_in_bio = 0.



       else
       n_in = 0.; n_in_ultra = 0.; n_in_dust = 0.; n_in_soot =
     & 0.; n_in_bio = 0.;
       endif

       if(temperature_K < 273.15 - 35.) then
       n_in_ultra = 1000.*(exp(0.1296*(SS_i*100.-10.))**0.3)*
     &FAC_CORRECT_RH/RHO_CFDC
       dnaux = 3.88
       naux=0.0


       RHI = (SS_i+1.)*100.
       if(RHI < 0.) RHI = 0.
       n_in_solO_star = 1000.e6*(7.7211e-5 * RHI - 9.2688e-3)/rho_AIDA


       n_in_dust_ultra = 0.;
       dn_in_dust = 0.0
       do ij = 1, ijstop_dust

       if (is_gocart) then
       mu = n_in_ultra*ALPHA_DUST*H_frac_dust*ahet_dust(ij)/
     &BASE_DUST_OMEGA
       else
       ahet_aux = PIE*D_grid_dust(ij)*D_grid_dust(ij)*D_grid_dust(ij)*
     &4.73*acorr_dust/6.0
       mu = n_in_ultra*ALPHA_DUST*H_frac_dust*ahet_aux/BASE_DUST_OMEGA
       end if



       naux = (1. - exp(-mu))*n_grid_dust(ij)
       n_in_dust_ultra = n_in_dust_ultra + naux
       dn_in_dust = max(mu*(n_grid_dust(ij)-naux)*(dnaux +dH_frac_dust),
     & 0.0) + dn_in_dust



       enddo



       n_in_soot_ultra = 0.0
       dn_in_soot = 0.0
       do ij = 1, ijstop_soot

       if (is_gocart) then
       mu = n_in_ultra*ALPHA_SOOT*H_frac_soot*ahet_bc/
     &BASE_SOOT_PHILIC_OMEGA
       else
       ahet_aux = PIE*D_grid_soot(ij)*D_grid_soot(ij)*D_grid_soot(ij)*
     &16.4*acorr_bc/6.0
       mu = n_in_ultra*ALPHA_SOOT*H_frac_soot*ahet_aux/
     &BASE_SOOT_PHILIC_OMEGA
       end if



       naux = (1. - exp(-mu))*n_grid_soot(ij)
       n_in_soot_ultra = n_in_soot_ultra + naux
       dn_in_soot = max(mu*(n_grid_soot(ij)-naux)*(dnaux +dH_frac_soot),
     & 0.0) + dn_in_soot

       enddo


       n_in_bio_ultra = 0.
       dn_in_bio = 0.0










       n_in_solO = Psi_solO*glass_frac*H_frac_solO*n_in_solO_star
       dn_in_solO =max(Psi_solO*glass_frac* (H_frac_solO*77211.0*100.0/
     &rho_AIDA + n_in_solO_star*dH_frac_solo), 0.0)


       else
       n_in_ultra = 0.; n_in_dust_ultra = 0.; n_in_soot_ultra =
     & 0.; n_in_solO = 0.; n_in_bio_ultra = 0.;
       endif



       n_in_dust = n_in_dust + n_in_dust_ultra;
       n_in_soot = n_in_soot + n_in_soot_ultra;
       n_in_bio = n_in_bio + n_in_bio_ultra;




! PROBLEM:  how to ensure that the frozen fraction does not exceed 1 ?

       if (.false.) then
       if(n_in_dust + n_in_bio + n_in_soot + n_in_solO > 0.) then

       CIHENC_dust = n_in_dust - nin_a_nuc_dust
       if(CIHENC_dust < 0.) CIHENC_dust = 0.

       CIHENC_soot = n_in_soot - nin_a_nuc_soot
       if(CIHENC_soot < 0.) CIHENC_soot = 0.

       CIHENC_bio = n_in_bio - nin_a_nuc_bio
       if(CIHENC_bio < 0.) CIHENC_bio = 0.

       CIHENC_solO = n_in_solO - nin_a_nuc_solO
       if(CIHENC_solO < 0.) CIHENC_solO = 0.


       n_iw = n_iw + CIHENC_dust
       nin_a_nuc_dust = nin_a_nuc_dust + CIHENC_dust
       num_ic_dust_imm = num_ic_dust_imm + CIHENC_dust

       n_iw = n_iw + CIHENC_soot
       nin_a_nuc_soot = nin_a_nuc_soot + CIHENC_soot
       num_ic_soot_imm = num_ic_soot_imm + CIHENC_soot

       n_iw = n_iw + CIHENC_bio
       nin_a_nuc_bio = nin_a_nuc_bio + CIHENC_bio
       num_ic_bio_imm = num_ic_bio_imm + CIHENC_bio

       n_iw = n_iw + CIHENC_solO
       nin_a_nuc_solO = nin_a_nuc_solO + CIHENC_solO
       num_ic_solO_imm = num_ic_solO_imm + CIHENC_solO


       endif
       end if
       endif
       endif

       n_iw = n_in_dust + n_in_bio + n_in_soot + n_in_solO
       dNall = dn_in_dust + dn_in_bio + dn_in_soot + dn_in_solO

       if (n_iw .gt. zero_par) then
       fdust_dep = DBLE(n_in_dust/n_iw)
       end if
       Nhet_dep = DBLE(n_in_dust)
       Nhet_dhf = DBLE(n_in_solO)




       if (( dNall > 0.) .and. (n_iw .gt. 0.0)) then
       Dsh=max(min(n_iw/dNall, SS_i), 0.005)
       else
       Dsh=0.005
       end if


      END SUBROUTINE EMPIRICAL_PARAM_PHILLIPS

      real function H_1(X, X_1, X_2, Hlo)
      real, intent(in) :: Hlo, X, X_1, X_2

      if(X >= X_2) H_1 = 1
      if(X <= X_1) H_1 = Hlo
      if(X > X_1 .and. X < X_2) H_1 = (X - X_1)/(X_2 - X_1)

      if( X_2 <= X_1) stop 91919

      return
      end function


      real function H_1_smooth(X, X_1, X_2, Hlo, Hhi,dH1smooth)
      real, intent(in) :: Hlo, Hhi, X, X_1, X_2
      real*4,intent(out) :: dH1smooth
      real :: a_0, a_1, a_2, a_3, A, B

      if(X >= X_2) H_1_smooth = Hhi
      if(X <= X_1) H_1_smooth = Hlo

      if(X >= X_2) dH1smooth = 0.0
      if(X <= X_1) dH1smooth = 0.0

      if(X > X_1 .and. X < X_2) then
       A = 6.*(Hlo - Hhi)/(X_2**3. - X_1**3. + 3.*(X_2*X_1*X_1 - X_1*
     &X_2*X_2) )
       a_3 = (A/3.)
       a_2 = -(A/2.)*(X_1 + X_2)
       a_1 = A*X_2*X_1
       B = Hlo + A*(X_1**3.)/6. - A*X_1*X_1*X_2/2.
       a_0 = B
       H_1_smooth = a_0 + a_1*X + a_2*X*X + a_3*X*X*X
       dH1smooth = a_1 + 2.0*a_2*X + 3.0*a_3*X*X
      endif

!H1smooth =min(dH1smooth , 1000000.0)
!H1smooth =max(dH1smooth , 0.0)

      if( X_2 <= X_1) stop 91919

      return
      end function





!    END ICE PARAMETERIZATION DONIF
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



       END MODULE aer_cloud






