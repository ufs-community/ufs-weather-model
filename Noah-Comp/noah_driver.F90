module noah_driver

  use machine, only: kind_phys

  use noah_loop, only: noah_loop_init, noah_loop_run
  use noah_type_mod



  implicit none
  private

  type (noah_type), public        ::   noah_pubinst

  public :: noah_loop_drv, init_driver

contains

  subroutine init_driver(procbounds)

    use proc_bounds,        only: procbounds_type
    use mpp_domains_mod,    only: domain2d, mpp_get_compute_domain
    use mpp_mod,            only: mpp_pe, mpp_root_pe
    use land_domain_mod,    only: domain_create
    use block_control_mod,  only: block_control_type, define_blocks_packed

    
    type(procbounds_type),  intent(in)    :: procbounds
    
    ! local
    !type (noah_type)        ::   noah
    type(domain2D) :: land_domain
    type (block_control_type), target   :: Lnd_block !  Block container
    integer :: isc, iec, jsc, jec
    
    integer                 ::   im         ! horiz dimension
    integer                 ::   isot       ! sfc soil type data source
    integer                 ::   ivegsrc    ! sfc veg type data source
    
    integer                 ::   blocksize
    logical, save        :: block_message = .true.
    
    character(len=128)      ::   errmsg     ! error messaging added to ccpp
    integer                 ::   errflg     ! error messaging added to ccpp

    integer                 :: i, j, nb, ix
    im = procbounds%im

    ! quick test, should be read in from nml
    isot = 1
    ivegsrc = 1
    blocksize = 64

    ! domain create with FMS:
    call domain_create(land_domain)

    ! Creat blocking a la FV3
    call mpp_get_compute_domain(land_domain,isc,iec,jsc,jec)
    write(*,*) "isc,iec,jsc,jec: ",isc,iec,jsc,jec
    
    call define_blocks_packed('land_model', Lnd_block, isc, iec, jsc, jec, 1, &
         blocksize, block_message)

    ! tmp debug
    if (mpp_pe() == mpp_root_pe()) then
       write(*,*) "block nblks, isc,iec,jsc,jec: ", Lnd_block%nblks, Lnd_block%isc, Lnd_block%iec, Lnd_block%jsc, Lnd_block%jec
       !write(*,*) "block blksz: ", Lnd_block%blksz

       do j=jsc,jec
          do i=isc,iec
             nb = Lnd_block%blkno(i,j)
             ix = Lnd_block%ixp(i,j)
             write(*,*) "i,j,nb,ix: ", i,j,nb,ix
          end do
       end do
    end if
    
    call noah_pubinst%Create(im)

    call noah_loop_init(0, isot, ivegsrc, 0 , errmsg, errflg)


  end subroutine init_driver

  subroutine noah_loop_drv(procbounds, noah_pubinst)

    use proc_bounds, only : procbounds_type

    use physcons, only :       &
         cp      => con_cp,    &
         eps     => con_eps,   &
         epsm1   => con_epsm1, &
         grav    => con_g,     &
         hvap    => con_hvap,  &
         rd      => con_rd,    &        
         rvrdm1  => con_fvirt, &
         tfreeze => con_t0c

    type(procbounds_type),  intent(in)    :: procbounds
    type(noah_type),        intent(inout)    :: noah_pubinst ! land model's variable type

    ! local
    integer                 :: i, de, gridbeg, gridend
    real(kind_phys)         :: foodata(procbounds%im)
    !
    real(kind_phys) :: ps        (procbounds%im) ! surface pressure (pa)                       im
    real(kind_phys) :: t1        (procbounds%im) ! surface layer mean temperature (k)          im
    real(kind_phys) :: q1        (procbounds%im) ! surface layer mean specific humidity        im
    integer         :: soiltyp   (procbounds%im) ! soil type (integer index)                   im
    integer         :: vegtype   (procbounds%im) ! vegetation type (integer index)             im
    real(kind_phys) :: sigmaf    (procbounds%im) ! areal fractional cover of green vegetation  im
    real(kind_phys) :: sfcemis   (procbounds%im) ! sfc lw emissivity ( fraction )              im
    real(kind_phys) :: dlwflx    (procbounds%im) ! total sky sfc downward lw flux ( w/m**2 )   im
    real(kind_phys) :: dswsfc    (procbounds%im) ! total sky sfc downward sw flux ( w/m**2 )   im
    real(kind_phys) :: snet      (procbounds%im) ! total sky sfc netsw flx into ground(w/m**2) im
    real(kind_phys) :: tg3       (procbounds%im) ! deep soil temperature (k)                   im
    real(kind_phys) :: cm        (procbounds%im) ! surface exchange coeff for momentum (m/s)   im
    real(kind_phys) :: ch        (procbounds%im) ! surface exchange coeff heat & moisture(m/s) im
    real(kind_phys) :: prsl1     (procbounds%im) ! sfc layer 1 mean pressure (pa)              im
    real(kind_phys) :: prslki    (procbounds%im) !                                             im
    real(kind_phys) :: zf        (procbounds%im) ! height of bottom layer (m)                  im
    logical         :: land      (procbounds%im) ! = T if a point with any land                im
    real(kind_phys) :: wind      (procbounds%im) ! wind speed (m/s)                            im
    integer         :: slopetyp  (procbounds%im) ! class of sfc slope (integer index)          im
    real(kind_phys) :: shdmin    (procbounds%im) ! min fractional coverage of green veg        im
    real(kind_phys) :: shdmax    (procbounds%im) ! max fractnl cover of green veg (not used)   im
    real(kind_phys) :: snoalb    (procbounds%im) ! upper bound on max albedo over deep snow    im
    real(kind_phys) :: sfalb     (procbounds%im) ! mean sfc diffused sw albedo (fractional)    im
    logical         :: flag_iter (procbounds%im) !                                             im
    logical         :: flag_guess(procbounds%im) !                                             im
    real(kind_phys) :: bexppert  (procbounds%im)
    real(kind_phys) :: xlaipert  (procbounds%im)
    real(kind_phys) :: vegfpert  (procbounds%im)
    real(kind_phys) :: weasd     (procbounds%im) ! water equivalent accumulated snow depth(mm) im
    real(kind_phys) :: snwdph    (procbounds%im) ! snow depth (water equiv) over land          im
    real(kind_phys) :: tskin     (procbounds%im) ! ground surface skin temperature ( k )       im
    real(kind_phys) :: tprcp     (procbounds%im) ! total precipitation                         im
    real(kind_phys) :: srflag    (procbounds%im) ! snow/rain flag for precipitation            im
    real(kind_phys) :: canopy    (procbounds%im) ! canopy moisture content (m)                 im
    real(kind_phys) :: trans     (procbounds%im) ! total plant transpiration (m/s)             im
    real(kind_phys) :: tsurf     (procbounds%im) ! surface skin temperature (after iteration)  im
    real(kind_phys) :: z0rl      (procbounds%im) ! surface roughness                           im
    real(kind_phys) :: sncovr1   (procbounds%im) ! snow cover over land (fractional)            im
    real(kind_phys) :: qsurf     (procbounds%im) ! specific humidity at sfc                     im
    real(kind_phys) :: gflux     (procbounds%im) ! soil heat flux (w/m**2)                      im
    real(kind_phys) :: drain     (procbounds%im) ! subsurface runoff (mm/s)                     im
    real(kind_phys) :: evap      (procbounds%im) ! evaperation from latent heat flux            im
    real(kind_phys) :: hflx      (procbounds%im) ! sensible heat flux                           im
    real(kind_phys) :: ep        (procbounds%im) ! potential evaporation                        im
    real(kind_phys) :: runoff    (procbounds%im) ! surface runoff (m/s)                         im
    real(kind_phys) :: cmm       (procbounds%im) !                                              im
    real(kind_phys) :: chh       (procbounds%im) !                                              im
    real(kind_phys) :: evbs      (procbounds%im) ! direct soil evaporation (m/s)                im
    real(kind_phys) :: evcw      (procbounds%im) ! canopy water evaporation (m/s)               im
    real(kind_phys) :: sbsno     (procbounds%im) ! sublimation/deposit from snopack (m/s)       im
    real(kind_phys) :: snowc     (procbounds%im) ! fractional snow cover                        im
    real(kind_phys) :: stm       (procbounds%im) ! total soil column moisture content (m)       im
    real(kind_phys) :: snohf     (procbounds%im) ! snow/freezing-rain latent heat flux (w/m**2) im
    real(kind_phys) :: smcwlt2   (procbounds%im) ! dry soil moisture threshold                  im
    real(kind_phys) :: smcref2   (procbounds%im) ! soil moisture threshold                      im
    real(kind_phys) :: wet1      (procbounds%im) ! normalized soil wetness                      im

    real(kind_phys) :: ustar     (procbounds%im)
    real(kind_phys) :: smc(procbounds%im,noah_pubinst%static%km) ! total soil moisture content (fractional)   im,km
    real(kind_phys) :: stc(procbounds%im,noah_pubinst%static%km) ! soil temp (k)                              im,km
    real(kind_phys) :: slc(procbounds%im,noah_pubinst%static%km) ! liquid soil moisture                       im,km

    ! tmp for testing. Missing imports for these
    real(kind_phys) :: prsik1(procbounds%im), z0pert(procbounds%im), &
                       ztpert(procbounds%im), stress(procbounds%im)

    ! tmp for testing. These should be coming from namelist
    real(kind_phys), parameter :: delt = 900.0_kind_phys
    integer, parameter :: ivegsrc        = 1
    integer, parameter :: isot           = 1
    logical, parameter :: lheatstrg      = .false.
    real(kind_phys), parameter :: pertvegf = 0.0_kind_phys
    ! outputs
    real(kind_phys) :: rb_lnd(procbounds%im), fm_lnd(procbounds%im), &
                       fh_lnd(procbounds%im), fm10_lnd(procbounds%im), fh2_lnd(procbounds%im)
    
    associate(foodata => noah_pubinst%model%foo_atm2lndfield,  &
         ps         => noah_pubinst%model%ps         ,&
         t1         => noah_pubinst%model%t1         ,&
         q1         => noah_pubinst%model%q1         ,&
         dlwflx     => noah_pubinst%model%dlwflx     ,&
         dswsfc     => noah_pubinst%model%dswsfc     ,&
         dswsfci    => noah_pubinst%model%dswsfci    ,&
         wind       => noah_pubinst%model%wind       ,&
         tprcp      => noah_pubinst%model%tprcp      ,&
         im         => procbounds%im                 ,& 
         km         => noah_pubinst%static%km        ,&
         ! delt       => noah_pubinst%static%delt      ,&
         ! isot       => noah_pubinst%static%isot      ,&
         ! ivegsrc    => noah_pubinst%static%ivegsrc   ,&
         ! pertvegf   => noah_pubinst%static%pertvegf  ,&  
         ! lheatstrg  => noah_pubinst%static%lheatstrg ,& 
         errmsg     => noah_pubinst%static%errmsg    ,& 
         errflg     => noah_pubinst%static%errflg    ,& 
         soiltyp    => noah_pubinst%model%soiltyp    ,&
         vegtype    => noah_pubinst%model%vegtype    ,&
         slopetyp   => noah_pubinst%model%slopetyp   ,&
         sigmaf     => noah_pubinst%model%sigmaf     ,&
         sfcemis    => noah_pubinst%model%sfcemis    ,&
         snet       => noah_pubinst%model%snet       ,&
         tg3        => noah_pubinst%model%tg3        ,&
         cm         => noah_pubinst%model%cm         ,&
         ch         => noah_pubinst%model%ch         ,&
         prsl1      => noah_pubinst%model%prsl1      ,&
         shdmin     => noah_pubinst%model%shdmin     ,&
         shdmax     => noah_pubinst%model%shdmax     ,&
         snoalb     => noah_pubinst%model%snoalb     ,&
         sfalb      => noah_pubinst%model%sfalb      ,&
         zf         => noah_pubinst%model%zf         ,&
         bexppert   => noah_pubinst%model%bexppert   ,&
         xlaipert   => noah_pubinst%model%xlaipert   ,&
         vegfpert   => noah_pubinst%model%vegfpert   ,&
         flag_iter  => noah_pubinst%model%flag_iter  ,&
         flag_guess => noah_pubinst%model%flag_guess ,&
         land       => noah_pubinst%model%land       ,&
         prslki     => noah_pubinst%model%prslki     ,&
         weasd      => noah_pubinst%model%weasd      ,&
         snwdph     => noah_pubinst%model%snwdph     ,&
         tskin      => noah_pubinst%model%tskin      ,&
         srflag     => noah_pubinst%model%srflag     ,&
         canopy     => noah_pubinst%model%canopy     ,&
         trans      => noah_pubinst%model%trans      ,&   
         tsurf      => noah_pubinst%model%tsurf      ,&
         z0rl       => noah_pubinst%model%z0rl       ,&
         ! smc        => noah_pubinst%model%smc        ,&
         ! stc        => noah_pubinst%model%stc        ,&
         ! slc        => noah_pubinst%model%slc        ,&
         sncovr1    => noah_pubinst%model%sncovr1    ,&
         qsurf      => noah_pubinst%model%qsurf      ,&
         gflux      => noah_pubinst%model%gflux      ,&
         drain      => noah_pubinst%model%drain      ,&
         evap       => noah_pubinst%model%evap       ,&
         hflx       => noah_pubinst%model%hflx       ,&
         ep         => noah_pubinst%model%ep         ,&
         runoff     => noah_pubinst%model%runoff     ,&
         cmm        => noah_pubinst%model%cmm        ,&
         chh        => noah_pubinst%model%chh        ,&
         evbs       => noah_pubinst%model%evbs       ,&
         evcw       => noah_pubinst%model%evcw       ,&
         sbsno      => noah_pubinst%model%sbsno      ,&
         snowc      => noah_pubinst%model%snowc      ,&
         stm        => noah_pubinst%model%stm        ,&
         snohf      => noah_pubinst%model%snohf      ,&
         smcwlt2    => noah_pubinst%model%smcwlt2    ,&
         smcref2    => noah_pubinst%model%smcref2    ,&
         wet1       => noah_pubinst%model%wet1       ,&
         !
         rb_lnd     => noah_pubinst%model%rb_lnd     ,&
         fm_lnd     => noah_pubinst%model%fm_lnd     ,&
         fh_lnd     => noah_pubinst%model%fh_lnd     ,&
         fm10_lnd   => noah_pubinst%model%fm10_lnd   ,&
         fh2_lnd    => noah_pubinst%model%fh2_lnd    ,&
         stress     => noah_pubinst%model%stress     ,&
         !
         ustar      => noah_pubinst%model%ustar       &
         )

      ! tmp for testing. Missing imports for these
      prsik1    = 0.1_kind_phys
      z0pert    = 0.0_kind_phys
      ztpert    = 0.0_kind_phys
      !ustar     = 0.1_kind_phys
      stress    = 0.1_kind_phys  ! think stress is not needed as input, only output from stability
      smc       = 0.65_kind_phys
      slc       = 0.65_kind_phys
      stc       = 296_kind_phys
      
      ! first test
      !write(*,*) 'NLP test: ', noah_pubinst%model%soiltyp ! get all zeros, good
      
      de = procbounds%de
      !im = procbounds%im
      gridbeg = procbounds%gridbeg
      gridend = procbounds%gridend

      ! tmp, debug
      !write(6,'("noah drv: dswsfc   - min/max/avg",3g16.6)') minval(dswsfc),   maxval(dswsfc),   sum(dswsfc)/size(dswsfc)

      call noah_loop_run( &
                                !! ARGS FROM NOAH
                                !  ---  inputs:
       im, km, grav, cp, hvap, rd, eps, epsm1, rvrdm1, ps,          &
       t1, q1, soiltyp, vegtype, sigmaf,                            &
       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,            &
       prsl1, prslki, zf, land, wind, slopetyp,                     &
       shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,        &
       lheatstrg, isot, ivegsrc,                                    &
       bexppert, xlaipert, vegfpert,pertvegf,                       & ! sfc perts, mgehne
                                !     ---  in/outs:
       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,          &
       canopy, trans, tsurf, z0rl,                                  &
                                !     ---  outputs:
       sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,        &
       cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,              &
       smcwlt2, smcref2, wet1,                                      &
                                !! ARGS FROM  stab_prep_lnd (minus those from noah
                                !  ---  inputs:
       prsik1,z0pert,ztpert,ustar,                                  &
                                !  ---  outputs:
                                !! ARGS FROM stability (minus those above)
                                !  ---  inputs:
                                !  ---  outputs:
       rb_lnd, fm_lnd, fh_lnd, fm10_lnd, fh2_lnd,                   &
       stress,                                                      &
                                !!
       errmsg, errflg                                               &
       )


      ! foodata(1:gridend-gridbeg+1) = noah_pubinst%model%foo_atm2lndfield(gridbeg:gridend)
      ! foodata = noah_pubinst%model%foo_atm2lndfield(gridbeg:gridend)
      ! foodata = noah_pubinst%model%foo_atm2lndfield


      ! do i = 1,im
      !    write(*,*) 'NLP2: ', de, gridbeg, gridend, foodata(i)
      ! end do

    end associate

  end subroutine noah_loop_drv

end module noah_driver
