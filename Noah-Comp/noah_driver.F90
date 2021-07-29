module noah_driver

  use machine, only: kind_phys

  use noah_loop, only: noah_loop_init, noah_loop_run
  use noah_type_mod
  use proc_bounds,        only: procbounds_type, control_init_type
  


  implicit none
  private

  type (noah_type),         public        ::   noah_model
  !type (noah_type),  dimension(:), allocatable, public   ::   noah_model
  type (control_init_type),                     public   ::   ctrl_init

  public :: noah_loop_drv, init_driver, noah_block_run

contains

  !subroutine init_driver(procbounds)
  subroutine init_driver(ctrl_init)

    use mpp_domains_mod,    only: domain2d, mpp_get_compute_domain
    use mpp_mod,            only: mpp_pe, mpp_root_pe
    use land_domain_mod,    only: domain_create
    use block_control_mod,  only: block_control_type, define_blocks_packed

    
    !type(procbounds_type),  intent(in)    :: procbounds
    !character(len=*), intent(out) :: gridchoice
    type(control_init_type), intent(out)  ::   ctrl_init
    
    ! ---------------
    ! local

    !type(control_init_type)  ::   ctrl_init
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
    integer                 :: nblks
    !im = procbounds%im

    ! ! setup ctrl_init
    
    call ctrl_init%init()

    if (mpp_pe() == mpp_root_pe()) then
       write(*,*) 'ctrl_init%grid: '     ,ctrl_init%grid
       write(*,*) 'ctrl_init%npx: '      ,ctrl_init%npx
       write(*,*) 'ctrl_init%npy: '      ,ctrl_init%npy
       write(*,*) 'ctrl_init%layout: '   ,ctrl_init%layout
       write(*,*) 'ctrl_init%ntiles: '   ,ctrl_init%ntiles
       write(*,*) 'ctrl_init%blocksize: ',ctrl_init%blocksize
       write(*,*) 'ctrl_init%ivegsrc: '  ,ctrl_init%ivegsrc
       write(*,*) 'ctrl_init%isot: '     ,ctrl_init%isot
    end if

    !gridchoice = ctrl_init%grid
    ! domain create with FMS:
    call domain_create(ctrl_init, land_domain)

    ! Creat blocking a la FV3
    call mpp_get_compute_domain(land_domain,isc,iec,jsc,jec)

    im = (iec-isc+1)*(jec-jsc+1)
    write(*,*) "isc,iec,jsc,jec, im: ",isc,iec,jsc,jec, im
    
    
    ! Create blocks, but not curretnly using
    call define_blocks_packed('land_model', Lnd_block, isc, iec, jsc, jec, 1, &
         ctrl_init%blocksize, block_message)

    ! tmp debug
    if (mpp_pe() == mpp_root_pe()) then
       write(*,*) "block nblks, isc,iec,jsc,jec: ", Lnd_block%nblks, Lnd_block%isc, Lnd_block%iec, Lnd_block%jsc, Lnd_block%jec
       !write(*,*) "block blksz: ", Lnd_block%blksz

    !    do j=jsc,jec
    !       do i=isc,iec
    !          nb = Lnd_block%blkno(i,j)
    !          ix = Lnd_block%ixp(i,j)
    !          write(*,*) "i,j,nb,ix: ", i,j,nb,ix
    !       end do
    !    end do
    end if

    !! IF were using blocking:
    !nblks = size(Lnd_block%blksz)
    !allocate(noah_model(nblks))
    ! do nb = 1,nblks
    !    ix = Lnd_block%blksz(nb)
    !    call noah_model(nb)%Create(ix)

    !    noah_model(nb)%control%nblks = Lnd_block%nblks
    !    noah_model(nb)%control%blksz = Lnd_block%blksz(nb)
    !    noah_model(nb)%control%isc = Lnd_block%isc
    !    noah_model(nb)%control%iec = Lnd_block%iec
    !    noah_model(nb)%control%jsc = Lnd_block%jsc
    !    noah_model(nb)%control%jec = Lnd_block%jec
    ! end do
    !! ELSE:
    noah_model%control%isc = isc
    noah_model%control%iec = iec
    noah_model%control%jsc = jsc
    noah_model%control%jec = jec
    noah_model%static%im  = im
    call noah_model%Create(im)
    
    call noah_loop_init(0, ctrl_init%isot, ctrl_init%ivegsrc, 0 , errmsg, errflg)


  end subroutine init_driver

  subroutine noah_block_run(procbounds, noah_model)

    type(procbounds_type),  intent(in)    :: procbounds
    type(noah_type),        intent(inout) :: noah_model(:) ! land model's variable type


    integer                 :: nb
    do nb=1, size(noah_model)
       call noah_loop_drv(procbounds, noah_model(nb))
    end do

  end subroutine noah_block_run

  subroutine noah_loop_drv(procbounds, noah_model)

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
    type(noah_type),        intent(inout)    :: noah_model ! land model's variable type

    ! local
    integer                 :: i, de, gridbeg, gridend
    real(kind_phys)         :: foodata(noah_model%static%im)
    !
    real(kind_phys) :: ps        (noah_model%static%im) ! surface pressure (pa)                       im
    real(kind_phys) :: t1        (noah_model%static%im) ! surface layer mean temperature (k)          im
    real(kind_phys) :: q1        (noah_model%static%im) ! surface layer mean specific humidity        im
    integer         :: soiltyp   (noah_model%static%im) ! soil type (integer index)                   im
    integer         :: vegtype   (noah_model%static%im) ! vegetation type (integer index)             im
    real(kind_phys) :: sigmaf    (noah_model%static%im) ! areal fractional cover of green vegetation  im
    real(kind_phys) :: sfcemis   (noah_model%static%im) ! sfc lw emissivity ( fraction )              im
    real(kind_phys) :: dlwflx    (noah_model%static%im) ! total sky sfc downward lw flux ( w/m**2 )   im
    real(kind_phys) :: dswsfc    (noah_model%static%im) ! total sky sfc downward sw flux ( w/m**2 )   im
    real(kind_phys) :: snet      (noah_model%static%im) ! total sky sfc netsw flx into ground(w/m**2) im
    real(kind_phys) :: tg3       (noah_model%static%im) ! deep soil temperature (k)                   im
    real(kind_phys) :: cm        (noah_model%static%im) ! surface exchange coeff for momentum (m/s)   im
    real(kind_phys) :: ch        (noah_model%static%im) ! surface exchange coeff heat & moisture(m/s) im
    real(kind_phys) :: prsl1     (noah_model%static%im) ! sfc layer 1 mean pressure (pa)              im
    real(kind_phys) :: prslki    (noah_model%static%im) !                                             im
    real(kind_phys) :: zf        (noah_model%static%im) ! height of bottom layer (m)                  im
    logical         :: land      (noah_model%static%im) ! = T if a point with any land                im
    real(kind_phys) :: wind      (noah_model%static%im) ! wind speed (m/s)                            im
    integer         :: slopetyp  (noah_model%static%im) ! class of sfc slope (integer index)          im
    real(kind_phys) :: shdmin    (noah_model%static%im) ! min fractional coverage of green veg        im
    real(kind_phys) :: shdmax    (noah_model%static%im) ! max fractnl cover of green veg (not used)   im
    real(kind_phys) :: snoalb    (noah_model%static%im) ! upper bound on max albedo over deep snow    im
    real(kind_phys) :: sfalb     (noah_model%static%im) ! mean sfc diffused sw albedo (fractional)    im
    logical         :: flag_iter (noah_model%static%im) !                                             im
    logical         :: flag_guess(noah_model%static%im) !                                             im
    real(kind_phys) :: bexppert  (noah_model%static%im)
    real(kind_phys) :: xlaipert  (noah_model%static%im)
    real(kind_phys) :: vegfpert  (noah_model%static%im)
    real(kind_phys) :: weasd     (noah_model%static%im) ! water equivalent accumulated snow depth(mm) im
    real(kind_phys) :: snwdph    (noah_model%static%im) ! snow depth (water equiv) over land          im
    real(kind_phys) :: tskin     (noah_model%static%im) ! ground surface skin temperature ( k )       im
    real(kind_phys) :: tprcp     (noah_model%static%im) ! total precipitation                         im
    real(kind_phys) :: srflag    (noah_model%static%im) ! snow/rain flag for precipitation            im
    real(kind_phys) :: canopy    (noah_model%static%im) ! canopy moisture content (m)                 im
    real(kind_phys) :: trans     (noah_model%static%im) ! total plant transpiration (m/s)             im
    real(kind_phys) :: tsurf     (noah_model%static%im) ! surface skin temperature (after iteration)  im
    real(kind_phys) :: z0rl      (noah_model%static%im) ! surface roughness                           im
    real(kind_phys) :: sncovr1   (noah_model%static%im) ! snow cover over land (fractional)            im
    real(kind_phys) :: qsurf     (noah_model%static%im) ! specific humidity at sfc                     im
    real(kind_phys) :: gflux     (noah_model%static%im) ! soil heat flux (w/m**2)                      im
    real(kind_phys) :: drain     (noah_model%static%im) ! subsurface runoff (mm/s)                     im
    real(kind_phys) :: evap      (noah_model%static%im) ! evaperation from latent heat flux            im
    real(kind_phys) :: hflx      (noah_model%static%im) ! sensible heat flux                           im
    real(kind_phys) :: ep        (noah_model%static%im) ! potential evaporation                        im
    real(kind_phys) :: runoff    (noah_model%static%im) ! surface runoff (m/s)                         im
    real(kind_phys) :: cmm       (noah_model%static%im) !                                              im
    real(kind_phys) :: chh       (noah_model%static%im) !                                              im
    real(kind_phys) :: evbs      (noah_model%static%im) ! direct soil evaporation (m/s)                im
    real(kind_phys) :: evcw      (noah_model%static%im) ! canopy water evaporation (m/s)               im
    real(kind_phys) :: sbsno     (noah_model%static%im) ! sublimation/deposit from snopack (m/s)       im
    real(kind_phys) :: snowc     (noah_model%static%im) ! fractional snow cover                        im
    real(kind_phys) :: stm       (noah_model%static%im) ! total soil column moisture content (m)       im
    real(kind_phys) :: snohf     (noah_model%static%im) ! snow/freezing-rain latent heat flux (w/m**2) im
    real(kind_phys) :: smcwlt2   (noah_model%static%im) ! dry soil moisture threshold                  im
    real(kind_phys) :: smcref2   (noah_model%static%im) ! soil moisture threshold                      im
    real(kind_phys) :: wet1      (noah_model%static%im) ! normalized soil wetness                      im

    real(kind_phys) :: ustar     (noah_model%static%im)
    real(kind_phys) :: smc(noah_model%static%im,noah_model%static%km) ! total soil moisture content (fractional)   im,km
    real(kind_phys) :: stc(noah_model%static%im,noah_model%static%km) ! soil temp (k)                              im,km
    real(kind_phys) :: slc(noah_model%static%im,noah_model%static%km) ! liquid soil moisture                       im,km

    ! tmp for testing. Missing imports for these
    real(kind_phys) :: prsik1(noah_model%static%im), z0pert(noah_model%static%im), &
                       ztpert(noah_model%static%im), stress(noah_model%static%im)

    ! tmp for testing. These should be coming from namelist
    real(kind_phys), parameter :: delt = 900.0_kind_phys
    integer, parameter :: ivegsrc        = 1
    integer, parameter :: isot           = 1
    logical, parameter :: lheatstrg      = .false.
    real(kind_phys), parameter :: pertvegf = 0.0_kind_phys
    ! outputs
    real(kind_phys) :: rb_lnd(noah_model%static%im), fm_lnd(noah_model%static%im), &
                       fh_lnd(noah_model%static%im), fm10_lnd(noah_model%static%im), fh2_lnd(noah_model%static%im)
    
    associate(foodata => noah_model%model%foo_atm2lndfield,  &
         ps         => noah_model%model%ps         ,&
         t1         => noah_model%model%t1         ,&
         q1         => noah_model%model%q1         ,&
         dlwflx     => noah_model%model%dlwflx     ,&
         dswsfc     => noah_model%model%dswsfc     ,&
         dswsfci    => noah_model%model%dswsfci    ,&
         wind       => noah_model%model%wind       ,&
         tprcp      => noah_model%model%tprcp      ,&
         im         => noah_model%static%im        ,& 
         km         => noah_model%static%km        ,&
         ! delt       => noah_model%static%delt      ,&
         ! isot       => noah_model%static%isot      ,&
         ! ivegsrc    => noah_model%static%ivegsrc   ,&
         ! pertvegf   => noah_model%static%pertvegf  ,&  
         ! lheatstrg  => noah_model%static%lheatstrg ,& 
         errmsg     => noah_model%static%errmsg    ,& 
         errflg     => noah_model%static%errflg    ,& 
         soiltyp    => noah_model%model%soiltyp    ,&
         vegtype    => noah_model%model%vegtype    ,&
         slopetyp   => noah_model%model%slopetyp   ,&
         sigmaf     => noah_model%model%sigmaf     ,&
         sfcemis    => noah_model%model%sfcemis    ,&
         snet       => noah_model%model%snet       ,&
         tg3        => noah_model%model%tg3        ,&
         cm         => noah_model%model%cm         ,&
         ch         => noah_model%model%ch         ,&
         prsl1      => noah_model%model%prsl1      ,&
         shdmin     => noah_model%model%shdmin     ,&
         shdmax     => noah_model%model%shdmax     ,&
         snoalb     => noah_model%model%snoalb     ,&
         sfalb      => noah_model%model%sfalb      ,&
         zf         => noah_model%model%zf         ,&
         bexppert   => noah_model%model%bexppert   ,&
         xlaipert   => noah_model%model%xlaipert   ,&
         vegfpert   => noah_model%model%vegfpert   ,&
         flag_iter  => noah_model%model%flag_iter  ,&
         flag_guess => noah_model%model%flag_guess ,&
         land       => noah_model%model%land       ,&
         prslki     => noah_model%model%prslki     ,&
         weasd      => noah_model%model%weasd      ,&
         snwdph     => noah_model%model%snwdph     ,&
         tskin      => noah_model%model%tskin      ,&
         srflag     => noah_model%model%srflag     ,&
         canopy     => noah_model%model%canopy     ,&
         trans      => noah_model%model%trans      ,&   
         tsurf      => noah_model%model%tsurf      ,&
         z0rl       => noah_model%model%z0rl       ,&
         ! smc        => noah_model%model%smc        ,&
         ! stc        => noah_model%model%stc        ,&
         ! slc        => noah_model%model%slc        ,&
         sncovr1    => noah_model%model%sncovr1    ,&
         qsurf      => noah_model%model%qsurf      ,&
         gflux      => noah_model%model%gflux      ,&
         drain      => noah_model%model%drain      ,&
         evap       => noah_model%model%evap       ,&
         hflx       => noah_model%model%hflx       ,&
         ep         => noah_model%model%ep         ,&
         runoff     => noah_model%model%runoff     ,&
         cmm        => noah_model%model%cmm        ,&
         chh        => noah_model%model%chh        ,&
         evbs       => noah_model%model%evbs       ,&
         evcw       => noah_model%model%evcw       ,&
         sbsno      => noah_model%model%sbsno      ,&
         snowc      => noah_model%model%snowc      ,&
         stm        => noah_model%model%stm        ,&
         snohf      => noah_model%model%snohf      ,&
         smcwlt2    => noah_model%model%smcwlt2    ,&
         smcref2    => noah_model%model%smcref2    ,&
         wet1       => noah_model%model%wet1       ,&
         !
         rb_lnd     => noah_model%model%rb_lnd     ,&
         fm_lnd     => noah_model%model%fm_lnd     ,&
         fh_lnd     => noah_model%model%fh_lnd     ,&
         fm10_lnd   => noah_model%model%fm10_lnd   ,&
         fh2_lnd    => noah_model%model%fh2_lnd    ,&
         stress     => noah_model%model%stress     ,&
         !
         ustar      => noah_model%model%ustar       &
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
      !write(*,*) 'NLP test: ', noah_model%model%soiltyp ! get all zeros, good
      
      ! de = procbounds%de
      ! !im = procbounds%im
      ! gridbeg = procbounds%gridbeg
      ! gridend = procbounds%gridend

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


      ! foodata(1:gridend-gridbeg+1) = noah_model%model%foo_atm2lndfield(gridbeg:gridend)
      ! foodata = noah_model%model%foo_atm2lndfield(gridbeg:gridend)
      ! foodata = noah_model%model%foo_atm2lndfield


      ! do i = 1,im
      !    write(*,*) 'NLP2: ', de, gridbeg, gridend, foodata(i)
      ! end do

    end associate

  end subroutine noah_loop_drv

end module noah_driver
