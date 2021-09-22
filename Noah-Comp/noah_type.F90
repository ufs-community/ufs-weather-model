

module noah_type_mod

  use machine, only: kind_phys

  implicit none
  save
  private

  !--- parameter constants used for default initializations
  real(kind_phys), parameter :: zero      = 0.0_kind_phys
  real(kind_phys), parameter :: clear_val = -9999_kind_phys

  type :: noah_control_type
     logical   :: first_time  ! flag for first time step
     integer   :: mype
     integer   :: nblks, blksz, isc, iec, jsc, jec
  end type noah_control_type

  type :: noah_static_type

     integer            ::   im         ! horiz dimension and num of used pts         1
     integer            ::   km         ! vertical soil layer dimension               1
     real(kind_phys)    ::   grav       ! constant added to call in ccpp
     real(kind_phys)    ::   cp         ! constant added to call in ccpp
     real(kind_phys)    ::   hvap       ! constant added to call in ccpp
     real(kind_phys)    ::   rd         ! constant added to call in ccpp
     real(kind_phys)    ::   eps        ! constant added to call in ccpp
     real(kind_phys)    ::   epsm1      ! constant added to call in ccpp
     real(kind_phys)    ::   rvrdm1     ! constant added to call in ccpp
     real(kind_phys)    ::   delt       ! time interval (second)                      1
     integer            ::   isot       ! sfc soil type data source zobler or statsgo
     integer            ::   ivegsrc    ! sfc veg type data source umd or igbp
     logical            ::   lheatstrg  ! flag for canopy heat storage parameterization  1
     character(len=128) ::   errmsg     ! error messaging added to ccpp
     integer            ::   errflg     ! error messaging added to ccpp
     real(kind_phys)    ::   pertvegf
     logical            ::   thsfc_loc ! this should be changed to match same FV3 var

  end type noah_static_type
  
  type sfcprop_type
     real(kind_phys), allocatable  :: landfrac(:)
     real(kind_phys), allocatable  :: slmsk(:)
     real(kind_phys), allocatable  :: tsfcl(:)  
     real(kind_phys), allocatable  :: weasd(:)  ! aka sheleg in sfc file
     real(kind_phys), allocatable  :: tg3(:)
     real(kind_phys), allocatable  :: zorll(:)  ! note, z0rl over land
     real(kind_phys), allocatable  :: alvsf(:)
     real(kind_phys), allocatable  :: alvwf(:)
     real(kind_phys), allocatable  :: alnsf(:)
     real(kind_phys), allocatable  :: alnwf(:)
     real(kind_phys), allocatable  :: facsf(:)
     real(kind_phys), allocatable  :: facwf(:)
     real(kind_phys), allocatable  :: vfrac(:)
     real(kind_phys), allocatable  :: canopy(:)
     real(kind_phys), allocatable  :: f10m(:)
     real(kind_phys), allocatable  :: t2m(:)
     real(kind_phys), allocatable  :: q2m(:)
     real(kind_phys), allocatable  :: vtype(:)
     real(kind_phys), allocatable  :: stype(:)
     real(kind_phys), allocatable  :: uustar(:)
     real(kind_phys), allocatable  :: ffmm(:)
     real(kind_phys), allocatable  :: ffhh(:)
     real(kind_phys), allocatable  :: hice(:)
     real(kind_phys), allocatable  :: fice(:)
     real(kind_phys), allocatable  :: tisfc(:)
     real(kind_phys), allocatable  :: tprcp(:)
     real(kind_phys), allocatable  :: srflag(:)
     real(kind_phys), allocatable  :: snowd(:)  ! aka snwdph in sfc file
     real(kind_phys), allocatable  :: shdmin(:)
     real(kind_phys), allocatable  :: shdmax(:)
     real(kind_phys), allocatable  :: slope(:)
     real(kind_phys), allocatable  :: snoalb(:)
     real(kind_phys), allocatable  :: sncovr(:)

     ! JP TODO: allocate these properly
     real(kind_phys), allocatable  :: stc(:,:)
     real(kind_phys), allocatable  :: smc(:,:)
     real(kind_phys), allocatable  :: slc(:,:)     
  end type sfcprop_type

  type  :: noah_model_type

     real(kind_phys), allocatable :: foo_atm2lndfield(:)
     ! from ufs-land-driver

     real(kind_phys), allocatable :: ps        (:) ! surface pressure (pa)                       im
     real(kind_phys), allocatable :: t1        (:) ! surface layer mean temperature (k)          im
     real(kind_phys), allocatable :: q1        (:) ! surface layer mean specific humidity        im
     integer        , allocatable :: soiltyp   (:) ! soil type (integer index)                   im
     integer        , allocatable :: vegtype   (:) ! vegetation type (integer index)             im
     real(kind_phys), allocatable :: sigmaf    (:) ! areal fractional cover of green vegetation  im
     real(kind_phys), allocatable :: sfcemis   (:) ! sfc lw emissivity ( fraction )              im
     real(kind_phys), allocatable :: dlwflx    (:) ! total sky sfc downward lw flux ( w/m**2 )   im
     real(kind_phys), allocatable :: dswsfc    (:) ! total sky sfc downward sw flux ( w/m**2 )   im
     real(kind_phys), allocatable :: dswsfci   (:) ! inst  sky sfc downward sw flux ( w/m**2 )   im
     real(kind_phys), allocatable :: snet      (:) ! total sky sfc netsw flx into ground(w/m**2) im
     real(kind_phys), allocatable :: tg3       (:) ! deep soil temperature (k)                   im
     real(kind_phys), allocatable :: cm        (:) ! surface exchange coeff for momentum (m/s)   im
     real(kind_phys), allocatable :: ch        (:) ! surface exchange coeff heat & moisture(m/s) im
     real(kind_phys), allocatable :: prsl1     (:) ! sfc layer 1 mean pressure (pa)              im
     real(kind_phys), allocatable :: prslki    (:) !                                             im
     real(kind_phys), allocatable :: zf        (:) ! height of bottom layer (m)                  im
     logical        , allocatable :: land      (:) ! = T if a point with any land                im
     real(kind_phys), allocatable :: wind      (:) ! wind speed (m/s)                            im
     integer        , allocatable :: slopetyp  (:) ! class of sfc slope (integer index)          im
     real(kind_phys), allocatable :: shdmin    (:) ! min fractional coverage of green veg        im
     real(kind_phys), allocatable :: shdmax    (:) ! max fractnl cover of green veg (not used)   im
     real(kind_phys), allocatable :: snoalb    (:) ! upper bound on max albedo over deep snow    im
     real(kind_phys), allocatable :: sfalb     (:) ! mean sfc diffused sw albedo (fractional)    im
     logical        , allocatable :: flag_iter (:) !                                             im
     logical        , allocatable :: flag_guess(:) !                                             im
     real(kind_phys), allocatable :: bexppert  (:)
     real(kind_phys), allocatable :: xlaipert  (:)
     real(kind_phys), allocatable :: vegfpert  (:)
     real(kind_phys), allocatable :: weasd     (:) ! water equivalent accumulated snow depth(mm) im
     real(kind_phys), allocatable :: snwdph    (:) ! snow depth (water equiv) over land          im
     real(kind_phys), allocatable :: tskin     (:) ! ground surface skin temperature ( k )       im
     real(kind_phys), allocatable :: tprcp     (:) ! total precipitation                         im
     real(kind_phys), allocatable :: srflag    (:) ! snow/rain flag for precipitation            im
     real(kind_phys), allocatable :: canopy    (:) ! canopy moisture content (m)                 im
     real(kind_phys), allocatable :: trans     (:) ! total plant transpiration (m/s)             im
     real(kind_phys), allocatable :: tsurf     (:) ! surface skin temperature (after iteration)  im
     real(kind_phys), allocatable :: z0rl      (:) ! surface roughness                           im
     real(kind_phys), allocatable :: sncovr1   (:) ! snow cover over land (fractional)            im
     real(kind_phys), allocatable :: qsurf     (:) ! specific humidity at sfc                     im
     real(kind_phys), allocatable :: gflux     (:) ! soil heat flux (w/m**2)                      im
     real(kind_phys), allocatable :: drain     (:) ! subsurface runoff (mm/s)                     im
     real(kind_phys), allocatable :: evap      (:) ! evaperation from latent heat flux            im
     real(kind_phys), allocatable :: hflx      (:) ! sensible heat flux                           im
     real(kind_phys), allocatable :: ep        (:) ! potential evaporation                        im
     real(kind_phys), allocatable :: runoff    (:) ! surface runoff (m/s)                         im
     real(kind_phys), allocatable :: cmm       (:) !                                              im
     real(kind_phys), allocatable :: chh       (:) !                                              im
     real(kind_phys), allocatable :: evbs      (:) ! direct soil evaporation (m/s)                im
     real(kind_phys), allocatable :: evcw      (:) ! canopy water evaporation (m/s)               im
     real(kind_phys), allocatable :: sbsno     (:) ! sublimation/deposit from snopack (m/s)       im
     real(kind_phys), allocatable :: snowc     (:) ! fractional snow cover                        im
     real(kind_phys), allocatable :: stm       (:) ! total soil column moisture content (m)       im
     real(kind_phys), allocatable :: snohf     (:) ! snow/freezing-rain latent heat flux (w/m**2) im
     real(kind_phys), allocatable :: smcwlt2   (:) ! dry soil moisture threshold                  im
     real(kind_phys), allocatable :: smcref2   (:) ! soil moisture threshold                      im
     real(kind_phys), allocatable :: wet1      (:) ! normalized soil wetness                      im
     real(kind_phys), allocatable :: prslk1    (:)

     ! JP TODO: allocate these properly
     real(kind_phys), allocatable :: smc(:,:) ! total soil moisture content (fractional)   im,km
     real(kind_phys), allocatable :: stc(:,:) ! soil temp (k)                              im,km
     real(kind_phys), allocatable :: slc(:,:) ! liquid soil moisture                       im,km

     ! rad
     real(kind_phys), allocatable :: albdvis_lnd (:)
     real(kind_phys), allocatable :: albdnir_lnd (:)
     real(kind_phys), allocatable :: albivis_lnd (:)
     real(kind_phys), allocatable :: albinir_lnd (:)
     real(kind_phys), allocatable :: adjvisbmd   (:)
     real(kind_phys), allocatable :: adjnirbmd   (:)
     real(kind_phys), allocatable :: adjvisdfd   (:)
     real(kind_phys), allocatable :: adjnirdfd   (:)

     ! from sfc_diff
     real(kind_phys), allocatable :: rb_lnd   (:)
     real(kind_phys), allocatable :: fm_lnd   (:)
     real(kind_phys), allocatable :: fh_lnd   (:)
     real(kind_phys), allocatable :: fm10_lnd (:)
     real(kind_phys), allocatable :: fh2_lnd  (:)
     real(kind_phys), allocatable :: stress   (:)  
     real(kind_phys), allocatable :: ustar    (:)
     real(kind_phys), allocatable :: garea    (:)

  end type noah_model_type


  type, public :: noah_type
     type(noah_static_type)  :: static
     type(noah_model_type)   :: model
     type(noah_control_type) :: control
     type(sfcprop_type)      :: sfcprop
   contains

     procedure, public  :: Create

  end type noah_type

contains

  subroutine Create(nh, im)

    implicit none

    class(noah_type)    :: nh
    integer, intent(in) :: im

    integer,parameter :: km = 4 ! tmp for testing. This should come from nml
    
    ! --------------------------------------------
    nh%control%first_time = .true.
    nh%control%mype       = clear_val

    ! --------------------------------------------
    !nh%static%im         = clear_val
    nh%static%km         = km 
    nh%static%grav       = clear_val
    nh%static%cp         = clear_val
    nh%static%hvap       = clear_val
    nh%static%rd         = clear_val
    nh%static%eps        = clear_val
    nh%static%epsm1      = clear_val
    nh%static%rvrdm1     = clear_val
    nh%static%delt       = clear_val
    nh%static%isot       = 1  ! TMP for testing. TODO: read in
    nh%static%ivegsrc    = 1  ! TMP for testing. TODO: read in
    nh%static%lheatstrg  = .false.
    nh%static%errmsg     = ""
    nh%static%errflg     = clear_val
    nh%static%pertvegf   = clear_val
    nh%static%thsfc_loc  = .true.
    ! --------------------------------------------
    ! --------------------------------------------    
    allocate(nh%model%foo_atm2lndfield        (im))
    allocate(nh%model%ps            (im))
    allocate(nh%model%t1            (im))
    allocate(nh%model%q1            (im))
    allocate(nh%model%soiltyp       (im))
    allocate(nh%model%vegtype       (im))
    allocate(nh%model%sigmaf        (im))
    allocate(nh%model%sfcemis       (im))
    allocate(nh%model%dlwflx        (im))
    allocate(nh%model%dswsfc        (im))
    allocate(nh%model%dswsfci       (im))
    allocate(nh%model%snet          (im))
    allocate(nh%model%tg3           (im))
    allocate(nh%model%cm            (im))
    allocate(nh%model%ch            (im))
    allocate(nh%model%prsl1         (im))
    allocate(nh%model%prslki        (im))
    allocate(nh%model%zf            (im))
    allocate(nh%model%land          (im))
    allocate(nh%model%wind          (im))
    allocate(nh%model%slopetyp      (im))
    allocate(nh%model%shdmin        (im))
    allocate(nh%model%shdmax        (im))
    allocate(nh%model%snoalb        (im))
    allocate(nh%model%sfalb         (im))
    allocate(nh%model%flag_iter     (im))
    allocate(nh%model%flag_guess    (im))
    allocate(nh%model%bexppert      (im))
    allocate(nh%model%xlaipert      (im))
    allocate(nh%model%vegfpert      (im))
    allocate(nh%model%weasd         (im))
    allocate(nh%model%snwdph        (im))
    allocate(nh%model%tskin         (im))
    allocate(nh%model%tprcp         (im))
    allocate(nh%model%srflag        (im))
    allocate(nh%model%canopy        (im))
    allocate(nh%model%trans         (im))
    allocate(nh%model%tsurf         (im))
    allocate(nh%model%z0rl          (im))
    allocate(nh%model%sncovr1       (im))
    allocate(nh%model%qsurf         (im))
    allocate(nh%model%gflux         (im))
    allocate(nh%model%drain         (im))
    allocate(nh%model%evap          (im))
    allocate(nh%model%hflx          (im))
    allocate(nh%model%ep            (im))
    allocate(nh%model%runoff        (im))
    allocate(nh%model%cmm           (im))
    allocate(nh%model%chh           (im))
    allocate(nh%model%evbs          (im))
    allocate(nh%model%evcw          (im))
    allocate(nh%model%sbsno         (im))
    allocate(nh%model%snowc         (im))
    allocate(nh%model%stm           (im))
    allocate(nh%model%snohf         (im))
    allocate(nh%model%smcwlt2       (im))
    allocate(nh%model%smcref2       (im))
    allocate(nh%model%wet1          (im))
    allocate(nh%model%albdvis_lnd   (im))
    allocate(nh%model%albdnir_lnd   (im))
    allocate(nh%model%albivis_lnd   (im))
    allocate(nh%model%albinir_lnd   (im))
    allocate(nh%model%adjvisbmd     (im))
    allocate(nh%model%adjnirbmd     (im))
    allocate(nh%model%adjvisdfd     (im))
    allocate(nh%model%adjnirdfd     (im))
    allocate(nh%model%prslk1        (im))
    allocate(nh%model%smc       (im,km))
    allocate(nh%model%stc       (im,km))
    allocate(nh%model%slc       (im,km))
    ! sfc_diff
    allocate(nh%model%rb_lnd        (im))
    allocate(nh%model%fm_lnd        (im))
    allocate(nh%model%fh_lnd        (im))
    allocate(nh%model%fm10_lnd      (im))
    allocate(nh%model%fh2_lnd       (im))
    allocate(nh%model%stress        (im))
    allocate(nh%model%ustar         (im))
    allocate(nh%model%garea         (im))
    
    !! Sfcprop -------------------------
    allocate(nh%sfcprop%landfrac (im))
    allocate(nh%sfcprop%slmsk    (im))      
    allocate(nh%sfcprop%tsfcl    (im))
    allocate(nh%sfcprop%weasd    (im))
    allocate(nh%sfcprop%tg3      (im))    
    allocate(nh%sfcprop%zorll    (im))     
    allocate(nh%sfcprop%alvsf    (im))      
    allocate(nh%sfcprop%alvwf    (im))      
    allocate(nh%sfcprop%alnsf    (im))      
    allocate(nh%sfcprop%alnwf    (im))      
    allocate(nh%sfcprop%facsf    (im))      
    allocate(nh%sfcprop%facwf    (im))      
    allocate(nh%sfcprop%vfrac    (im))      
    allocate(nh%sfcprop%canopy   (im))       
    allocate(nh%sfcprop%f10m     (im))     
    allocate(nh%sfcprop%t2m      (im))    
    allocate(nh%sfcprop%q2m      (im))    
    allocate(nh%sfcprop%vtype    (im))      
    allocate(nh%sfcprop%stype    (im))      
    allocate(nh%sfcprop%uustar   (im))       
    allocate(nh%sfcprop%ffmm     (im))     
    allocate(nh%sfcprop%ffhh     (im))     
    allocate(nh%sfcprop%hice     (im))     
    allocate(nh%sfcprop%fice     (im))     
    allocate(nh%sfcprop%tisfc    (im))      
    allocate(nh%sfcprop%tprcp    (im))      
    allocate(nh%sfcprop%srflag   (im))       
    allocate(nh%sfcprop%snowd    (im))
    allocate(nh%sfcprop%shdmin   (im))       
    allocate(nh%sfcprop%shdmax   (im))       
    allocate(nh%sfcprop%slope    (im))      
    allocate(nh%sfcprop%snoalb   (im))       
    allocate(nh%sfcprop%sncovr   (im))       
    
    allocate(nh%sfcprop%smc      (im,km))
    allocate(nh%sfcprop%stc      (im,km))
    allocate(nh%sfcprop%slc      (im,km))

    ! --------------------------------------------------------

    nh%model%foo_atm2lndfield  = clear_val
    nh%model%ps         = clear_val
    nh%model%t1         = clear_val
    nh%model%q1         = clear_val
    nh%model%soiltyp    = clear_val
    nh%model%vegtype    = clear_val
    nh%model%sigmaf     = clear_val
    nh%model%sfcemis    = clear_val
    nh%model%dlwflx     = clear_val
    nh%model%dswsfc     = clear_val
    nh%model%dswsfci    = clear_val
    nh%model%snet       = clear_val
    nh%model%tg3        = clear_val
    nh%model%cm         = clear_val
    nh%model%ch         = clear_val
    nh%model%prsl1      = clear_val
    nh%model%prslki     = clear_val
    nh%model%zf         = clear_val
    nh%model%land       = .false.
    nh%model%wind       = clear_val
    nh%model%slopetyp   = clear_val
    nh%model%shdmin     = clear_val
    nh%model%shdmax     = clear_val
    nh%model%snoalb     = clear_val
    nh%model%sfalb      = clear_val
    nh%model%flag_iter  = .false.
    nh%model%flag_guess = .false.
    nh%model%bexppert   = clear_val
    nh%model%xlaipert   = clear_val
    nh%model%vegfpert   = clear_val
    nh%model%weasd      = clear_val
    nh%model%snwdph     = clear_val
    nh%model%tskin      = clear_val
    nh%model%tprcp      = clear_val
    nh%model%srflag     = clear_val
    nh%model%canopy     = clear_val
    nh%model%trans      = clear_val
    nh%model%tsurf      = clear_val
    nh%model%z0rl       = clear_val
    nh%model%sncovr1    = clear_val
    nh%model%qsurf      = clear_val
    nh%model%gflux      = clear_val
    nh%model%drain      = clear_val
    nh%model%evap       = clear_val
    nh%model%hflx       = clear_val
    nh%model%ep         = clear_val
    nh%model%runoff     = clear_val
    nh%model%cmm        = clear_val
    nh%model%chh        = clear_val
    nh%model%evbs       = clear_val
    nh%model%evcw       = clear_val
    nh%model%sbsno      = clear_val
    nh%model%snowc      = clear_val
    nh%model%stm        = clear_val
    nh%model%snohf      = clear_val
    nh%model%smcwlt2    = clear_val
    nh%model%smcref2    = clear_val
    nh%model%wet1       = clear_val
    nh%model%smc        = clear_val
    nh%model%stc        = clear_val
    nh%model%slc        = clear_val     

    nh%model%rb_lnd     = clear_val
    nh%model%fm_lnd     = clear_val
    nh%model%fh_lnd     = clear_val
    nh%model%fm10_lnd   = clear_val
    nh%model%fh2_lnd    = clear_val
    nh%model%stress     = clear_val
    nh%model%ustar      = clear_val

    !! Surf Prop
    nh%sfcprop%landfrac    = clear_val
    nh%sfcprop%landfrac    = clear_val
    nh%sfcprop%slmsk       = clear_val
    nh%sfcprop%tsfcl       = clear_val
    nh%sfcprop%weasd       = clear_val
    nh%sfcprop%tg3         = clear_val
    nh%sfcprop%zorll       = clear_val
    nh%sfcprop%alvsf       = clear_val
    nh%sfcprop%alvwf       = clear_val
    nh%sfcprop%alnsf       = clear_val
    nh%sfcprop%alnwf       = clear_val
    nh%sfcprop%facsf       = clear_val
    nh%sfcprop%facwf       = clear_val
    nh%sfcprop%vfrac       = clear_val
    nh%sfcprop%canopy      = clear_val
    nh%sfcprop%f10m        = clear_val
    nh%sfcprop%t2m         = clear_val
    nh%sfcprop%q2m         = clear_val
    nh%sfcprop%vtype       = clear_val
    nh%sfcprop%stype       = clear_val
    nh%sfcprop%uustar      = clear_val
    nh%sfcprop%ffmm        = clear_val
    nh%sfcprop%ffhh        = clear_val
    nh%sfcprop%hice        = clear_val
    nh%sfcprop%fice        = clear_val
    nh%sfcprop%tisfc       = clear_val
    nh%sfcprop%tprcp       = clear_val
    nh%sfcprop%srflag      = clear_val
    nh%sfcprop%snowd       = clear_val
    nh%sfcprop%shdmin      = clear_val
    nh%sfcprop%shdmax      = clear_val
    nh%sfcprop%slope       = clear_val
    nh%sfcprop%snoalb      = clear_val
    nh%sfcprop%sncovr      = clear_val

    nh%sfcprop%smc         = clear_val
    nh%sfcprop%stc         = clear_val
    nh%sfcprop%slc         = clear_val     

    nh%model%garea       = zero
    
  end subroutine Create




end module noah_type_mod
