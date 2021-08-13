

module noah_type_mod

  use machine, only: kind_phys
  
  implicit none
  save
  private

  !--- parameter constants used for default initializations
  real(kind_phys), parameter :: zero      = 0.0_kind_phys
  real(kind_phys), parameter :: clear_val = zero

type :: noah_control_type
   logical   :: first_time  ! flag for first time step
   integer   :: mype
   integer   :: nblks, blksz, isc, iec, jsc, jec
end type noah_control_type
  
type :: noah_static_type

  integer                 ::   im         ! horiz dimension and num of used pts         1
  integer                 ::   km         ! vertical soil layer dimension               1
  real(kind_phys)    ::   grav       ! constant added to call in ccpp
  real(kind_phys)    ::   cp         ! constant added to call in ccpp
  real(kind_phys)    ::   hvap       ! constant added to call in ccpp
  real(kind_phys)    ::   rd         ! constant added to call in ccpp
  real(kind_phys)    ::   eps        ! constant added to call in ccpp
  real(kind_phys)    ::   epsm1      ! constant added to call in ccpp
  real(kind_phys)    ::   rvrdm1     ! constant added to call in ccpp
  real(kind_phys)    ::   delt       ! time interval (second)                      1
  integer                 ::   isot       ! sfc soil type data source zobler or statsgo
  integer                 ::   ivegsrc    ! sfc veg type data source umd or igbp
  logical                 ::   lheatstrg  ! flag for canopy heat storage parameterization  1
  character(len=128)      ::   errmsg     ! error messaging added to ccpp
  integer                 ::   errflg     ! error messaging added to ccpp
  real(kind_phys)    ::   pertvegf

end type noah_static_type

type sfcprop_type
  real(kind_phys), allocatable      :: landfrac(:)
end type sfcprop_type

type  :: noah_model_type

  real(kind_phys), allocatable :: foo_atm2lndfield(:)
   ! from ufs-land-driver

  real(kind_phys), allocatable :: ps        (:) ! surface pressure (pa)                       im
  real(kind_phys), allocatable :: t1        (:) ! surface layer mean temperature (k)          im
  real(kind_phys), allocatable :: q1        (:) ! surface layer mean specific humidity        im
  integer             , allocatable :: soiltyp   (:) ! soil type (integer index)                   im
  integer             , allocatable :: vegtype   (:) ! vegetation type (integer index)             im
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
  logical             , allocatable :: land      (:) ! = T if a point with any land                im
  real(kind_phys), allocatable :: wind      (:) ! wind speed (m/s)                            im
  integer             , allocatable :: slopetyp  (:) ! class of sfc slope (integer index)          im
  real(kind_phys), allocatable :: shdmin    (:) ! min fractional coverage of green veg        im
  real(kind_phys), allocatable :: shdmax    (:) ! max fractnl cover of green veg (not used)   im
  real(kind_phys), allocatable :: snoalb    (:) ! upper bound on max albedo over deep snow    im
  real(kind_phys), allocatable :: sfalb     (:) ! mean sfc diffused sw albedo (fractional)    im
  logical             , allocatable :: flag_iter (:) !                                             im
  logical             , allocatable :: flag_guess(:) !                                             im
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

  real(kind_phys), allocatable :: smc(:,:) ! total soil moisture content (fractional)   im,km
  real(kind_phys), allocatable :: stc(:,:) ! soil temp (k)                              im,km
  real(kind_phys), allocatable :: slc(:,:) ! liquid soil moisture                       im,km

  ! from sfc_diff
  real(kind_phys), allocatable :: rb_lnd   (:)
  real(kind_phys), allocatable :: fm_lnd   (:)
  real(kind_phys), allocatable :: fh_lnd   (:)
  real(kind_phys), allocatable :: fm10_lnd (:)
  real(kind_phys), allocatable :: fh2_lnd  (:)
  real(kind_phys), allocatable :: stress   (:)  
  real(kind_phys), allocatable :: ustar    (:)  
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

    class(noah_type)            :: nh
    integer,                intent(in) :: im

    allocate(nh%model%foo_atm2lndfield        (im))
    allocate(nh%model%ps        (im))
    allocate(nh%model%t1        (im))
    allocate(nh%model%q1        (im))
    allocate(nh%model%soiltyp   (im))
    allocate(nh%model%vegtype   (im))
    allocate(nh%model%sigmaf    (im))
    allocate(nh%model%sfcemis   (im))
    allocate(nh%model%dlwflx    (im))
    allocate(nh%model%dswsfc    (im))
    allocate(nh%model%dswsfci   (im))
    allocate(nh%model%snet      (im))
    allocate(nh%model%tg3       (im))
    allocate(nh%model%cm        (im))
    allocate(nh%model%ch        (im))
    allocate(nh%model%prsl1     (im))
    allocate(nh%model%prslki    (im))
    allocate(nh%model%zf        (im))
    allocate(nh%model%land      (im))
    allocate(nh%model%wind      (im))
    allocate(nh%model%slopetyp  (im))
    allocate(nh%model%shdmin    (im))
    allocate(nh%model%shdmax    (im))
    allocate(nh%model%snoalb    (im))
    allocate(nh%model%sfalb     (im))
    allocate(nh%model%flag_iter (im))
    allocate(nh%model%flag_guess(im))
    allocate(nh%model%bexppert  (im))
    allocate(nh%model%xlaipert  (im))
    allocate(nh%model%vegfpert  (im))
    allocate(nh%model%weasd     (im))
    allocate(nh%model%snwdph    (im))
    allocate(nh%model%tskin     (im))
    allocate(nh%model%tprcp     (im))
    allocate(nh%model%srflag    (im))
    allocate(nh%model%canopy    (im))
    allocate(nh%model%trans     (im))
    allocate(nh%model%tsurf     (im))
    allocate(nh%model%z0rl      (im))
    allocate(nh%model%sncovr1   (im))
    allocate(nh%model%qsurf     (im))
    allocate(nh%model%gflux     (im))
    allocate(nh%model%drain     (im))
    allocate(nh%model%evap      (im))
    allocate(nh%model%hflx      (im))
    allocate(nh%model%ep        (im))
    allocate(nh%model%runoff    (im))
    allocate(nh%model%cmm       (im))
    allocate(nh%model%chh       (im))
    allocate(nh%model%evbs      (im))
    allocate(nh%model%evcw      (im))
    allocate(nh%model%sbsno     (im))
    allocate(nh%model%snowc     (im))
    allocate(nh%model%stm       (im))
    allocate(nh%model%snohf     (im))
    allocate(nh%model%smcwlt2   (im))
    allocate(nh%model%smcref2   (im))
    allocate(nh%model%wet1      (im))

    !allocate(nh%model%smc  (im,num_soil_levels))
    !allocate(nh%model%stc  (im,num_soil_levels))
    !allocate(nh%model%slc  (im,num_soil_levels))
    ! sfc_diff
    allocate(nh%model%rb_lnd  (im))
    allocate(nh%model%fm_lnd  (im))
    allocate(nh%model%fh_lnd  (im))
    allocate(nh%model%fm10_lnd(im))
    allocate(nh%model%fh2_lnd (im))
    allocate(nh%model%stress  (im))
    allocate(nh%model%ustar  (im))


    allocate(nh%sfcprop%landfrac (im))

    
    ! --------------------------------------------------------
    nh%control%first_time = .true.
    nh%control%mype       = -999

    !nh%static%im         = -999
    nh%static%km         = 4 ! tmp for testing. This should come from nml
    nh%static%grav       = zero !huge
    nh%static%cp         = zero !huge
    nh%static%hvap       = zero !huge
    nh%static%rd         = zero !huge
    nh%static%eps        = zero !huge
    nh%static%epsm1      = zero !huge
    nh%static%rvrdm1     = zero !huge
    nh%static%delt       = zero !huge
    nh%static%isot       = zero !huge
    nh%static%ivegsrc    = zero !huge
    nh%static%lheatstrg  = .false.
    nh%static%errmsg     = ""
    nh%static%errflg     = zero !huge
    nh%static%pertvegf   = zero !huge

    nh%model%foo_atm2lndfield         = zero !huge
    nh%model%ps         = zero !huge
    nh%model%t1         = zero !huge
    nh%model%q1         = zero !huge
    nh%model%soiltyp    = zero !huge
    nh%model%vegtype    = zero !huge
    nh%model%sigmaf     = zero !huge
    nh%model%sfcemis    = zero !huge
    nh%model%dlwflx     = zero !huge
    nh%model%dswsfc     = zero !huge
    nh%model%dswsfci    = zero !huge
    nh%model%snet       = zero !huge
    nh%model%tg3        = zero !huge
    nh%model%cm         = zero !huge
    nh%model%ch         = zero !huge
    nh%model%prsl1      = zero !huge
    nh%model%prslki     = zero !huge
    nh%model%zf         = zero !huge
    nh%model%land       = .false.
    nh%model%wind       = zero !huge
    nh%model%slopetyp   = zero !huge
    nh%model%shdmin     = zero !huge
    nh%model%shdmax     = zero !huge
    nh%model%snoalb     = zero !huge
    nh%model%sfalb      = zero !huge
    nh%model%flag_iter  = .false.
    nh%model%flag_guess = .false.
    nh%model%bexppert   = zero !huge
    nh%model%xlaipert   = zero !huge
    nh%model%vegfpert   = zero !huge
    nh%model%weasd      = zero !huge
    nh%model%snwdph     = zero !huge
    nh%model%tskin      = zero !huge
    nh%model%tprcp      = zero !huge
    nh%model%srflag     = zero !huge
    nh%model%canopy     = zero !huge
    nh%model%trans      = zero !huge
    nh%model%tsurf      = zero !huge
    nh%model%z0rl       = zero !huge
    nh%model%sncovr1    = zero !huge
    nh%model%qsurf      = zero !huge
    nh%model%gflux      = zero !huge
    nh%model%drain      = zero !huge
    nh%model%evap       = zero !huge
    nh%model%hflx       = zero !huge
    nh%model%ep         = zero !huge
    nh%model%runoff     = zero !huge
    nh%model%cmm        = zero !huge
    nh%model%chh        = zero !huge
    nh%model%evbs       = zero !huge
    nh%model%evcw       = zero !huge
    nh%model%sbsno      = zero !huge
    nh%model%snowc      = zero !huge
    nh%model%stm        = zero !huge
    nh%model%snohf      = zero !huge
    nh%model%smcwlt2    = zero !huge
    nh%model%smcref2    = zero !huge
    nh%model%wet1       = zero !huge
    nh%model%smc        = zero !huge
    nh%model%stc        = zero !huge
    nh%model%slc        = zero !huge

    nh%model%rb_lnd     = zero
    nh%model%fm_lnd     = zero
    nh%model%fh_lnd     = zero
    nh%model%fm10_lnd   = zero
    nh%model%fh2_lnd    = zero
    nh%model%stress     = zero
    nh%model%ustar      = zero

    nh%sfcprop%landfrac  = zero

  end subroutine Create




end module noah_type_mod
