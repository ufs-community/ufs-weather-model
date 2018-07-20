      MODULE module_microphysics
!
      USE MACHINE , ONLY : kind_phys
      USE FUNCPHYS
      USE PHYSCONS, CP => con_CP, RD => con_RD, RV => con_RV            &
     &,             T0C => con_T0C, HVAP => con_HVAP, HFUS => con_HFUS  &
     &,             EPS => con_EPS, EPSM1 => con_EPSM1                  &
     &,             EPS1 => con_FVirt, pi => con_pi, grav => con_g 
      implicit none
!
!--- Common block of constants used in column microphysics
!
      real,private ::  ABFR, CBFR, CIACW, CIACR, C_N0r0,                &
     &CN0r0, CN0r_DMRmin, CN0r_DMRmax, CRACW, CRAUT, ESW0,              &
     &QAUTx, RFmax, RQR_DR1, RQR_DR2, RQR_DR3, RQR_DRmin,               &
     &RQR_DRmax, RR_DRmin, RR_DR1, RR_DR2, RR_DR3, RR_DRmax
!
      integer, private :: mic_step
!
!--- Common block for lookup table used in calculating growth rates of
!    nucleated ice crystals growing in water saturated conditions
!--- Discretized growth rates of small ice crystals after their nucleation
!     at 1 C intervals from -1 C to -35 C, based on calculations by Miller
!     and Young (1979, JAS) after 600 s of growth.  Resultant growth rates
!     are multiplied by physics time step in GSMCONST.
!
      INTEGER, PRIVATE,PARAMETER :: MY_T1=1, MY_T2=35
      REAL,PRIVATE,DIMENSION(MY_T1:MY_T2) :: MY_GROWTH
!
!--- Parameters for ice lookup tables, which establish the range of mean ice
!    particle diameters; from a minimum mean diameter of 0.05 mm (DMImin) to a
!    maximum mean diameter of 1.00 mm (DMImax).  The tables store solutions
!    at 1 micron intervals (DelDMI) of mean ice particle diameter.
!
      REAL, PRIVATE,PARAMETER :: DMImin=.05e-3,      DMImax=1.e-3,      &
     &                           XMImin=1.e6*DMImin, XMImax=1.e6*DMImax,&
     &                           DelDMI=1.e-6
      INTEGER, PRIVATE,PARAMETER :: MDImin=XMImin, MDImax=XMImax
!
!--- Various ice lookup tables
!
      REAL, PRIVATE,DIMENSION(MDImin:MDImax) ::                         &
     &      ACCRI,MASSI,SDENS,VSNOWI,VENTI1,VENTI2
!
!--- Mean rain drop diameters varying from 50 microns (0.05 mm) to 450 microns
!      (0.45 mm), assuming an exponential size distribution.
!
      REAL, PRIVATE,PARAMETER :: DMRmin=.05e-3,      DMRmax=.45e-3,     &
     &                           XMRmin=1.e6*DMRmin, XMRmax=1.e6*DMRmax,&
     &                           DelDMR=1.e-6,       NLImin=100.
!    &,                          NLImin=100., NLImax=20.E3
      INTEGER, PRIVATE,PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax
!
!--- Factor of 1.5 for RECImin, RESNOWmin, & RERAINmin accounts for
!    integrating exponential distributions for effective radius
!    (i.e., the r**3/r**2 moments).
!
!     INTEGER, PRIVATE, PARAMETER :: INDEXSmin=300
!!    INTEGER, PRIVATE, PARAMETER :: INDEXSmin=200
      INTEGER, PRIVATE, PARAMETER :: INDEXSmin=100
      REAL, PRIVATE, PARAMETER :: RERAINmin=1.5*XMRmin                  &
!    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=8.0
!    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=7.5
     &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=10.
!    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=15.
!    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=5.

!
!--- Various rain lookup tables
!--- Rain lookup tables for mean rain drop diameters from DMRmin to DMRmax,
!      assuming exponential size distributions for the rain drops
!
      REAL, PRIVATE,DIMENSION(MDRmin:MDRmax)::                          &
     &      ACCRR,MASSR,RRATE,VRAIN,VENTR1,VENTR2
!
!--- Common block for riming tables
!--- VEL_RF - velocity increase of rimed particles as functions of crude
!      particle size categories (at 0.1 mm intervals of mean ice particle
!      sizes) and rime factor (different values of Rime Factor of 1.1**N,
!      where N=0 to Nrime).
!
      INTEGER, PRIVATE,PARAMETER :: Nrime=40
      REAL, DIMENSION(2:9,0:Nrime),PRIVATE :: VEL_RF
!
!--- The following variables are for microphysical statistics
!
      INTEGER, PARAMETER :: ITLO=-60, ITHI=40
      INTEGER  NSTATS(ITLO:ITHI,4)
      REAL     QMAX(ITLO:ITHI,5),  QTOT(ITLO:ITHI,22)
!
      REAL, PRIVATE,  PARAMETER ::                                      &
!    &  T_ICE=-10., T_ICE_init=-5.      !- Ver1
!!!  &, T_ICE=-20.                      !- Ver2
     &  T_ICE=-40., T_ICE_init=-15.     !- Ver2
!    &  T_ICE=-30., T_ICE_init=-5.      !- Ver2
!
!     Some other miscellaneous parameters
!
      REAL, PRIVATE, PARAMETER :: Thom=T_ICE, TNW=50., TOLER=1.0E-20    &
!     REAL, PRIVATE, PARAMETER :: Thom=T_ICE, TNW=50., TOLER=5.E-7
!     REAL, PRIVATE, PARAMETER :: Thom=-35., TNW=50., TOLER=5.E-7

!    &, emisCU=.75/1.66       ! Used for convective cloud l/w emissivity

! Assume fixed cloud ice effective radius
     &, RECICE=RECImin                                                  &
     &, EPSQ=1.0E-20                                                    &
!    &, EPSQ=1.E-12                                                     &
     &, FLG0P1=0.1, FLG0P2=0.2, FLG1P0=1.0
!
!
      CONTAINS
!
!#######################################################################
!------- Initialize constants & lookup tables for microphysics ---------
!#######################################################################
!
      SUBROUTINE GSMCONST (DTPG,mype,first)
!
      implicit none
!-------------------------------------------------------------------------------
!---  SUBPROGRAM DOCUMENTATION BLOCK
!   PRGRMMR: Ferrier         ORG: W/NP22     DATE: February 2001
!-------------------------------------------------------------------------------
! ABSTRACT:
!   * Reads various microphysical lookup tables used in COLUMN_MICRO
!   * Lookup tables were created "offline" and are read in during execution
!   * Creates lookup tables for saturation vapor pressure w/r/t water & ice
!-------------------------------------------------------------------------------
!
! USAGE: CALL GSMCONST FROM SUBROUTINE GSMDRIVE AT MODEL START TIME
!
!   INPUT ARGUMENT LIST:
!       DTPH - physics time step (s)
!
!   OUTPUT ARGUMENT LIST:
!     NONE
!
!   OUTPUT FILES:
!     NONE
!
!
!   SUBROUTINES:
!     MY_GROWTH_RATES - lookup table for growth of nucleated ice
!
!   UNIQUE: NONE
!
!   LIBRARY: NONE
!
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!
      integer mype
      real    dtpg
      logical first
!
!--- Parameters & data statement for local calculations
!
      REAL, PARAMETER :: C1=1./3., DMR1=.1E-3, DMR2=.2E-3, DMR3=.32E-3, &
     & N0r0=8.E6, N0s0=4.E6, RHOL=1000., RHOS=100.,                     &
     & XMR1=1.e6*DMR1, XMR2=1.e6*DMR2, XMR3=1.e6*DMR3
      INTEGER, PARAMETER :: MDR1=XMR1, MDR2=XMR2, MDR3=XMR3
!
      real dtph, bbfr
      integer i
!
!--- Added on 5/16/01 for Moorthi
!
      logical, parameter :: read_lookup=.false., write_lookup=.false.
!
!------------------------------------------------------------------------
!  *************  Parameters used in ETA model -- Not used in Global Model *****
!
!--- DPHD, DLMD are delta latitude and longitude at the model (NOT geodetic) equator
!    => "DX" is the hypotenuse of the model zonal & meridional grid increments.
!
!     DX=111.*(DPHD**2+DLMD**2)**.5         ! Resolution at MODEL equator (km)
!     DX=MIN(100., MAX(5., DX) )
!
!--- Assume the following functional relationship for key constants that
!    depend on grid resolution from DXmin (5 km) to DXmax (100 km) resolution:
!
!     DXmin=5.
!     DXmax=100.
!     DX=MIN(DXmax, MAX(DXmin, DX) )
!
!--- EXtune determines the degree to which the coefficients change with resolution.
!    The larger EXtune is, the more sensitive the parameter.
!
!     EXtune=1.

!
!--- FXtune ==> F(DX) is the grid-resolution tuning parameter (from 0 to 1)
!
!     FXtune=((DXmax-DX)/(DXmax-DXmin))**EXtune
!     FXtune=MAX(0., MIN(1., FXtune))
!
!--- Calculate grid-averaged RH for the onset of condensation (RHgrd) based on
!    simple ***ASSUMED*** (user-specified) values at DXmax and at DXmin.
!
!     RH_DXmax=.90              !-- 90% RH at DXmax=100 km
!     RH_DXmin=.98              !-- 98% RH at DXmin=5 km
!
!--- Note that RHgrd is right now fixed throughout the domain!!
!
!     RHgrd=RH_DXmax+(RH_DXmin-RH_DXmax)*FXtune
!   ********************************************************************************
!
!
      if (first) then
!
!--- Read in various lookup tables
!
      if ( read_lookup ) then
        OPEN (UNIT=1,FILE="eta_micro_lookup.dat",FORM="UNFORMATTED")
        READ(1) VENTR1
        READ(1) VENTR2
        READ(1) ACCRR
        READ(1) MASSR
        READ(1) VRAIN
        READ(1) RRATE
        READ(1) VENTI1
        READ(1) VENTI2
        READ(1) ACCRI
        READ(1) MASSI
        READ(1) VSNOWI
        READ(1) VEL_RF
!       read(1) my_growth    ! Applicable only for DTPH=180 s for offline testing
        CLOSE (1)
      else
        CALL ICE_LOOKUP                   ! Lookup tables for ice
        CALL RAIN_LOOKUP                  ! Lookup tables for rain
        if (write_lookup) then
          open(unit=1,file='micro_lookup.dat',form='unformatted')
          write(1) ventr1
          write(1) ventr2
          write(1) accrr
          write(1) massr
          write(1) vrain
          write(1) rrate
          write(1) venti1
          write(1) venti2
          write(1) accri
          write(1) massi
          write(1) vsnowi
          write(1) vel_rf
!         write(1) my_growth    ! Applicable only for DTPH=180 s ????
          CLOSE (1)
        endif
      endif
!!
!--- Constants associated with Biggs (1953) freezing of rain, as parameterized
!    following Lin et al. (JCAM, 1983) & Reisner et al. (1998, QJRMS).
!
      ABFR=-0.66
      BBFR=100.
      CBFR=20.*PI*PI*BBFR*RHOL*1.E-21
!
!--- QAUT0 is the threshold cloud content for autoconversion to rain
!      needed for droplets to reach a diameter of 20 microns (following
!      Manton and Cotton, 1977; Banta and Hanson, 1987, JCAM).  It is
!      **STRONGLY** affected by the assumed droplet number concentrations
!     XNCW!  For example, QAUT0=1.2567, 0.8378, or 0.4189 g/m**3 for
!     droplet number concentrations of 300, 200, and 100 cm**-3, respectively.
!
!--- Calculate grid-averaged XNCW based on simple ***ASSUMED*** (user-specified)
!    values at DXmax and at DXmin.
!
!     XNCW_DXmax=50.E6          !--  50 /cm**3 at DXmax=100 km
!     XNCW_DXmin=200.E6         !-- 200 /cm**3 at DXmin=5 km
!
!--- Note that XNCW is right now fixed throughout the domain!!
!
!     XNCW=XNCW_DXmax+(XNCW_DXmin-XNCW_DXmax)*FXtune
!
!     QAUT0=PI*RHOL*XNCW*(20.E-6)**3/6.
      QAUTx=PI*RHOL*1.0E6*(20.E-6)**3/6.
!
!--- Based on rain lookup tables for mean diameters from 0.05 to 0.45 mm
!    * Four different functional relationships of mean drop diameter as
!      a function of rain rate (RR), derived based on simple fits to
!      mass-weighted fall speeds of rain as functions of mean diameter
!      from the lookup tables.
!
      RR_DRmin=N0r0*RRATE(MDRmin)     ! RR for mean drop diameter of .05 mm
      RR_DR1=N0r0*RRATE(MDR1)         ! RR for mean drop diameter of .10 mm
      RR_DR2=N0r0*RRATE(MDR2)         ! RR for mean drop diameter of .20 mm
      RR_DR3=N0r0*RRATE(MDR3)         ! RR for mean drop diameter of .32 mm
      RR_DRmax=N0r0*RRATE(MDRmax)     ! RR for mean drop diameter of .45 mm
!
      RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
      RQR_DR1=N0r0*MASSR(MDR1)        ! Rain content for mean drop diameter of .10 mm
      RQR_DR2=N0r0*MASSR(MDR2)        ! Rain content for mean drop diameter of .20 mm
      RQR_DR3=N0r0*MASSR(MDR3)        ! Rain content for mean drop diameter of .32 mm
      RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
      C_N0r0=PI*RHOL*N0r0
      CN0r0=1.E6/C_N0r0**.25
      CN0r_DMRmin=1./(PI*RHOL*DMRmin**4)
      CN0r_DMRmax=1./(PI*RHOL*DMRmax**4)
!
      endif                     !  If (first) then loop ends here
!
!     Find out what microphysics time step should be
!
      mic_step = max(1, int(dtpg/600.0+0.5))
!     mic_step = max(1, int(dtpg/300.0+0.5))
      dtph     = dtpg / mic_step
      if (mype == 0) print *,' DTPG=',DTPG,' mic_step=',mic_step        &
     &,                ' dtph=',dtph
!
!--- Calculates coefficients for growth rates of ice nucleated in water
!    saturated conditions, scaled by physics time step (lookup table)
!
      CALL MY_GROWTH_RATES (DTPH)
!
!--- CIACW is used in calculating riming rates
!      The assumed effective collection efficiency of cloud water rimed onto
!      ice is =0.5 below:
!
!Moor CIACW=DTPH*0.25*PI*0.5*(1.E5)**C1   ! commented on 20050422
!      ice is =0.1 below:
      CIACW=DTPH*0.25*PI*0.1*(1.E5)**C1
!     CIACW = 0.0      ! Brad's suggestion 20040614
!
!--- CIACR is used in calculating freezing of rain colliding with large ice
!      The assumed collection efficiency is 1.0
!
      CIACR=PI*DTPH
!
!--- CRACW is used in calculating collection of cloud water by rain (an
!      assumed collection efficiency of 1.0)
!
!Moor CRACW=DTPH*0.25*PI*1.0                 ! commented on 20050422
!
!      assumed collection efficiency of 0.1)
      CRACW=DTPH*0.25*PI*0.1
!     CRACW = 0.0      ! Brad's suggestion 20040614
!
      ESW0=FPVSL(T0C)           ! Saturation vapor pressure at 0C
      RFmax=1.1**Nrime          ! Maximum rime factor allowed
!
!------------------------------------------------------------------------
!--------------- Constants passed through argument list -----------------
!------------------------------------------------------------------------
!
!--- Important parameters for self collection (autoconversion) of
!    cloud water to rain.
!
!--- CRAUT is proportional to the rate that cloud water is converted by
!      self collection to rain (autoconversion rate)
!
      CRAUT=1.-EXP(-1.E-3*DTPH)
!
!     IF (MYPE == 0)
!    & WRITE(6,"(A, A,F6.2,A, A,F5.4, A,F7.3,A, A,F6.2,A, A,F6.3,A)")
!    &   'KEY MICROPHYSICAL PARAMETERS FOR '
!    &  ,'DX=',DX,' KM:'
!    &  ,'   FXtune=',FXtune
!    &  ,'   RHgrd=',100.*RHgrd,' %'
!    &  ,'   NCW=',1.E-6*XNCW,' /cm**3'
!    &  ,'   QAUT0=',1.E3*QAUT0,' g/kg'
!
!--- For calculating snow optical depths by considering bulk density of
!      snow based on emails from Q. Fu (6/27-28/01), where optical
!      depth (T) = 1.5*SWP/(Reff*DENS), SWP is snow water path, Reff
!      is effective radius, and DENS is the bulk density of snow.
!
!    SWP (kg/m**2)=(1.E-3 kg/g)*SWPrad, SWPrad in g/m**2 used in radiation
!    T = 1.5*1.E3*SWPrad/(Reff*DENS)
!
!    See derivation for MASSI(INDEXS), note equal to RHO*QSNOW/NSNOW
!
!    SDENS=1.5e3/DENS, DENS=MASSI(INDEXS)/[PI*(1.E-6*INDEXS)**3]
!
      DO I=MDImin,MDImax
!MoorthiSDENS(I)=PI*1.5E-15*FLOAT(I*I*I)/MASSI(I)
        SDENS(I)=PI*1.0E-15*FLOAT(I*I*I)/MASSI(I)
      ENDDO
!
!-----------------------------------------------------------------------
!
      END subroutine gsmconst

!
!#######################################################################
!--- Sets up lookup table for calculating initial ice crystal growth ---
!#######################################################################
!
      SUBROUTINE MY_GROWTH_RATES (DTPH)
!
      implicit none
!
!--- Below are tabulated values for the predicted mass of ice crystals
!    after 600 s of growth in water saturated conditions, based on 
!    calculations from Miller and Young (JAS, 1979).  These values are
!    crudely estimated from tabulated curves at 600 s from Fig. 6.9 of
!    Young (1993).  Values at temperatures colder than -27C were 
!    assumed to be invariant with temperature.  
!
!--- Used to normalize Miller & Young (1979) calculations of ice growth
!    over large time steps using their tabulated values at 600 s.
!    Assumes 3D growth with time**1.5 following eq. (6.3) in Young (1993).
!
      real dtph, dt_ice
      REAL MY_600(MY_T1:MY_T2)
!
!-- 20090714: These values are in g and need to be converted to kg below
      DATA MY_600 /                                                     &
     & 5.5e-8,  1.4E-7,  2.8E-7, 6.E-7,   3.3E-6,                       & !  -1 to  -5 deg C
     & 2.E-6,   9.E-7,   8.8E-7, 8.2E-7,  9.4e-7,                       & !  -6 to -10 deg C
     & 1.2E-6,  1.85E-6, 5.5E-6, 1.5E-5,  1.7E-5,                       & ! -11 to -15 deg C
     & 1.5E-5,  1.E-5,   3.4E-6, 1.85E-6, 1.35E-6,                      & ! -16 to -20 deg C
     & 1.05E-6, 1.E-6,   9.5E-7, 9.0E-7 , 9.5E-7,                       &  ! -21 to -25 deg C
     & 9.5E-7,  9.E-7,   9.E-7,  9.E-7,   9.E-7,                        &  ! -26 to -30 deg C
     & 9.E-7,   9.E-7,   9.E-7,  9.E-7,   9.E-7 /                         ! -31 to -35 deg C
!
!-----------------------------------------------------------------------
!
      DT_ICE=(DTPH/600.)**1.5
!     MY_GROWTH=DT_ICE*MY_600          ! original version
      MY_GROWTH=DT_ICE*MY_600*1.E-3    !-- 20090714: Convert from g to kg
!
!-----------------------------------------------------------------------
!
      END subroutine MY_GROWTH_RATES
!
!#######################################################################
!--------------- Creates lookup tables for ice processes ---------------
!#######################################################################
!
      subroutine ice_lookup
!
      implicit none
!-----------------------------------------------------------------------------------
!
!---- Key diameter values in mm
!
!-----------------------------------------------------------------------------------
!
!---- Key concepts:
!       - Actual physical diameter of particles (D)
!       - Ratio of actual particle diameters to mean diameter (x=D/MD)
!       - Mean diameter of exponentially distributed particles, which is the
!         same as 1./LAMDA of the distribution (MD)
!       - All quantitative relationships relating ice particle characteristics as
!         functions of their diameter (e.g., ventilation coefficients, normalized
!         accretion rates, ice content, and mass-weighted fall speeds) are a result
!         of using composite relationships for ice crystals smaller than 1.5 mm
!         diameter merged with analogous relationships for larger sized aggregates.
!         Relationships are derived as functions of mean ice particle sizes assuming
!         exponential size spectra and assuming the properties of ice crystals at
!         sizes smaller than 1.5 mm and aggregates at larger sizes.  
!
!-----------------------------------------------------------------------------------
!
!---- Actual ice particle diameters for which integrated distributions are derived
!       - DminI - minimum diameter for integration (.02 mm, 20 microns)
!       - DmaxI - maximum diameter for integration (2 cm)
!       - DdelI - interval for integration (1 micron)
!
      real, parameter :: DminI=.02e-3, DmaxI=20.e-3, DdelI=1.e-6,       &
     &  XImin=1.e6*DminI, XImax=1.e6*DmaxI
      integer, parameter :: IDImin=XImin, IDImax=XImax
!
!---- Meaning of the following arrays:
!        - diam - ice particle diameter (m)
!        - mass - ice particle mass (kg)
!        - vel  - ice particle fall speeds (m/s)
!        - vent1 - 1st term in ice particle ventilation factor
!        - vent2 - 2nd term in ice particle ventilation factor
!
      real diam(IDImin:IDImax),mass(IDImin:IDImax),vel(IDImin:IDImax),  &
     & vent1(IDImin:IDImax),vent2(IDImin:IDImax)
!
!-----------------------------------------------------------------------------------
!
!---- Found from trial & error that the m(D) & V(D) mass & velocity relationships
!       between the ice crystals and aggregates overlapped & merged near a particle
!       diameter sizes of 1.5 mm.  Thus, ice crystal relationships are used for
!       sizes smaller than 1.5 mm and aggregate relationships for larger sizes.
!
      real, parameter :: d_crystal_max=1.5
!
!---- The quantity xmax represents the maximum value of "x" in which the
!       integrated values are calculated.  For xmax=20., this means that
!       integrated ventilation, accretion, mass, and precipitation rates are
!       calculated for ice particle sizes less than 20.*mdiam, the mean particle diameter.
!
      real, parameter :: xmax=20.
!
!-----------------------------------------------------------------------------------
!
!---- Meaning of the following arrays:
!        - mdiam - mean diameter (m)
!        - VENTI1 - integrated quantity associated w/ ventilation effects
!                   (capacitance only) for calculating vapor deposition onto ice
!        - VENTI2 - integrated quantity associated w/ ventilation effects
!                   (with fall speed) for calculating vapor deposition onto ice
!        - ACCRI  - integrated quantity associated w/ cloud water collection by ice
!        - MASSI  - integrated quantity associated w/ ice mass 
!        - VSNOWI - mass-weighted fall speed of snow, used to calculate precip rates
!
!--- Mean ice-particle diameters varying from 50 microns to 1000 microns (1 mm), 
!      assuming an exponential size distribution.  
!
      real mdiam
!
!-----------------------------------------------------------------------------------
!------------- Constants & parameters for ventilation factors of ice ---------------
!-----------------------------------------------------------------------------------
!
!---- These parameters are used for calculating the ventilation factors for ice
!       crystals between 0.2 and 1.5 mm diameter (Hall and Pruppacher, JAS, 1976).  
!       From trial & error calculations, it was determined that the ventilation
!       factors of smaller ice crystals could be approximated by a simple linear
!       increase in the ventilation coefficient from 1.0 at 50 microns (.05 mm) to 
!       1.1 at 200 microns (0.2 mm), rather than using the more complex function of
!       1.0 + .14*(Sc**.33*Re**.5)**2 recommended by Hall & Pruppacher.
!
      real, parameter :: cvent1i=.86, cvent2i=.28
!
!---- These parameters are used for calculating the ventilation factors for larger
!       aggregates, where D>=1.5 mm (see Rutledge and Hobbs, JAS, 1983; 
!       Thorpe and Mason, 1966).
!
      real, parameter :: cvent1a=.65, cvent2a=.44
!
      real m_agg,m_bullet,m_column,m_ice,m_plate
!
!---- Various constants
!
      real, parameter :: c1=2./3., cexp=1./3.
!
      logical :: iprint
      logical, parameter :: print_diag=.false.
!
!-----------------------------------------------------------------------------------
!- Constants & parameters for calculating the increase in fall speed of rimed ice --
!-----------------------------------------------------------------------------------
!
!---- Constants & arrays for estimating increasing fall speeds of rimed ice.
!     Based largely on theory and results from Bohm (JAS, 1989, 2419-2427).
!
!-------------------- Standard atmosphere conditions at 1000 mb --------------------
!
      real, parameter :: t_std=288., dens_std=1000.e2/(287.04*288.)
!
!---- These "bulk densities" are the actual ice densities in the ice portion of the 
!     lattice.  They are based on text associated w/ (12) on p. 2425 of Bohm (JAS, 
!     1989).  Columns, plates, & bullets are assumed to have an average bulk density 
!     of 850 kg/m**3.  Aggregates were assumed to have a slightly larger bulk density 
!     of 600 kg/m**3 compared with dendrites (i.e., the least dense, most "lacy" & 
!     tenous ice crystal, which was assumed to be ~500 kg/m**3 in Bohm).  
!
      real, parameter :: dens_crystal=850., dens_agg=600.
!
!--- A value of Nrime=40 for a logarithmic ratio of 1.1 yields a maximum rime factor
!      of 1.1**40 = 45.26 that is resolved in these tables.  This allows the largest
!      ice particles with a mean diameter of MDImax=1000 microns to achieve bulk 
!      densities of 900 kg/m**3 for rimed ice.  
!
!     integer, parameter :: Nrime=40
      real m_rime,                                                      &
     &     rime_factor(0:Nrime), rime_vel(0:Nrime),                     &
     &     vel_rime(IDImin:IDImax,Nrime), ivel_rime(MDImin:MDImax,Nrime)
!
      integer i, j, jj, k, icount
      real c2,      cbulk, cbulk_ice, px, dynvis_std, crime1            &
     &,    crime2,  crime3, crime4, crime5, d, c_avg, c_agg             &
     &,    c_bullet, c_column, c_plate, cl_agg, cl_bullet               &
     &,    cl_column, cl_plate, v_agg, v_bullet, v_column               &
     &,    v_plate,   wd,       ecc_column                              &
     &,    cvent1,    cvent2, crime_best, rime_m1, rime_m2              &
     &,    x_rime,    re_rime, smom3, pratei, expf                      &
     &,    bulk_dens, xmass,  xmdiam, ecc_plate, dx
!
!-----------------------------------------------------------------------------------
!----------------------------- BEGIN EXECUTION -------------------------------------
!-----------------------------------------------------------------------------------
!
!
      c2=1./sqrt(3.)
!     pi=acos(-1.)
      cbulk=6./pi
      cbulk_ice=900.*pi/6.    ! Maximum bulk ice density allowed of 900 kg/m**3
      px=.4**cexp             ! Convert fall speeds from 400 mb (Starr & Cox) to 1000 mb
!
!--------------------- Dynamic viscosity (1000 mb, 288 K) --------------------------
!
      dynvis_std=1.496e-6*t_std**1.5/(t_std+120.)
      crime1=pi/24.
      crime2=8.*9.81*dens_std/(pi*dynvis_std**2)
      crime3=crime1*dens_crystal
      crime4=crime1*dens_agg
      crime5=dynvis_std/dens_std
      do i=0,Nrime
        rime_factor(i)=1.1**i
      enddo
!
!#######################################################################
!      Characteristics as functions of actual ice particle diameter 
!#######################################################################
!
!----   M(D) & V(D) for 3 categories of ice crystals described by Starr 
!----   & Cox (1985). 
!
!----   Capacitance & characteristic lengths for Reynolds Number calculations
!----   are based on Young (1993; p. 144 & p. 150).  c-axis & a-axis 
!----   relationships are from Heymsfield (JAS, 1972; Table 1, p. 1351).
!
      icount=60
!
      if (print_diag)                                                   & 
     &  write(7,"(2a)") '---- Increase in fall speeds of rimed ice',    &
     &    ' particles as function of ice particle diameter ----'
      do i=IDImin,IDImax
        if (icount == 60 .and. print_diag) then
          write(6,"(/2a/3a)") 'Particle masses (mg), fall speeds ',     &
     &      '(m/s), and ventilation factors',                           &
     &      '  D(mm)  CR_mass   Mass_bull   Mass_col  Mass_plat ',      &
     &      '  Mass_agg   CR_vel  V_bul CR_col CR_pla Aggreg',          &
     &      '    Vent1      Vent2 '                               
          write(7,"(3a)") '        <----------------------------------',&
     &      '---------------  Rime Factor  --------------------------', &
     &      '--------------------------->'
          write(7,"(a,23f5.2)") '  D(mm)',(rime_factor(k), k=1,5),      &
     &       (rime_factor(k), k=6,40,2)
          icount=0
        endif
        d=(float(i)+.5)*1.e-3    ! in mm
        c_avg=0.
        c_agg=0.
        c_bullet=0.
        c_column=0.
        c_plate=0.
        cl_agg=0.
        cl_bullet=0.
        cl_column=0.
        cl_plate=0.
        m_agg=0.
        m_bullet=0.
        m_column=0.
        m_plate=0.
        v_agg=0.
        v_bullet=0.
        v_column=0.
        v_plate=0.
        if (d < d_crystal_max) then
!
!---- This block of code calculates bulk characteristics based on average
!     characteristics of bullets, plates, & column ice crystals <1.5 mm size
!
!---- Mass-diameter relationships from Heymsfield (1972) & used
!       in Starr & Cox (1985), units in mg
!---- "d" is maximum dimension size of crystal in mm, 
!
! Mass of pure ice for spherical particles, used as an upper limit for the
!   mass of small columns (<~ 80 microns) & plates (<~ 35 microns)
!
          m_ice=.48*d**3   ! Mass of pure ice for spherical particle
!
          m_bullet=min(.044*d**3, m_ice)
          m_column=min(.017*d**1.7, m_ice)
          m_plate=min(.026*d**2.5, m_ice)
!
          mass(i)=m_bullet+m_column+m_plate
!
!---- These relationships are from Starr & Cox (1985), applicable at 400 mb
!---- "d" is maximum dimension size of crystal in mm, dx in microns
!
          dx=1000.*d            ! Convert from mm to microns
          if (dx <= 200.) then
            v_column=8.114e-5*dx**1.585
            v_bullet=5.666e-5*dx**1.663
            v_plate=1.e-3*dx
          else if (dx <= 400.) then
            v_column=4.995e-3*dx**.807
            v_bullet=3.197e-3*dx**.902
            v_plate=1.48e-3*dx**.926
          else if (dx <= 600.) then
            v_column=2.223e-2*dx**.558
            v_bullet=2.977e-2*dx**.529
            v_plate=9.5e-4*dx
          else if (dx <= 800.) then
            v_column=4.352e-2*dx**.453
            v_bullet=2.144e-2*dx**.581
            v_plate=3.161e-3*dx**.812
          else 
            v_column=3.833e-2*dx**.472
            v_bullet=3.948e-2*dx**.489
            v_plate=7.109e-3*dx**.691
          endif
!
!---- Reduce fall speeds from 400 mb to 1000 mb
!
          v_column=px*v_column
          v_bullet=px*v_bullet
          v_plate=px*v_plate
!
!---- DIFFERENT VERSION!  CALCULATES MASS-WEIGHTED CRYSTAL FALL SPEEDS
!
          vel(i)=(m_bullet*v_bullet+m_column*v_column+m_plate*v_plate)/ &
     &           mass(i)
          mass(i)=mass(i)/3.
!
!---- Shape factor and characteristic length of various ice habits,
!     capacitance is equal to 4*PI*(Shape factor)
!       See Young (1993, pp. 143-152 for guidance)
!
!---- Bullets:
!
!---- Shape factor for bullets (Heymsfield, 1975)
          c_bullet=.5*d
!---- Length-width functions for bullets from Heymsfield (JAS, 1972)
          if (d > 0.3) then
            wd=.25*d**.7856     ! Width (mm); a-axis
          else
            wd=.185*d**.552
          endif
!---- Characteristic length for bullets (see first multiplicative term on right
!       side of eq. 7 multiplied by crystal width on p. 821 of Heymsfield, 1975)
          cl_bullet=.5*pi*wd*(.25*wd+d)/(d+wd)
!
!---- Plates:
!
!---- Length-width function for plates from Heymsfield (JAS, 1972)
          wd=.0449*d**.449      ! Width or thickness (mm); c-axis
!---- Eccentricity & shape factor for thick plates following Young (1993, p. 144)
          ecc_plate=sqrt(1.-wd*wd/(d*d))         ! Eccentricity
          c_plate=d*ecc_plate/asin(ecc_plate)    ! Shape factor
!---- Characteristic length for plates following Young (1993, p. 150, eq. 6.6)
          cl_plate=d+2.*wd      ! Characteristic lengths for plates
!
!---- Columns:
!
!---- Length-width function for columns from Heymsfield (JAS, 1972)
          if (d > 0.2) then
            wd=.1973*d**.414    ! Width (mm); a-axis
          else
            wd=.5*d             ! Width (mm); a-axis
          endif
!---- Eccentricity & shape factor for columns following Young (1993, p. 144)
          ecc_column=sqrt(1.-wd*wd/(d*d))                     ! Eccentricity
          c_column=ecc_column*d/log((1.+ecc_column)*d/wd)     ! Shape factor
!---- Characteristic length for columns following Young (1993, p. 150, eq. 6.7)
          cl_column=(wd+2.*d)/(c1+c2*d/wd)       ! Characteristic lengths for columns
!
!---- Convert shape factor & characteristic lengths from mm to m for 
!       ventilation calculations
!
          c_bullet=.001*c_bullet
          c_plate=.001*c_plate
          c_column=.001*c_column
          cl_bullet=.001*cl_bullet
          cl_plate=.001*cl_plate
          cl_column=.001*cl_column
!
!---- Make a smooth transition between a ventilation coefficient of 1.0 at 50 microns
!       to 1.1 at 200 microns
!
          if (d > 0.2) then
            cvent1=cvent1i
            cvent2=cvent2i/3.
          else
            cvent1=1.0+.1*max(0., d-.05)/.15
            cvent2=0.
          endif
!
!---- Ventilation factors for ice crystals:
!
          vent1(i)=cvent1*(c_bullet+c_plate+c_column)/3.
          vent2(i)=cvent2*(c_bullet*sqrt(v_bullet*cl_bullet)            &
     &                    +c_plate*sqrt(v_plate*cl_plate)               &
     &                    +c_column*sqrt(v_column*cl_column) )
          crime_best=crime3     ! For calculating Best No. of rimed ice crystals
        else
!
!---- This block of code calculates bulk characteristics based on average
!     characteristics of unrimed aggregates >= 1.5 mm using Locatelli & 
!     Hobbs (JGR, 1974, 2185-2197) data.
!
!----- This category is a composite of aggregates of unrimed radiating 
!-----   assemblages of dendrites or dendrites; aggregates of unrimed
!-----   radiating assemblages of plates, side planes, bullets, & columns;
!-----   aggregates of unrimed side planes (mass in mg, velocity in m/s)
!
          m_agg=(.073*d**1.4+.037*d**1.9+.04*d**1.4)/3.
          v_agg=(.8*d**.16+.69*d**.41+.82*d**.12)/3.
          mass(i)=m_agg
          vel(i)=v_agg
!
!---- Assume spherical aggregates
!
!---- Shape factor is the same as for bullets, = D/2
          c_agg=.001*.5*d         ! Units of m
!---- Characteristic length is surface area divided by perimeter
!       (.25*PI*D**2)/(PI*D**2) = D/4
          cl_agg=.5*c_agg         ! Units of m
!
!---- Ventilation factors for aggregates:
!
          vent1(i)=cvent1a*c_agg
          vent2(i)=cvent2a*c_agg*sqrt(v_agg*cl_agg)
          crime_best=crime4     ! For calculating Best No. of rimed aggregates
        endif
!
!---- Convert from shape factor to capacitance for ventilation factors
!
        vent1(i)=4.*pi*vent1(i)
        vent2(i)=4.*pi*vent2(i)
        diam(i)=1.e-3*d             ! Convert from mm to m
        mass(i)=1.e-6*mass(i)       ! Convert from mg to kg
!
!---- Calculate increase in fall speeds of individual rimed ice particles
!
        do k=0,Nrime
!---- Mass of rimed ice particle associated with rime_factor(k)
          rime_m1=rime_factor(k)*mass(i)
          rime_m2=cbulk_ice*diam(i)**3
          m_rime=min(rime_m1, rime_m2)
!---- Best Number (X) of rimed ice particle combining eqs. (8) & (12) in Bohm
          x_rime=crime2*m_rime*(crime_best/m_rime)**.25
!---- Reynolds Number for rimed ice particle using eq. (11) in Bohm
          re_rime=8.5*(sqrt(1.+.1519*sqrt(x_rime))-1.)**2
          rime_vel(k)=crime5*re_rime/diam(i)
        enddo
        do k=1,Nrime
          vel_rime(i,k)=rime_vel(k)/rime_vel(0)
        enddo
        if (print_diag) then
   !
   !---- Determine if statistics should be printed out.
   !
          iprint=.false.
          if (d <= 1.) then
            if (mod(i,10) == 0) iprint=.true.
          else
            if (mod(i,100) == 0) iprint=.true.
          endif
          if (iprint) then
            write(6,"(f7.4,5e11.4,1x,5f7.4,1x,2e11.4)")                 & 
     &        d,1.e6*mass(i),m_bullet,m_column,m_plate,m_agg,           &
     &        vel(i),v_bullet,v_column,v_plate,v_agg,                   &
     &        vent1(i),vent2(i)
            write(7,"(f7.4,23f5.2)") d,(vel_rime(i,k), k=1,5),          &
     &        (vel_rime(i,k), k=6,40,2)
            icount=icount+1
          endif
        endif
      enddo
!
!#######################################################################
!      Characteristics as functions of mean particle diameter
!#######################################################################
!
      VENTI1=0.
      VENTI2=0.
      ACCRI=0.
      MASSI=0.
      VSNOWI=0.
      VEL_RF=0.
      ivel_rime=0.
      icount=0
      if (print_diag) then
        icount=60
        write(6,"(/2a)") '------------- Statistics as functions of ',   &
     &    'mean particle diameter -------------'
        write(7,"(/2a)") '------ Increase in fall speeds of rimed ice', &
     &    ' particles as functions of mean particle diameter -----'
      endif
      do j=MDImin,MDImax
        if (icount == 60 .AND. print_diag) then
          write(6,"(/2a)") 'D(mm)    Vent1      Vent2    ',             &
     &       'Accrete       Mass     Vel  Dens  '
          write(7,"(3a)") '      <----------------------------------',  &
     &      '---------------  Rime Factor  --------------------------', &
     &      '--------------------------->'
          write(7,"(a,23f5.2)") 'D(mm)',(rime_factor(k), k=1,5),        &
     &       (rime_factor(k), k=6,40,2)
          icount=0
        endif
        mdiam=DelDMI*float(j)       ! in m
        smom3=0.
        pratei=0.
        rime_vel=0.                 ! Note that this array is being reused!
        do i=IDImin,IDImax
          dx=diam(i)/mdiam
          if (dx <= xmax) then      ! To prevent arithmetic underflows
            expf=exp(-dx)*DdelI
            VENTI1(J)=VENTI1(J)+vent1(i)*expf
            VENTI2(J)=VENTI2(J)+vent2(i)*expf
            ACCRI(J)=ACCRI(J)+diam(i)*diam(i)*vel(i)*expf
            xmass=mass(i)*expf
            do k=1,Nrime
              rime_vel(k)=rime_vel(k)+xmass*vel_rime(i,k)
            enddo
            MASSI(J)=MASSI(J)+xmass
            pratei=pratei+xmass*vel(i)
            smom3=smom3+diam(i)**3*expf
          else
            exit
          endif
        enddo
   !
   !--- Increased fall velocities functions of mean diameter (j),
   !      normalized by ice content, and rime factor (k) 
   !
        do k=1,Nrime
          ivel_rime(j,k)=rime_vel(k)/MASSI(J)
        enddo
   !
   !--- Increased fall velocities functions of ice content at 0.1 mm
   !      intervals (j_100) and rime factor (k); accumulations here
   !
        jj=j/100
        if (jj >= 2 .AND. jj <= 9) then
          do k=1,Nrime
            VEL_RF(jj,k)=VEL_RF(jj,k)+ivel_rime(j,k)
          enddo
        endif
        bulk_dens=cbulk*MASSI(J)/smom3
        VENTI1(J)=VENTI1(J)/mdiam
        VENTI2(J)=VENTI2(J)/mdiam
        ACCRI(J)=ACCRI(J)/mdiam
        VSNOWI(J)=pratei/MASSI(J)
        MASSI(J)=MASSI(J)/mdiam
        if (mod(j,10) == 0 .AND. print_diag) then
          xmdiam=1.e3*mdiam
          write(6,"(f6.3,4e11.4,f6.3,f8.3)") xmdiam,VENTI1(j),VENTI2(j),&
     &      ACCRI(j),MASSI(j),VSNOWI(j),bulk_dens
          write(7,"(f6.3,23f5.2)") xmdiam,(ivel_rime(j,k), k=1,5),      &
     &       (ivel_rime(j,k), k=6,40,2)
          icount=icount+1
        endif
      enddo
!
!--- Average increase in fall velocities rimed ice as functions of mean
!      particle diameter (j, only need 0.1 mm intervals) and rime factor (k)
!
      if (print_diag) then
        write(7,"(/2a)") ' ------- Increase in fall speeds of rimed ',  &
     &    'ice particles at reduced, 0.1-mm intervals  --------'
        write(7,"(3a)") '        <----------------------------------',  &
     &    '---------------  Rime Factor  --------------------------',   &
     &    '--------------------------->'
        write(7,"(a,23f5.2)") 'D(mm)',(rime_factor(k), k=1,5),          &
     &    (rime_factor(k), k=6,40,2)
      endif
      do j=2,9
        VEL_RF(j,0)=1.
        do k=1,Nrime
          VEL_RF(j,k)=.01*VEL_RF(j,k)
        enddo
        if (print_diag) write(7,"(f4.1,2x,23f5.2)") 0.1*j,              &
     &    (VEL_RF(j,k), k=1,5),(VEL_RF(j,k), k=6,40,2)
      enddo
!
!-----------------------------------------------------------------------------------
!
      end subroutine ice_lookup
!
!#######################################################################
!-------------- Creates lookup tables for rain processes ---------------
!#######################################################################
!
      subroutine rain_lookup
      implicit none
!
!--- Parameters & arrays for fall speeds of rain as a function of rain drop
!      diameter.  These quantities are integrated over exponential size
!      distributions of rain drops at 1 micron intervals (DdelR) from minimum 
!      drop sizes of .05 mm (50 microns, DminR) to maximum drop sizes of 10 mm 
!      (DmaxR). 
!
      real, parameter :: DminR=.05e-3, DmaxR=10.e-3, DdelR=1.e-6,       &
     & XRmin=1.e6*DminR, XRmax=1.e6*DmaxR
      integer, parameter :: IDRmin=XRmin, IDRmax=XRmax
      real diam(IDRmin:IDRmax), vel(IDRmin:IDRmax)
!
!--- Parameters rain lookup tables, which establish the range of mean drop
!      diameters; from a minimum mean diameter of 0.05 mm (DMRmin) to a 
!      maximum mean diameter of 0.45 mm (DMRmax).  The tables store solutions
!      at 1 micron intervals (DelDMR) of mean drop diameter.  
!
      real mdiam, mass
!
      logical, parameter :: print_diag=.false.
!
      real d, cmass, pi2, expf
      integer i, j, i1, i2
!
!-----------------------------------------------------------------------
!------- Fall speeds of rain as function of rain drop diameter ---------
!-----------------------------------------------------------------------
!
      do i=IDRmin,IDRmax
        diam(i)=float(i)*DdelR
        d=100.*diam(i)         ! Diameter in cm
        if (d <= .42) then
   !
   !--- Rutledge & Hobbs (1983); vel (m/s), d (cm)
   !
          vel(i)=max(0., -.267+51.5*d-102.25*d*d+75.7*d**3)
        else if (d > 0.42 .and. d <= .58) then
   !
   !--- Linear interpolation of Gunn & Kinzer (1949) data
   !
          vel(i)=8.92+.25/(.58-.42)*(d-.42)
        else
          vel(i)=9.17
        endif
      enddo
      do i=1,100
        i1=(i-1)*100+IDRmin
        i2=i1+90
   !
   !--- Print out rain fall speeds only for D<=5.8 mm (.58 cm)
   !
        if (diam(i1) > .58e-2) exit
        if (print_diag) then
          write(6,"(1x)")
          write(6,"('D(mm)->  ',10f7.3)") (1000.*diam(j), j=i1,i2,10)
          write(6,"('V(m/s)-> ',10f7.3)") (vel(j), j=i1,i2,10)
        endif
      enddo
!
!-----------------------------------------------------------------------
!------------------- Lookup tables for rain processes ------------------
!-----------------------------------------------------------------------
!
!     pi=acos(-1.)
      pi2=2.*pi
      cmass=1000.*pi/6.
      if (print_diag) then
        write(6,"(/'Diam - Mean diameter (mm)'                          &
     &          /'VENTR1 - 1st ventilation coefficient (m**2)'          &
     &          /'VENTR2 - 2nd ventilation coefficient (m**3/s**.5)'    &
     &          /'ACCRR - accretion moment (m**4/s)'                    &
     &          /'RHO*QR - mass content (kg/m**3) for N0r=8e6'          &
     &          /'RRATE - rain rate moment (m**5/s)'                    &
     &          /'VR - mass-weighted rain fall speed (m/s)'             &
     &    /' Diam      VENTR1      VENTR2       ACCRR      ',           &
     &    'RHO*QR       RRATE    VR')")
      endif
      do j=MDRmin,MDRmax
        mdiam=float(j)*DelDMR
        VENTR2(J)=0.
        ACCRR(J)=0.
        MASSR(J)=0.
        RRATE(J)=0.
        do i=IDRmin,IDRmax
          expf=exp(-diam(i)/mdiam)*DdelR
          VENTR2(J)=VENTR2(J)+diam(i)**1.5*vel(i)**.5*expf
          ACCRR(J)=ACCRR(J)+diam(i)*diam(i)*vel(i)*expf
          MASSR(J)=MASSR(J)+diam(i)**3*expf
          RRATE(J)=RRATE(J)+diam(i)**3*vel(i)*expf
        enddo
   !
   !--- Derived based on ventilation, F(D)=0.78+.31*Schmidt**(1/3)*Reynold**.5,
   !      where Reynold=(V*D*rho/dyn_vis), V is velocity, D is particle diameter,
   !      rho is air density, & dyn_vis is dynamic viscosity.  Only terms 
   !      containing velocity & diameter are retained in these tables.  
   !
        VENTR1(J)=.78*pi2*mdiam**2
        VENTR2(J)=.31*pi2*VENTR2(J)
   !
        MASSR(J)=cmass*MASSR(J)
        RRATE(J)=cmass*RRATE(J)
        VRAIN(J)=RRATE(J)/MASSR(J)
        if (print_diag) write(6,"(f6.3,5g12.5,f6.3)") 1000.*mdiam,      &
     &    ventr1(j),ventr2(j),accrr(j),8.e6*massr(j),rrate(j),vrain(j)
      enddo
!
!-----------------------------------------------------------------------
!
      end subroutine rain_lookup
!
!###############################################################################
! ***** VERSION OF MICROPHYSICS DESIGNED FOR HIGHER RESOLUTION MESO ETA MODEL
!       (1) Represents sedimentation by preserving a portion of the precipitation
!           through top-down integration from cloud-top.  Modified procedure to
!           Zhao and Carr (1997).
!       (2) Microphysical equations are modified to be less sensitive to time
!           steps by use of Clausius-Clapeyron equation to account for changes in
!           saturation mixing ratios in response to latent heating/cooling.  
!       (3) Prevent spurious temperature oscillations across 0C due to 
!           microphysics.
!       (4) Uses lookup tables for: calculating two different ventilation 
!           coefficients in condensation and deposition processes; accretion of
!           cloud water by precipitation; precipitation mass; precipitation rate
!           (and mass-weighted precipitation fall speeds).
!       (5) Assumes temperature-dependent variation in mean diameter of large ice
!           (Houze et al., 1979; Ryan et al., 1996).
!        -> 8/22/01: This relationship has been extended to colder temperatures
!           to parameterize smaller large-ice particles down to mean sizes of MDImin,
!           which is 50 microns reached at -55.9C.
!       (6) Attempts to differentiate growth of large and small ice, mainly for
!           improved transition from thin cirrus to thick, precipitating ice
!           anvils.
!        -> 8/22/01: This feature has been diminished by effectively adjusting to
!           ice saturation during depositional growth at temperatures colder than
!           -10C.  Ice sublimation is calculated more explicitly.  The logic is
!           that sources of are either poorly understood (e.g., nucleation for NWP) 
!           or are not represented in the Eta model (e.g., detrainment of ice from 
!           convection).  Otherwise the model is too wet compared to the radiosonde
!           observations based on 1 Feb - 18 March 2001 retrospective runs.  
!       (7) Top-down integration also attempts to treat mixed-phase processes,
!           allowing a mixture of ice and water.  Based on numerous observational
!           studies, ice growth is based on nucleation at cloud top &
!           subsequent growth by vapor deposition and riming as the ice particles 
!           fall through the cloud.  Effective nucleation rates are a function
!           of ice supersaturation following Meyers et al. (JAM, 1992).  
!        -> 8/22/01: The simulated relative humidities were far too moist compared 
!           to the rawinsonde observations.  This feature has been substantially 
!           diminished, limited to a much narrower temperature range of 0 to -10C.  
!       (8) Depositional growth of newly nucleated ice is calculated for large time
!           steps using Fig. 8 of Miller and Young (JAS, 1979), at 1 deg intervals
!           using their ice crystal masses calculated after 600 s of growth in water
!           saturated conditions.  The growth rates are normalized by time step
!           assuming 3D growth with time**1.5 following eq. (6.3) in Young (1993).
!        -> 8/22/01: This feature has been effectively limited to 0 to -10C.  
!       (9) Ice precipitation rates can increase due to increase in response to
!           cloud water riming due to (a) increased density & mass of the rimed
!           ice, and (b) increased fall speeds of rimed ice.
!        -> 8/22/01: This feature has been effectively limited to 0 to -10C.  
!###############################################################################
!###############################################################################
!
      SUBROUTINE GSMCOLUMN ( ARAING, ASNOWG, DTPG, I_index, J_index,    &
     & LSFC, P_col, QI_col, QR_col, QV_col, QW_col, RimeF_col, T_col,   &
     & THICK_col, WC_col, LM, RHC_col, XNCW, FLGmin, PRINT_diag, psfc)
!
      implicit none
!
!###############################################################################
!###############################################################################
!
!-------------------------------------------------------------------------------
!----- NOTE:  In this version of the Code threading should be done outside!  
!-------------------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:  Grid-scale microphysical processes - condensation & precipitation
!   PRGRMMR: Ferrier         ORG: W/NP22     DATE: 08-2001
!   Updated: Moorthi for GFS application
!-------------------------------------------------------------------------------
! ABSTRACT:
!   * Merges original GSCOND & PRECPD subroutines.   
!   * Code has been substantially streamlined and restructured.
!   * Exchange between water vapor & small cloud condensate is calculated using
!     the original Asai (1965, J. Japan) algorithm.  See also references to
!     Yau and Austin (1979, JAS), Rutledge and Hobbs (1983, JAS), and Tao et al.
!     (1989, MWR).  This algorithm replaces the Sundqvist et al. (1989, MWR)
!     parameterization.  
!-------------------------------------------------------------------------------
!     
! USAGE: 
!   * CALL GSMCOLUMN FROM SUBROUTINE GSMDRIVE
!   * SUBROUTINE GSMDRIVE CALLED FROM MAIN PROGRAM EBU
!
! INPUT ARGUMENT LIST:
!   DTPH       - physics time step (s)
!   I_index    - I index
!   J_index    - J index
!   LSFC       - Eta level of level above surface, ground
!   P_col      - vertical column of model pressure (Pa)
!   QI_col     - vertical column of model ice mixing ratio (kg/kg)
!   QR_col     - vertical column of model rain ratio (kg/kg)
!   QV_col     - vertical column of model water vapor specific humidity (kg/kg)
!   QW_col     - vertical column of model cloud water mixing ratio (kg/kg)
!   RimeF_col  - vertical column of rime factor for ice in model (ratio, defined below)
!   T_col      - vertical column of model temperature (deg K)
!   THICK_col  - vertical column of model mass thickness (density*height increment)
!   WC_col     - vertical column of model mixing ratio of total condensate (kg/kg)
!   
!
! OUTPUT ARGUMENT LIST: 
!   ARAIN      - accumulated rainfall at the surface (kg)
!   ASNOW      - accumulated snowfall at the surface (kg)
!   QV_col     - vertical column of model water vapor specific humidity (kg/kg)
!   WC_col     - vertical column of model mixing ratio of total condensate (kg/kg)
!   QW_col     - vertical column of model cloud water mixing ratio (kg/kg)
!   QI_col     - vertical column of model ice mixing ratio (kg/kg)
!   QR_col     - vertical column of model rain ratio (kg/kg)
!   RimeF_col  - vertical column of rime factor for ice in model (ratio, defined below)
!   T_col      - vertical column of model temperature (deg K)
!     
! OUTPUT FILES:
!     NONE
!     
! Subprograms & Functions called:
!   * Real Function CONDENSE  - cloud water condensation
!   * Real Function DEPOSIT   - ice deposition (not sublimation)
!
! UNIQUE: NONE
!  
! LIBRARY: NONE
!  
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!
!------------------------------------------------------------------------- 
!--------------- Arrays & constants in argument list --------------------- 
!------------------------------------------------------------------------- 
!
      integer lm
      REAL ARAING, ASNOWG, P_col(LM), QI_col(LM), QR_col(LM), QV_col(LM)&
     &,    QW_col(LM), RimeF_col(LM), T_col(LM), THICK_col(LM),         &
     &     WC_col(LM), RHC_col(LM), XNCW(LM), ARAIN, ASNOW, dtpg, psfc
      real flgmin
!
      INTEGER I_index, J_index, LSFC
!
!
!------------------------------------------------------------------------- 
!
!--- Mean ice-particle diameters varying from 50 microns to 1000 microns
!      (1 mm), assuming an exponential size distribution.  
!
!---- Meaning of the following arrays: 
!        - mdiam - mean diameter (m)
!        - VENTI1 - integrated quantity associated w/ ventilation effects 
!                   (capacitance only) for calculating vapor deposition onto ice
!        - VENTI2 - integrated quantity associated w/ ventilation effects 
!                   (with fall speed) for calculating vapor deposition onto ice
!        - ACCRI  - integrated quantity associated w/ cloud water collection by ice
!        - MASSI  - integrated quantity associated w/ ice mass 
!        - VSNOWI - mass-weighted fall speed of snow (large ice), used to calculate 
!                   precipitation rates
!
      REAL,    PARAMETER :: DMImin=.05e-3, DMImax=1.e-3, DelDMI=1.e-6,  &
     &                      XMImin=1.e6*DMImin, XMImax=1.e6*DMImax
      INTEGER, PARAMETER :: MDImin=XMImin, MDImax=XMImax
!
!------------------------------------------------------------------------- 
!------- Key parameters, local variables, & important comments ---------
!-----------------------------------------------------------------------
!
!--- KEY Parameters:
!
!---- Comments on 14 March 2002
!    * Set EPSQ to the universal value of 1.e-12 throughout the code
!      condensate.  The value of EPSQ will need to be changed in the other 
!      subroutines in order to make it consistent throughout the Eta code.  
!    * Set CLIMIT=10.*EPSQ as the lower limit for the total mass of 
!      condensate in the current layer and the input flux of condensate
!      from above (TOT_ICE, TOT_ICEnew, TOT_RAIN, and TOT_RAINnew).
!
!-- NLImax - maximum number concentration of large ice crystals (20,000 /m**3, 20 per liter)
!-- NLImin - minimum number concentration of large ice crystals (100 /m**3, 0.1 per liter)
!
      REAL, PARAMETER ::   RHOL=1000.,  XLS=HVAP+HFUS                   &

!    &, T_ICE=-10.          !- Ver1
!    &, T_ICE_init=-5.      !- Ver1
!!!  &, T_ICE=-20.          !- Ver2
!    &, T_ICE=-40.          !- Ver2
!    &, T_ICE_init=-15.,    !- Ver2
!
!    & CLIMIT=10.*EPSQ, EPS1=RV/RD-1., RCP=1./CP,

     &,CLIMIT=10.*EPSQ, RCP=1./CP,                                      &
     & RCPRV=RCP/RV, RRHOL=1./RHOL, XLS1=XLS*RCP, XLS2=XLS*XLS*RCPRV,   &
     & XLS3=XLS*XLS/RV,                                                 &
     & C1=1./3., C2=1./6., C3=3.31/6.,                                  &
     & DMR1=.1E-3, DMR2=.2E-3, DMR3=.32E-3, N0r0=8.E6, N0rmin=1.e4,     &
     & N0s0=4.E6, RHO0=1.194, XMR1=1.e6*DMR1, XMR2=1.e6*DMR2,           &
     & XMR3=1.e6*DMR3, Xratio=.025
      INTEGER, PARAMETER :: MDR1=XMR1, MDR2=XMR2, MDR3=XMR3
!
!--- If BLEND=1:
!      precipitation (large) ice amounts are estimated at each level as a 
!      blend of ice falling from the grid point above and the precip ice
!      present at the start of the time step (see TOT_ICE below).
!--- If BLEND=0:
!      precipitation (large) ice amounts are estimated to be the precip
!      ice present at the start of the time step.
!
!--- Extended to include sedimentation of rain on 2/5/01 
!
      REAL, PARAMETER :: BLEND=1.
!
!--- This variable is for debugging purposes (if .true.)
!
      LOGICAL  PRINT_diag
!
!--- Local variables
!
      REAL    EMAIRI, N0r,         NLICE,       NSmICE, NLImax, pfac
      LOGICAL CLEAR,  ICE_logical, DBG_logical, RAIN_logical
 
      integer lbef, ipass, ixrf, ixs, itdx, idr                         &
     &,       index_my, indexr, indexr1, indexs                         &
     &,       i, j, k, l, ntimes, item
!    &,       i, j, k, my_600, i1, i2, l, ntimes

      real flimass,  xlimass, vsnow,   qi_min, dum,    piloss           &
     &,    tot_ice,  xsimass, vel_inc, vrimef, rimef1, dum1             &
     &,    dum2,     fws,     denomi,  dwv                              &
     &,    xrf,      qw0,     dli,     xli,    fsmall                   &
     &,    prevp,    tk2,     dtph                                      &
     &,    pievp,    picnd,   piacr,   pracw                            &
     &,    praut,    pimlt,   qtice,   qlice                            &
     &,    gammar,   flarge,  wvqw,    dynvis                           &
     &,    tfactor,  denom,   gammas,  diffus, therm_cond               &
     &,    wvnew,    delv,    tnew,    tot_icenew, rimef                &
     &,    deli,     fwr,     crevp,   ventr,      delt                 &
     &,    delw,     fir,     delr,    qsinew,     qswnew               &
     &,    budget,   wsnew,   vrain2,  tot_rainnew                      &
     &,    qtnew,    qt,      wcnew,   abw                              &
     &,    aievp,    tcc,     denomf,  abi                              &
     &,    sfactor,  pidep_max,        didep,      ventis, ventil       &
     &,    dievp,    rqr,     rfactor, dwvr,       rr,     tot_rain     &
     &,    dwv0,     qsw0,    prloss,  qtrain,     vrain1               &
     &,    qsw,      ws,      esi,     esw, wv, wc, rhgrd, rho          &
     &,    rrho,     dtrho,   wsgrd,   qsi, qswgrd, qsigrd              &
     &,    tk,       tc,      pp,      bldtrh                           &
     &,    xlv,      xlv1,    xlf,     xlf1,  xlv2, denomw, denomwi     &
     &,    qwnew,    pcond,   pidep,   qrnew, qi,   qr,     qw          &
     &,    piacw,    piacwi,  piacwr,  qv,    dwvi                      &
     &,    arainnew, thick,   asnownew                                  &
     &,    qinew,    qi_min_0c, QSW_l, QSI_l, QSW0_l, SCHMIT_FAC
    
!
!
!#######################################################################
!########################## Begin Execution ############################
!#######################################################################
!
      DTPH   = DTPG / mic_step
      ARAING = 0.    ! Total Accumulated rainfall at surface (kg/m**2)
      ASNOWG = 0.    ! Total Accumulated snowfall at surface (kg/m**2)
!
      do ntimes =1,mic_step
!
        QI_min_0C = 10.E3*MASSI(MDImin)   !- Ver5
        ARAIN = 0.   ! Accumulated rainfall at surface for this step (kg/m**2)
        ASNOW = 0.   ! Accumulated snowfall at surface for this step (kg/m**2)

        INDEXR   = MDRmin
!
!-----------------------------------------------------------------------
!
        DO L=1,LSFC      !      Loop from top (L=1) to surface (L=LSFC)

!---      Skip this level and go to the next lower level if no condensate 
!         and very low specific humidities
!
          IF (QV_col(L) > EPSQ .OR. WC_col(L) > EPSQ) THEN
!
!-----------------------------------------------------------------------
!------------ Proceed with cloud microphysics calculations -------------
!-----------------------------------------------------------------------
!
            TK = T_col(L)         ! Temperature (deg K)
            TC = TK-T0C           ! Temperature (deg C)
            PP = P_col(L)         ! Pressure (Pa)
            QV = QV_col(L)        ! Specific humidity of water vapor (kg/kg)
!           WV = QV/(1.-QV)       ! Water vapor mixing ratio (kg/kg)
            WV = QV               ! Water vapor specific humidity (kg/kg)
            WC = WC_col(L)        ! Grid-scale mixing ratio of total condensate
                                ! (water or ice; kg/kg)
!           WC = WC/(1.-WC)
            RHgrd = RHC_col(L)

!
!   Pressure dependen scaling factor for flgmin (tunable)
!
!!!         pfac = max(0.5, (min(1.0, pp*0.00002))**2)   ! commented on 02182011
!go         pfac = max(0.5, (sqrt(min(1.0, pp*0.00004))))
            pfac = 1.0
!
            CLEAR = .TRUE.
!    
!--- Check grid-scale saturation when no condensate is present
!    
            ESW = min(PP, FPVSL(TK))     ! Saturation vapor pressure w/r/t water
!           QSW = EPS*ESW/(PP-ESW)       ! Saturation mixing ratio w/r/t water
            QSW = EPS*ESW/(PP+epsm1*ESW) ! Saturation specific humidity  w/r/t water
            WS  = QSW                    ! General saturation mixing ratio (water/ice)
            QSI = QSW
            IF (TC < 0.) THEN
              ESI = min(PP,FPVSI(TK))      ! Saturation vapor pressure w/r/t ice
!             QSI = EPS*ESI/(PP-ESI)       ! Saturation mixing ratio w/r/t water
              QSI = EPS*ESI/(PP+epsm1*ESI) ! Saturation specific humidity w/r/t water
              WS  = QSI                    ! General saturation mixing ratio (water/ice)
              if (pp <= esi) ws = wv / rhgrd
            ENDIF
!
            dum  = min(PP, ESW0)
            QSW0 = EPS*dum/(PP+epsm1*dum)  ! Saturation specific Humidity at 0C
!
!--- Effective grid-scale Saturation mixing ratios
!
            QSWgrd = RHgrd*QSW
            QSIgrd = RHgrd*QSI
            WSgrd  = RHgrd*WS
            QSW_l  = QSWgrd
            QSI_l  = QSIgrd
            QSW0_l = QSW0*RHgrd
!
!--- Check if air is subsaturated and w/o condensate
!
            IF (WV > WSgrd .OR. WC > EPSQ) CLEAR = .FALSE.  ! Cloudy case
            IF (ARAIN > CLIMIT) THEN ! If any rain is falling into layer from above
              CLEAR = .FALSE.
            ELSE
              ARAIN = 0.
            ENDIF
!
!--- Check if any ice is falling into layer from above
!
!--- NOTE that "SNOW" in variable names is synonomous with 
!    large, precipitation ice particles
!
            IF (ASNOW > CLIMIT) THEN
              CLEAR = .FALSE.
            ELSE
              ASNOW = 0.
            ENDIF
!
!-----------------------------------------------------------------------
!-- Loop to the end if in clear, subsaturated air free of condensate ---
!-----------------------------------------------------------------------
!
            IF (.not. CLEAR) THEN
!
!-----------------------------------------------------------------------
!--------- Initialize RHO, THICK & microphysical processes -------------
!-----------------------------------------------------------------------
!
!
!--- Virtual temperature, TV=T*(1./EPS-1)*Q, Q is specific humidity;
!    (see pp. 63-65 in Fleagle & Businger, 1963)
!
              RHO    = PP/(RD*TK*(1.+EPS1*QV)) ! Air density (kg/m**3)
              RRHO   = 1./RHO                  ! Reciprocal of air density
              DTRHO  = DTPH*RHO                ! Time step * air density
              BLDTRH = BLEND*DTRHO             ! Blend parameter * time step * air density
              THICK  = THICK_col(L)   ! Layer thickness = RHO*DZ = -DP/G
!
              ARAINnew = 0.           ! Updated accumulated rainfall at surface
              ASNOWnew = 0.           ! Updated accumulated snowfall at surface
              QI       = QI_col(L)    ! Ice mixing ratio
              QInew    = 0.           ! Updated ice mixing ratio
              QR       = QR_col(L)    ! Rain mixing ratio
              QRnew    = 0.           ! Updated rain ratio
              QW       = QW_col(L)    ! Cloud water mixing ratio
              QWnew    = 0.           ! Updated cloud water ratio
!
              PCOND    = 0.       ! Condensation (>0) or evaporation (<0)
                                  ! of cloud water (kg/kg)
              PIDEP    = 0.       ! Deposition (>0) or sublimation (<0)
                                  ! of ice crystals (kg/kg)
              PIACW   = 0.        ! Cloud water collection (riming)
                                  ! by precipitation ice (kg/kg; >0)
              PIACWI  = 0.        ! Growth of precip ice by riming (kg/kg; >0)
              PIACWR  = 0.        ! Shedding of accreted cloud water
                                  ! to form rain (kg/kg; >0)
              PIACR   = 0.        ! Freezing of rain onto large ice
                                  ! at supercooled temps (kg/kg; >0)
              PICND   = 0.        ! Condensation (>0) onto wet, melting
                                  ! ice (kg/kg)
              PIEVP   = 0.        ! Evaporation (<0) from wet, melting
                                  ! ice (kg/kg)
              PIMLT   = 0.        ! Melting ice (kg/kg; >0)
              PRAUT   = 0.        ! Cloud water autoconversion to rain (kg/kg; >0)
              PRACW   = 0.        ! Cloud water collection (accretion) by rain (kg/kg; >0)
              PREVP   = 0.        ! Rain evaporation (kg/kg; <0)
!
!---------------------------------------------------------------------------
!--- Double check input hydrometeor mixing ratios
!
!             DUM  = WC - (QI+QW+QR)
!             DUM1 = ABS(DUM)
!             DUM2 = TOLER * MIN(WC, QI+QW+QR)
!             IF (DUM1 >  DUM2) THEN
!               WRITE(6,"(/2(a,i4),a,i2)") '{@ i=',I_index,' j=',J_index,
!     &                                     ' L=',L
!               WRITE(6,"(4(a12,g11.4,1x))") 
!     & '{@ TCold=',TC,'P=',.01*PP,'DIFF=',DUM,'WCold=',WC,
!     & '{@ QIold=',QI,'QWold=',QW,'QRold=',QR
!             ENDIF
!
!***********************************************************************
!*********** MAIN MICROPHYSICS CALCULATIONS NOW FOLLOW! ****************
!***********************************************************************
!
!--- Calculate a few variables, which are used more than once below
!
!--- Latent heat of vaporization as a function of temperature from
!      Bolton (1980, JAS)
!
              XLV    = 3.148E6 - 2370*TK     ! Latent heat of vaporization (Lv)
              XLF    = XLS-XLV               ! Latent heat of fusion (Lf)
              XLV1   = XLV*RCP               ! Lv/Cp
              XLF1   = XLF*RCP               ! Lf/Cp
              TK2    = 1./(TK*TK)            ! 1./TK**2
              XLV2   = XLV*XLV*QSW_l*TK2/RV  ! Lv**2*Qsw_l/(Rv*TK**2)
              DENOMW = 1. + XLV2*RCP         ! Denominator term, Clausius-Clapeyron correction
!
!--- Basic thermodynamic quantities
!      *      DYNVIS     - dynamic viscosity           [ kg/(m*s) ]
!      *      THERM_COND - thermal conductivity        [ J/(m*s*K) ]
!      *      DIFFUS     - diffusivity of water vapor  [ m**2/s ]
!
!             TFACTOR    = TK**1.5/(TK+120.)
              TFACTOR    = TK*sqrt(TK)/(TK+120.)
              DYNVIS     = 1.496E-6*TFACTOR
              THERM_COND = 2.116E-3*TFACTOR
              DIFFUS     = 8.794E-5*TK**1.81/PP
              SCHMIT_FAC = (RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
!
!--- Air resistance term for the fall speed of ice following the
!      basic research by Heymsfield, Kajikawa, others 
!
              GAMMAS = (1.E5/PP)**C1
!
!--- Air resistance for rain fall speed (Beard, 1985, JAOT, p.470)
!
              GAMMAR = (RHO0/RHO)**0.4
!
!----------------------------------------------------------------------
!-------------  IMPORTANT MICROPHYSICS DECISION TREE  -----------------
!----------------------------------------------------------------------
!
!--- Determine if conditions supporting ice are present
!
              IF (TC < 0. .OR. QI > EPSQ .OR. ASNOW > CLIMIT) THEN
                ICE_logical = .TRUE.
              ELSE
                ICE_logical = .FALSE.
                QLICE = 0.
                QTICE = 0.
              ENDIF
!
!--- Determine if rain is present
!
              RAIN_logical = .FALSE.
              IF (ARAIN > CLIMIT .OR. QR > EPSQ) RAIN_logical = .TRUE.
!
              IF (ICE_logical) THEN
!
!--- IMPORTANT:  Estimate time-averaged properties.
!
!---
!  * FLARGE  - ratio of number of large ice to total (large & small) ice
!  * FSMALL  - ratio of number of small ice crystals to large ice particles
!  ->  Small ice particles are assumed to have a mean diameter of 50 microns.
!  * XSIMASS - used for calculating small ice mixing ratio
!---
!  * TOT_ICE - total mass (small & large) ice before microphysics,
!              which is the sum of the total mass of large ice in the 
!              current layer and the input flux of ice from above
!  * PILOSS  - greatest loss (<0) of total (small & large) ice by
!              sublimation, removing all of the ice falling from above
!              and the ice within the layer
!  * RimeF1  - Rime Factor, which is the mass ratio of total (unrimed & rimed) 
!              ice mass to the unrimed ice mass (>=1)
!  * VrimeF  - the velocity increase due to rime factor or melting (ratio, >=1)
!  * VSNOW   - Fall speed of rimed snow w/ air resistance correction
!  * EMAIRI  - equivalent mass of air associated layer and with fall of snow into layer
!  * XLIMASS - used for calculating large ice mixing ratio
!  * FLIMASS - mass fraction of large ice
!  * QTICE   - time-averaged mixing ratio of total ice
!  * QLICE   - time-averaged mixing ratio of large ice
!  * NLICE   - time-averaged number concentration of large ice
!  * NSmICE  - number concentration of small ice crystals at current level
!---
!--- Assumed number fraction of large ice particles to total (large & small) 
!    ice particles, which is based on a general impression of the literature.
!
                WVQW = WV + QW                ! Water vapor + cloud water
!
!--- 6/19/03 - Deleted some code here ....
!
!  *********************************************************

!               IF (TC >= 0. .OR. WVQW < QSIgrd) THEN
!  !
!  !--- Eliminate small ice particle contributions for melting & sublimation
!  !
!                 FLARGE = FLARGE1
!               ELSE
!  !
!  !--- Enhanced number of small ice particles during depositional growth
!  !    (effective only when 0C > T >= T_ice [-10C] )
!  !
!                 FLARGE = FLARGE2
!  !
!  !--- Larger number of small ice particles due to rime splintering
!  !
!                 IF (TC >= -8. .AND. TC <= -3.) FLARGE=.5*FLARGE
!
!               ENDIF            ! End IF (TC >= 0. .OR. WVQW < QSIgrd)
!               FSMALL=(1.-FLARGE)/FLARGE
!               XSIMASS=RRHO*MASSI(MDImin)*FSMALL
!  *********************************************************
!
                IF (QI <= EPSQ .AND. ASNOW <= CLIMIT) THEN
                  INDEXS  = MDImin
                  FLARGE  = 0.                   !--- Begin 6/19/03 changes
                  FSMALL  = 1.
                  XSIMASS = RRHO*MASSI(MDImin)   !--- End 6/19/03 changes
                  TOT_ICE = 0.
                  PILOSS  = 0.
                  RimeF1  = 1.
                  VrimeF  = 1.
                  VEL_INC = GAMMAS
                  VSNOW   = 0.
                  EMAIRI  = THICK
                  XLIMASS = RRHO*RimeF1*MASSI(INDEXS)
                  FLIMASS = XLIMASS/(XLIMASS+XSIMASS)
                  QLICE   = 0.
                  QTICE   = 0.
                  NLICE   = 0.
                  NSmICE  = 0.
                ELSE
   !
   !--- For T<0C mean particle size follows Houze et al. (JAS, 1979, p. 160), 
   !    converted from Fig. 5 plot of LAMDAs.  Similar set of relationships 
   !    also shown in Fig. 8 of Ryan (BAMS, 1996, p. 66).
   !
   !--- Begin 6/19/03 changes => allow NLImax to increase & FLARGE to 
   !    decrease at COLDER temperatures; set FLARGE to zero (i.e., only small
   !    ice) if the ice mixing ratio is below QI_min

!                 DUM    = MAX(0.05, MIN(1., EXP(.0536*TC)) )
                  DUM    = MAX(0.05, MIN(1., EXP(.0564*TC)) )
                  INDEXS = MIN(MDImax, MAX(MDImin, INT(XMImax*DUM) ) )
!
                  DUM    = MAX(FLGmin*pfac, DUM)

                  QI_min = QI_min_0C * dum  !- Ver5    ----Not used ----
!!                QI_min = QI_min_0C        !- Ver5
!!!               QI_min = QI_min_0C/DUM    !- Ver5

                  NLImax = 10.E3/sqrt(DUM)  !- Ver3
                  IF (TC < 0.) THEN
                    FLARGE = DUM            !- Ver4
                  ELSE
                    FLARGE = 1.
                  ENDIF
                  FSMALL  = (1.-FLARGE)/FLARGE
                  XSIMASS = RRHO*MASSI(MDImin)*FSMALL
                  TOT_ICE = THICK*QI + BLEND*ASNOW
                  PILOSS  = -TOT_ICE/THICK
                  LBEF    = MAX(1,L-1)
                  RimeF1  = (RimeF_col(L)*THICK*QI                      &
     &                    + RimeF_col(LBEF)*BLEND*ASNOW)/TOT_ICE
                  RimeF1  = MIN(RimeF1, RFmax)
                  VSNOW   = 0.0
                  DO IPASS=0,1
                    IF (RimeF1 .LE. 1.) THEN
                      RimeF1 = 1.
                      VrimeF = 1.
                    ELSE
                      IXS  = MAX(2, MIN(INDEXS/100, 9))
                      XRF  = 10.492*LOG(RimeF1)
                      IXRF = MAX(0, MIN(INT(XRF), Nrime))
                      IF (IXRF .GE. Nrime) THEN
                        VrimeF = VEL_RF(IXS,Nrime)
                      ELSE
                        VrimeF = VEL_RF(IXS,IXRF)+(XRF-FLOAT(IXRF))*    &
     &                          (VEL_RF(IXS,IXRF+1)-VEL_RF(IXS,IXRF))
                      ENDIF
                    ENDIF            ! End IF (RimeF1 <= 1.)
                    VEL_INC = GAMMAS*VrimeF
                    VSNOW   = VEL_INC*VSNOWI(INDEXS)
                    EMAIRI  = THICK + BLDTRH*VSNOW
                    XLIMASS = RRHO*RimeF1*MASSI(INDEXS)
                    FLIMASS = XLIMASS/(XLIMASS+XSIMASS)
                    QTICE   = TOT_ICE/EMAIRI
                    QLICE   = FLIMASS*QTICE
                    NLICE   = QLICE/XLIMASS
                    NSmICE  = Fsmall*NLICE
   !
                    IF ( (NLICE >= NLImin .AND. NLICE <= NLImax)        & 
     &                    .OR. IPASS == 1) THEN
                      EXIT
                    ELSE
                      IF(TC < 0) THEN
                        XLI = RHO*(QTICE/DUM-XSIMASS)/RimeF1
                        IF (XLI <= MASSI(MDImin) ) THEN
                          INDEXS = MDImin
                        ELSE IF (XLI <= MASSI(450) ) THEN
                          DLI    = 9.5885E5*XLI**.42066       ! DLI in microns
                          INDEXS = MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                        ELSE IF (XLI <= MASSI(MDImax) ) THEN
                          DLI    = 3.9751E6*XLI**.49870       ! DLI in microns
                          INDEXS = MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                        ELSE 
                          INDEXS = MDImax
                        ENDIF             ! End IF (XLI <= MASSI(MDImin) ) 
                      ENDIF               ! End IF (TC < 0)
!
!--- Reduce excessive accumulation of ice at upper levels
!    associated with strong grid-resolved ascent
!
!--- Force NLICE to be between NLImin and NLImax
!
!--- 8/22/01: Increase density of large ice if maximum limits
!    are reached for number concentration (NLImax) and mean size
!    (MDImax).  Done to increase fall out of ice.
!
!

                      DUM = MAX(NLImin, MIN(NLImax, NLICE) )
                      IF (DUM >= NLImax .AND. INDEXS >= MDImax)         &
     &                 RimeF1 = RHO*(QTICE/NLImax-XSIMASS)/MASSI(INDEXS)
!
!                WRITE(6,"(4(a12,g11.4,1x))") 
!     & '{$ TC=',TC,'P=',.01*PP,'NLICE=',NLICE,'DUM=',DUM,
!     & '{$ XLI=',XLI,'INDEXS=',FLOAT(INDEXS),'RHO=',RHO,'QTICE=',QTICE,
!     & '{$ XSIMASS=',XSIMASS,'RimeF1=',RimeF1

                    ENDIF    ! End IF ( (NLICE >=NLImin .AND. NLICE >= NLImax)
                  ENDDO      ! End DO IPASS=0,1
                ENDIF        ! End IF (QI <= EPSQ .AND. ASNOW <= CLIMIT)
              ENDIF          ! End IF (ICE_logical)
!
!----------------------------------------------------------------------
!--------------- Calculate individual processes -----------------------
!----------------------------------------------------------------------
!
!--- Cloud water autoconversion to rain and collection by rain
!
              IF (QW > EPSQ .AND. TC >= T_ICE) THEN
   !
   !--- QW0 could be modified based on land/sea properties, 
   !      presence of convection, etc.  This is why QAUT0 and CRAUT
   !      are passed into the subroutine as externally determined
   !      parameters.  Can be changed in the future if desired.
   !
!               QW0   = QAUT0*RRHO
                QW0   = QAUTx*RRHO*XNCW(L)
                PRAUT = MAX(0., QW-QW0)*CRAUT
                IF (QLICE  >  EPSQ) THEN
      !
      !--- Collection of cloud water by large ice particles ("snow")
      !    PIACWI=PIACW for riming, PIACWI=0 for shedding
      !
!Moor              FWS   = MIN(1., CIACW*VEL_INC*NLICE*ACCRI(INDEXS)/PP**C1) ! 20050422
                   FWS   = MIN(0.1, CIACW*VEL_INC*NLICE*ACCRI(INDEXS)   &
     &                                                  /PP**C1)
                   PIACW = FWS*QW
                   IF (TC  < 0.) PIACWI = PIACW    ! Large ice riming

                ENDIF             ! End IF (QLICE > EPSQ)
              ENDIF               ! End IF (QW > EPSQ .AND. TC >= T_ICE)
!
!----------------------------------------------------------------------
!--- Loop around some of the ice-phase processes if no ice should be present
!----------------------------------------------------------------------
!
              IF (ICE_logical) THEN
!
!--- Now the pretzel logic of calculating ice deposition
!
                IF (TC < T_ICE .AND. (WV > QSIgrd .OR. QW > EPSQ)) THEN
!
!--- Adjust to ice saturation at T<T_ICE (-10C) if supersaturated.
!    Sources of ice due to nucleation and convective detrainment are
!    either poorly understood, poorly resolved at typical NWP 
!    resolutions, or are not represented (e.g., no detrained 
!    condensate in BMJ Cu scheme).
!    
                  PCOND = -QW
                  DUM1  = TK + XLV1*PCOND              ! Updated (dummy) temperature (deg K)
                  DUM2  = WV+QW                        ! Updated (dummy) water vapor mixing ratio
                  DUM   = min(pp,FPVSI(DUM1))          ! Updated (dummy) saturation vapor pressure w/r/t ice
                  DUM   = RHgrd*EPS*DUM/(pp+epsm1*dum) ! Updated (dummy) saturation specific humidity w/r/t ice
!                 DUM   = RHgrd*EPS*DUM/(PP-DUM)       ! Updated (dummy) saturation mixing ratio w/r/t ice

                  IF (DUM2 > DUM) PIDEP = DEPOSIT(PP, RHgrd, DUM1, DUM2)

                  DWVi = 0.                            ! Used only for debugging
!
                ELSE IF (TC < 0.) THEN
!
!--- These quantities are handy for ice deposition/sublimation
!    PIDEP_max - max deposition or minimum sublimation to ice saturation
!
                  DENOMI    = 1. + XLS2*QSI_l*TK2
                  DWVi      = MIN(WVQW,QSW_l)-QSI_l
                  PIDEP_max = MAX(PILOSS, DWVi/DENOMI)
                  IF (QTICE > 0.) THEN
!
!--- Calculate ice deposition/sublimation
!      * SFACTOR - [VEL_INC**.5]*[Schmidt**(1./3.)]*[(RHO/DYNVIS)**.5],
!        where Schmidt (Schmidt Number) =DYNVIS/(RHO*DIFFUS)
!      * Units: SFACTOR - s**.5/m ;  ABI - m**2/s ;  NLICE - m**-3 ;
!               VENTIL, VENTIS - m**-2 ;  VENTI1 - m ;  
!               VENTI2 - m**2/s**.5 ; DIDEP - unitless
!
!                   SFACTOR = VEL_INC**.5*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
                    SFACTOR = sqrt(VEL_INC)*SCHMIT_FAC
                    ABI     = 1./(RHO*XLS3*QSI*TK2/THERM_COND+1./DIFFUS)
!
!--- VENTIL - Number concentration * ventilation factors for large ice
!--- VENTIS - Number concentration * ventilation factors for small ice
!
!--- Variation in the number concentration of ice with time is not
!      accounted for in these calculations (could be in the future).
!
                    VENTIL = (VENTI1(INDEXS) + SFACTOR*VENTI2(INDEXS))  &
     &                                       * NLICE
                    VENTIS = (VENTI1(MDImin) + SFACTOR*VENTI2(MDImin))  &
     &                                       * NSmICE
                    DIDEP  = ABI*(VENTIL+VENTIS)*DTPH
!
!--- Account for change in water vapor supply w/ time
!
                    IF (DIDEP >= Xratio)                                & 
     &                DIDEP = (1.-EXP(-DIDEP*DENOMI))/DENOMI
                    IF (DWVi > 0.) THEN
                      PIDEP = MIN(DWVi*DIDEP, PIDEP_max)
                    ELSE IF (DWVi < 0.) THEN
                      PIDEP = MAX(DWVi*DIDEP, PIDEP_max)
                    ENDIF
!
                  ELSE IF (WVQW > QSI_l .AND. TC <= T_ICE_init) THEN
!
!--- Ice nucleation in near water-saturated conditions.  Ice crystal
!    growth during time step calculated using Miller & Young (1979, JAS).
!--- These deposition rates could drive conditions below water saturation,
!    which is the basis of these calculations.  Intended to approximate
!    more complex & computationally intensive calculations.
!
                    INDEX_MY = MAX(MY_T1, MIN( INT(.5-TC), MY_T2 ) )
!
!--- DUM1 is the supersaturation w/r/t ice at water-saturated conditions
!
!--- DUM2 is the number of ice crystals nucleated at water-saturated 
!    conditions based on Meyers et al. (JAM, 1992).
!
!--- Prevent unrealistically large ice initiation (limited by PIDEP_max)
!      if DUM2 values are increased in future experiments
!
                    DUM1  = QSW/QSI - 1.      
                    DUM2  = 1.E3*EXP(12.96*DUM1-0.639)
                    PIDEP = MIN(PIDEP_max,DUM2*MY_GROWTH(INDEX_MY)*RRHO)
!
                  ENDIF ! End IF (QTICE > 0.)
!
                ENDIF   ! End IF (TC < T_ICE .AND. (WV > QSIgrd .OR. QW > EPSQ))
!
!------------------------------------------------------------------------
!
              ENDIF     ! End of IF (ICE_logical)then loop
!
!------------------------------------------------------------------------
!
!--- Cloud water condensation
!
              IF (TC >= T_ICE .AND. (QW > EPSQ .OR. WV > QSWgrd)) THEN
                IF (PIACWI == 0. .AND. PIDEP == 0.) THEN
                  PCOND = CONDENSE (PP, QW, RHgrd, TK, WV)
                ELSE  !-- Modify cloud condensation in response to ice processes
                  DUM     = XLV*QSWgrd*RCPRV*TK2
                  DENOMWI = 1. + XLS*DUM
                  DENOMF  = XLF*DUM
                  DUM     = MAX(0., PIDEP)
                  PCOND   = (WV-QSWgrd-DENOMWI*DUM-DENOMF*PIACWI)/DENOMW
                  DUM1    = -QW
                  DUM2    = PCOND - PIACW
                  IF (DUM2  <  DUM1) THEN    !--- Limit cloud water sinks
                    DUM    = DUM1/DUM2
                    PCOND  = DUM*PCOND
                    PIACW  = DUM*PIACW
                    PIACWI = DUM*PIACWI
                  ENDIF ! EndIF (DUM2 <  DUM1)
                ENDIF   ! EndIF (PIACWI == 0. .AND. PIDEP == 0.)
              ENDIF     ! EndIF (TC >= T_ICE .AND. (QW > EPSQ .OR. WV > QSWgrd))
!
!--- Limit freezing of accreted rime to prevent temperature oscillations,
!    a crude Schumann-Ludlam limit (p. 209 of Young, 1993). 
!
              TCC = TC + XLV1*PCOND + XLS1*PIDEP + XLF1*PIACWI
              IF (TCC  >  0.) THEN
                PIACWI = 0.
                TCC    = TC + XLV1*PCOND + XLS1*PIDEP
              ENDIF
!
              IF (TC > 0. .AND. TCC > 0. .AND. ICE_logical) THEN
!
!--- Calculate melting and evaporation/condensation
!      * Units: SFACTOR - s**.5/m ;  ABI - m**2/s ;  NLICE - m**-3 ;
!               VENTIL - m**-2 ;  VENTI1 - m ;  
!               VENTI2 - m**2/s**.5 ; CIEVP - /s
!
!               SFACTOR = VEL_INC**.5*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
                SFACTOR = sqrt(VEL_INC)*SCHMIT_FAC
                VENTIL  = NLICE*(VENTI1(INDEXS)+SFACTOR*VENTI2(INDEXS))
                AIEVP   = VENTIL*DIFFUS*DTPH
                IF (AIEVP  <  Xratio) THEN
                  DIEVP = AIEVP
                ELSE
                  DIEVP = 1. - EXP(-AIEVP)
                ENDIF
!               QSW0 = EPS*ESW0/(PP-ESW0)
!               QSW0 = EPS*ESW0/(PP+epsm1*ESW0)
!!              dum  = min(PP, ESW0)
!!              QSW0 = EPS*dum/(PP+epsm1*dum)
!               DWV0 = MIN(WV,QSW)-QSW0
                DWV0 = MIN(WV,QSW_l)-QSW0_l
                DUM  = QW + PCOND
                IF (WV < QSW_l .AND. DUM <= EPSQ) THEN
   !
   !--- Evaporation from melting snow (sink of snow) or shedding
   !    of water condensed onto melting snow (source of rain)
   !
                  DUM   = DWV0*DIEVP
                  PIEVP = MAX( MIN(0., DUM), PILOSS)
                  PICND = MAX(0., DUM)
                ENDIF            ! End IF (WV < QSW_l .AND. DUM <= EPSQ)
                PIMLT = THERM_COND*TCC*VENTIL*RRHO*DTPH/XLF
   !
   !--- Limit melting to prevent temperature oscillations across 0C
   !
                DUM1  = MAX( 0., (TCC+XLV1*PIEVP)/XLF1 )
                PIMLT = MIN(PIMLT, DUM1)
   !
   !--- Limit loss of snow by melting (>0) and evaporation
   !
                DUM = PIEVP - PIMLT
                IF (DUM < PILOSS) THEN
                  DUM1  = PILOSS/DUM
                  PIMLT = PIMLT*DUM1
                  PIEVP = PIEVP*DUM1
                ENDIF       ! End IF (DUM  > QTICE)
              ENDIF         ! End IF (TC > 0. .AND. TCC > 0. .AND. ICE_logical)
!
!--- IMPORTANT:  Estimate time-averaged properties.
!
!  * TOT_RAIN - total mass of rain before microphysics, which is the sum of
!               the total mass of rain in the current layer and the input 
!               flux of rain from above
!  * VRAIN1   - fall speed of rain into grid from above (with air resistance correction)
!  * QTRAIN   - time-averaged mixing ratio of rain (kg/kg)
!  * PRLOSS   - greatest loss (<0) of rain, removing all rain falling from
!               above and the rain within the layer
!  * RQR      - rain content (kg/m**3)
!  * INDEXR   - mean size of rain drops to the nearest 1 micron in size
!  * N0r      - intercept of rain size distribution (typically 10**6 m**-4)
!
              TOT_RAIN = 0.
              VRAIN1   = 0.
              QTRAIN   = 0.
              PRLOSS   = 0.
              RQR      = 0.
              N0r      = 0.
              INDEXR   = MDRmin
              INDEXR1  = INDEXR    ! For debugging only
              IF (RAIN_logical) THEN
                IF (ARAIN <= 0.) THEN
                  INDEXR = MDRmin
                  VRAIN1 = 0.
                ELSE
   !
   !--- INDEXR (related to mean diameter) & N0r could be modified 
   !      by land/sea properties, presence of convection, etc.
   !
   !--- Rain rate normalized to a density of 1.194 kg/m**3
   !
                  RR = ARAIN / (DTPH*GAMMAR)
   !
                  IF (RR <= RR_DRmin) THEN
        !
        !--- Assume fixed mean diameter of rain (0.2 mm) for low rain rates, 
        !      instead vary N0r with rain rate
        !
                    INDEXR = MDRmin
                  ELSE IF (RR <= RR_DR1) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.05 and 0.10 mm:
        !      V(Dr)=5.6023e4*Dr**1.136, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*5.6023e4*Dr**(4+1.136) = 1.408e15*Dr**5.136,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.123e-3*RR**.1947 -> Dr (microns) = 1.123e3*RR**.1947
        !
                    INDEXR = INT( 1.123E3*RR**.1947 + .5 )
                    INDEXR = MAX( MDRmin, MIN(INDEXR, MDR1) )

                  ELSE IF (RR <= RR_DR2) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.10 and 0.20 mm:
        !      V(Dr)=1.0867e4*Dr**.958, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*1.0867e4*Dr**(4+.958) = 2.731e14*Dr**4.958,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.225e-3*RR**.2017 -> Dr (microns) = 1.225e3*RR**.2017
        !
                    INDEXR = INT( 1.225E3*RR**.2017 + .5 )
                    INDEXR = MAX( MDR1, MIN(INDEXR, MDR2) )

                  ELSE IF (RR <= RR_DR3) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.20 and 0.32 mm:
        !      V(Dr)=2831.*Dr**.80, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*2831.*Dr**(4+.80) = 7.115e13*Dr**4.80, 
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.3006e-3*RR**.2083 -> Dr (microns) = 1.3006e3*RR**.2083
        !
                    INDEXR = INT( 1.3006E3*RR**.2083 + .5 )
                    INDEXR = MAX( MDR2, MIN(INDEXR, MDR3) )

                  ELSE IF (RR <= RR_DRmax) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.32 and 0.45 mm:
        !      V(Dr)=944.8*Dr**.6636, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*944.8*Dr**(4+.6636) = 2.3745e13*Dr**4.6636,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.355e-3*RR**.2144 -> Dr (microns) = 1.355e3*RR**.2144
        !
                    INDEXR = INT( 1.355E3*RR**.2144 + .5 )
                    INDEXR = MAX( MDR3, MIN(INDEXR, MDRmax) )
                  ELSE 
        !
        !--- Assume fixed mean diameter of rain (0.45 mm) for high rain rates, 
        !      instead vary N0r with rain rate
        !
                    INDEXR = MDRmax
                  ENDIF               ! End IF (RR <= RR_DRmin) etc. 
!
                  VRAIN1 = GAMMAR*VRAIN(INDEXR)
                ENDIF                 ! End IF (ARAIN <= 0.)
!
                INDEXR1  = INDEXR     ! For debugging only
                TOT_RAIN = THICK*QR+BLEND*ARAIN
                QTRAIN   = TOT_RAIN/(THICK+BLDTRH*VRAIN1)
                PRLOSS   = -TOT_RAIN/THICK
                RQR      = RHO*QTRAIN
   !
   !--- RQR - time-averaged rain content (kg/m**3)
   !
                IF (RQR <= RQR_DRmin) THEN
                  N0r    = MAX(N0rmin, CN0r_DMRmin*RQR)
                  INDEXR = MDRmin
                ELSE IF (RQR >= RQR_DRmax) THEN
                  N0r    = CN0r_DMRmax*RQR
                  INDEXR = MDRmax
                ELSE
                  N0r    = N0r0
!                 INDEXR = MAX( XMRmin, MIN(CN0r0*RQR**.25, XMRmax) )
                  item   = CN0r0*sqrt(sqrt(RQR))               ! Moorthi 07/31/08
                  INDEXR = MAX( MDRmin, MIN(item, MDRmax) )    ! Moorthi 07/31/08
                ENDIF
   !
                IF (TC < T_ICE) THEN
                  PIACR = -PRLOSS
                ELSE
                  DWVr = WV - PCOND - QSW_l
                  DUM  = QW + PCOND
                  IF (DWVr < 0. .AND. DUM <= EPSQ) THEN
!
!--- Rain evaporation
!
!    * RFACTOR - [GAMMAR**.5]*[Schmidt**(1./3.)]*[(RHO/DYNVIS)**.5],
!        where Schmidt (Schmidt Number) =DYNVIS/(RHO*DIFFUS)
!
!    * Units: RFACTOR - s**.5/m ;  ABW - m**2/s ;  VENTR - m**-2 ;  
!             N0r - m**-4 ;  VENTR1 - m**2 ;  VENTR2 - m**3/s**.5 ;
!             CREVP - unitless
!
!                   RFACTOR = GAMMAR**.5*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
                    RFACTOR = sqrt(GAMMAR)*SCHMIT_FAC
                    ABW     = 1./(RHO*XLV2/THERM_COND+1./DIFFUS)
!
!--- Note that VENTR1, VENTR2 lookup tables do not include the 
!      1/Davg multiplier as in the ice tables
!
                    VENTR = N0r*(VENTR1(INDEXR)+RFACTOR*VENTR2(INDEXR))
                    CREVP = ABW*VENTR*DTPH
                    IF (CREVP < Xratio) THEN
                      DUM = DWVr*CREVP
                    ELSE
                      DUM = DWVr*(1.-EXP(-CREVP*DENOMW))/DENOMW
                    ENDIF
                    PREVP = MAX(DUM, PRLOSS)
                  ELSE IF (QW > EPSQ) THEN
                    FWR   = CRACW*GAMMAR*N0r*ACCRR(INDEXR)
!Moor               PRACW = MIN(1.,FWR)*QW               ! 20050422
                    PRACW = MIN(0.1,FWR)*QW
                  ENDIF           ! End IF (DWVr < 0. .AND. DUM <= EPSQ)
!
                  IF (TC < 0. .AND. TCC < 0.) THEN
!
!--- Biggs (1953) heteorogeneous freezing (e.g., Lin et al., 1983)
!   - Rescaled mean drop diameter from microns (INDEXR) to mm (DUM) to prevent underflow
!
                    DUM   = .001*FLOAT(INDEXR)
                    dum1  = dum * dum
                    DUM   = (EXP(ABFR*TC)-1.)*DUM1*DUM1*DUM1*DUM
                    PIACR = MIN(CBFR*N0r*RRHO*DUM, QTRAIN)
                    IF (QLICE > EPSQ) THEN
            !
            !--- Freezing of rain by collisions w/ large ice
            !
                      DUM  = GAMMAR*VRAIN(INDEXR)
                      DUM1 = DUM-VSNOW
            !
            !--- DUM2 - Difference in spectral fall speeds of rain and
            !      large ice, parameterized following eq. (48) on p. 112 of 
            !      Murakami (J. Meteor. Soc. Japan, 1990)
            !
                      DUM2 = (DUM1*DUM1+.04*DUM*VSNOW)**.5
                      DUM1 = 5.E-12*INDEXR*INDEXR+2.E-12*INDEXR*INDEXS      &
     &                     +.5E-12*INDEXS*INDEXS
                      FIR = MIN(1., CIACR*NLICE*DUM1*DUM2)
            !
            !--- Future?  Should COLLECTION BY SMALL ICE SHOULD BE INCLUDED???
            !
                      PIACR = MIN(PIACR+FIR*QTRAIN, QTRAIN)
                    ENDIF        ! End IF (QLICE >  EPSQ)
                    DUM = PREVP - PIACR
                    If (DUM < PRLOSS) THEN
                      DUM1  = PRLOSS/DUM
                      PREVP = DUM1*PREVP
                      PIACR = DUM1*PIACR
                    ENDIF        ! End If (DUM < PRLOSS)
                  ENDIF          ! End IF (TC < 0. .AND. TCC < 0.)
                ENDIF            ! End IF (TC < T_ICE)
              ENDIF              ! End IF (RAIN_logical) 
!
!----------------------------------------------------------------------
!---------------------- Main Budget Equations -------------------------
!----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!--- Update fields, determine characteristics for next lower layer ----
!-----------------------------------------------------------------------
!
!--- Carefully limit sinks of cloud water
!
              DUM1 = PIACW + PRAUT + PRACW - MIN(0.,PCOND)
              IF (DUM1 > QW) THEN
                DUM    = QW/DUM1
                PIACW  = DUM*PIACW
                PIACWI = DUM*PIACWI
                PRAUT  = DUM*PRAUT
                PRACW  = DUM*PRACW
                IF (PCOND < 0.) PCOND=DUM*PCOND
              ENDIF
              PIACWR = PIACW - PIACWI          ! TC >= 0C
!
!--- QWnew - updated cloud water mixing ratio
!
              DELW  = PCOND - PIACW - PRAUT - PRACW
              QWnew = QW+DELW
              IF (QWnew <=  EPSQ) QWnew = 0.
              IF (QW > 0. .AND. QWnew /= 0.) THEN
                DUM = QWnew/QW
                IF (DUM  < TOLER) QWnew = 0.
              ENDIF
!
!--- Update temperature and water vapor mixing ratios
!
              DELT = XLV1 * (PCOND+PIEVP+PICND+PREVP)                   &
     &             + XLS1 * PIDEP + XLF1*(PIACWI+PIACR-PIMLT)
              Tnew = TK + DELT
!
              DELV  = -PCOND - PIDEP - PIEVP - PICND - PREVP
              WVnew = WV + DELV
!
!--- Update ice mixing ratios
!
!---
!  * TOT_ICEnew - total mass (small & large) ice after microphysics,
!                 which is the sum of the total mass of large ice in the 
!                 current layer and the flux of ice out of the grid box below
!  * RimeF      - Rime Factor, which is the mass ratio of total (unrimed & 
!                 rimed) ice mass to the unrimed ice mass (>=1)
!  * QInew      - updated mixing ratio of total (large & small) ice in layer
!      -> TOT_ICEnew=QInew*THICK+BLDTRH*QLICEnew*VSNOW
!        -> But QLICEnew=QInew*FLIMASS, so
!      -> TOT_ICEnew=QInew*(THICK+BLDTRH*FLIMASS*VSNOW)
!  * ASNOWnew   - updated accumulation of snow at bottom of grid cell
!---
!
              DELI  = 0.
              RimeF = 1.
              IF (ICE_logical) THEN
                DELI       = PIDEP + PIEVP + PIACWI + PIACR - PIMLT
                TOT_ICEnew = TOT_ICE + THICK*DELI
                IF (TOT_ICE > 0. .AND. TOT_ICEnew /= 0.) THEN
                  DUM = TOT_ICEnew/TOT_ICE
                  IF (DUM  < TOLER) TOT_ICEnew = 0.
                ENDIF
                IF (TOT_ICEnew <= CLIMIT) THEN
                  TOT_ICEnew = 0.
                  RimeF      = 1.
                  QInew      = 0.
                  ASNOWnew   = 0.
                ELSE
      !
      !--- Update rime factor if appropriate
      !
                  DUM = PIACWI + PIACR
                  IF (DUM <= EPSQ .AND. PIDEP <= EPSQ) THEN
                    RimeF = RimeF1
                  ELSE
         !
         !--- Rime Factor, RimeF = (Total ice mass)/(Total unrimed ice mass)
         !      DUM1 - Total ice mass, rimed & unrimed
         !      DUM2 - Estimated mass of *unrimed* ice
         !
                    DUM1 = TOT_ICE + THICK*(PIDEP+DUM)
                    DUM2 = TOT_ICE/RimeF1 + THICK*PIDEP
                    IF (DUM2 <= 0.) THEN
                      RimeF = RFmax
                    ELSE
                      RimeF = MIN(RFmax, MAX(1., DUM1/DUM2) )
                    ENDIF
                  ENDIF       ! End IF (DUM <= EPSQ .AND. PIDEP <= EPSQ)
                  QInew = TOT_ICEnew/(THICK+BLDTRH*FLIMASS*VSNOW)
                  IF (QInew  <=  EPSQ) QInew = 0.
                  IF (QI > 0. .AND. QInew /= 0.) THEN
                    DUM = QInew/QI
                    IF (DUM < TOLER) QInew = 0.
                  ENDIF
                  ASNOWnew = BLDTRH*FLIMASS*VSNOW*QInew
                  IF (ASNOW > 0. .AND. ASNOWnew /= 0.) THEN
                    DUM = ASNOWnew/ASNOW
                    IF (DUM < TOLER) ASNOWnew = 0.
                  ENDIF
                ENDIF         ! End IF (TOT_ICEnew <= CLIMIT)
              ENDIF           ! End IF (ICE_logical)
!
!--- Update rain mixing ratios
!
!---
! * TOT_RAINnew - total mass of rain after microphysics
!                 current layer and the input flux of ice from above
! * VRAIN2      - time-averaged fall speed of rain in grid and below 
!                 (with air resistance correction)
! * QRnew       - updated rain mixing ratio in layer
!      -> TOT_RAINnew=QRnew*(THICK+BLDTRH*VRAIN2)
!  * ARAINnew  - updated accumulation of rain at bottom of grid cell
!---
!
              DELR        = PRAUT + PRACW + PIACWR - PIACR + PIMLT      &
     &                    + PREVP + PICND
              TOT_RAINnew = TOT_RAIN+THICK*DELR
              IF (TOT_RAIN > 0. .AND. TOT_RAINnew /= 0.) THEN
                DUM = TOT_RAINnew/TOT_RAIN
                IF (DUM < TOLER) TOT_RAINnew = 0.
              ENDIF
              IF (TOT_RAINnew <= CLIMIT) THEN
                TOT_RAINnew = 0.
                VRAIN2      = 0.
                QRnew       = 0.
                ARAINnew    = 0.
              ELSE
   !
   !--- 1st guess time-averaged rain rate at bottom of grid box
   !
                RR = TOT_RAINnew/(DTPH*GAMMAR)
   !
   !--- Use same algorithm as above for calculating mean drop diameter
   !      (IDR, in microns), which is used to estimate the time-averaged
   !      fall speed of rain drops at the bottom of the grid layer.  This
   !      isn't perfect, but the alternative is solving a transcendental 
   !      equation that is numerically inefficient and nasty to program
   !      (coded in earlier versions of GSMCOLUMN prior to 8-22-01).
   !
                IF (RR <= RR_DRmin) THEN
                  IDR = MDRmin
                ELSE IF (RR <= RR_DR1) THEN
                  IDR = INT( 1.123E3*RR**.1947 + .5 )
                  IDR = MAX( MDRmin, MIN(IDR, MDR1) )
                ELSE IF (RR <= RR_DR2) THEN
                  IDR = INT( 1.225E3*RR**.2017 + .5 )
                  IDR = MAX( MDR1, MIN(IDR, MDR2) )
                ELSE IF (RR <= RR_DR3) THEN
                  IDR = INT( 1.3006E3*RR**.2083 + .5 )
                  IDR = MAX( MDR2, MIN(IDR, MDR3) )
                ELSE IF (RR <= RR_DRmax) THEN
                  IDR = INT( 1.355E3*RR**.2144 + .5 )
                  IDR = MAX( MDR3, MIN(IDR, MDRmax) )
                ELSE 
                  IDR = MDRmax
                ENDIF              ! End IF (RR <= RR_DRmin)
                VRAIN2 = GAMMAR*VRAIN(IDR)
                QRnew  = TOT_RAINnew / (THICK+BLDTRH*VRAIN2)
                IF (QRnew <= EPSQ) QRnew = 0.
                IF (QR > 0. .AND. QRnew /= 0.) THEN
                  DUM = QRnew / QR
                  IF (DUM < TOLER) QRnew = 0.
                ENDIF
                ARAINnew = BLDTRH*VRAIN2*QRnew
                IF (ARAIN > 0. .AND. ARAINnew /= 0.) THEN
                  DUM = ARAINnew/ARAIN
                  IF (DUM < TOLER) ARAINnew = 0.
                ENDIF
              ENDIF                ! End IF (TOT_RAINnew < CLIMIT)
!
              WCnew = QWnew + QRnew + QInew
!
!----------------------------------------------------------------------
!-------------- Begin debugging & verification ------------------------
!----------------------------------------------------------------------
!
!--- QT, QTnew - total water (vapor & condensate) before & after microphysics, resp.
!
!             QT=THICK*(QV+WC_col(l))+ARAIN+ASNOW
!             QTnew  = THICK*(WVnew/(1.+WVnew)+WCnew/(1.+wcnew))
!    &               + ARAINnew + ASNOWnew

              QT     = THICK*(WV+WC)       + ARAIN    + ASNOW
              QTnew  = THICK*(WVnew+WCnew) + ARAINnew + ASNOWnew
              BUDGET = QT-QTnew
!
!--- Additional check on budget preservation, accounting for truncation effects
!
              DBG_logical=.FALSE.
!             DUM=ABS(BUDGET)
!             IF (DUM > TOLER) THEN
!               DUM=DUM/MIN(QT, QTnew)
!               IF (DUM > TOLER) DBG_logical=.TRUE.
!             ENDIF
!
!             DUM=(RHgrd+.001)*QSInew
!             IF ( (QWnew > EPSQ .OR. QRnew > EPSQ .OR. WVnew > DUM)
!     &           .AND. TC < T_ICE )  DBG_logical=.TRUE.
!
!             IF (TC > 5. .AND. QInewr > EPSQ) DBG_logical=.TRUE.
!
              IF ((WVnew < EPSQ .OR. DBG_logical) .AND. PRINT_diag) THEN
!
                WRITE(6,"(/2(a,i4),2(a,i2))") '{} i=',I_index,' j=',    &
     &                                J_index, ' L=',L,' LSFC=',LSFC
!
                ESW    = min(PP, FPVSL(Tnew))
!               QSWnew = EPS*ESW/(PP-ESW)
                QSWnew = EPS*ESW/(PP+epsm1*ESW)
                IF (TC < 0. .OR. Tnew  <  0.) THEN
                  ESI    = min(PP, FPVSI(Tnew))
!                 QSInew = EPS*ESI/(PP-ESI)
                  QSInew = EPS*ESI/(PP+epsm1*ESI)
                ELSE
                  QSI    = QSW
                  QSInew = QSWnew
                ENDIF
                WSnew = QSInew
                WRITE(6,"(4(a12,g11.4,1x))")                            &
     & '{} TCold=',TC,'TCnew=',Tnew-T0C,'P=',.01*PP,'RHO=',RHO,         &
     & '{} THICK=',THICK,'RHold=',WV/WS,'RHnew=',WVnew/WSnew,           &
     &   'RHgrd=',RHgrd,                                                &
     & '{} RHWold=',WV/QSW,'RHWnew=',WVnew/QSWnew,'RHIold=',WV/QSI,     &
     &   'RHInew=',WVnew/QSInew,                                        &
     & '{} QSWold=',QSW,'QSWnew=',QSWnew,'QSIold=',QSI,'QSInew=',QSInew,&
     & '{} WSold=',WS,'WSnew=',WSnew,'WVold=',WV,'WVnew=',WVnew,        &
     & '{} WCold=',WC,'WCnew=',WCnew,'QWold=',QW,'QWnew=',QWnew,        &
     & '{} QIold=',QI,'QInew=',QInew,'QRold=',QR,'QRnew=',QRnew,        &
     & '{} ARAINold=',ARAIN,'ARAINnew=',ARAINnew,'ASNOWold=',ASNOW,     &
     &   'ASNOWnew=',ASNOWnew,                                          &
     & '{} TOT_RAIN=',TOT_RAIN,'TOT_RAINnew=',TOT_RAINnew,              &
     &   'TOT_ICE=',TOT_ICE,'TOT_ICEnew=',TOT_ICEnew,                   &
     & '{} BUDGET=',BUDGET,'QTold=',QT,'QTnew=',QTnew
!
                WRITE(6,"(4(a12,g11.4,1x))")                            &
     & '{} DELT=',DELT,'DELV=',DELV,'DELW=',DELW,'DELI=',DELI,          &
     & '{} DELR=',DELR,'PCOND=',PCOND,'PIDEP=',PIDEP,'PIEVP=',PIEVP,    &
     & '{} PICND=',PICND,'PREVP=',PREVP,'PRAUT=',PRAUT,'PRACW=',PRACW,  &
     & '{} PIACW=',PIACW,'PIACWI=',PIACWI,'PIACWR=',PIACWR,'PIMLT=',    &
     &    PIMLT,                                                        &
     & '{} PIACR=',PIACR
!
                IF (ICE_logical) WRITE(6,"(4(a12,g11.4,1x))")           &
     & '{} RimeF1=',RimeF1,'GAMMAS=',GAMMAS,'VrimeF=',VrimeF,           &
     &   'VSNOW=',VSNOW,                                                &
     & '{} INDEXS=',FLOAT(INDEXS),'FLARGE=',FLARGE,'FSMALL=',FSMALL,    &
     &   'FLIMASS=',FLIMASS,                                            &
     & '{} XSIMASS=',XSIMASS,'XLIMASS=',XLIMASS,'QLICE=',QLICE,         &
     &   'QTICE=',QTICE,                                                &
     & '{} NLICE=',NLICE,'NSmICE=',NSmICE,'PILOSS=',PILOSS,             &
     &   'EMAIRI=',EMAIRI,                                              &
     & '{} RimeF=',RimeF
!
                IF (TOT_RAIN > 0. .OR. TOT_RAINnew > 0.)                &
     &            WRITE(6,"(4(a12,g11.4,1x))")                          &
     & '{} INDEXR1=',FLOAT(INDEXR1),'INDEXR=',FLOAT(INDEXR),            &
     &   'GAMMAR=',GAMMAR,'N0r=',N0r,                                   &
     & '{} VRAIN1=',VRAIN1,'VRAIN2=',VRAIN2,'QTRAIN=',QTRAIN,'RQR=',RQR,&
     & '{} PRLOSS=',PRLOSS,'VOLR1=',THICK+BLDTRH*VRAIN1,                &
     &   'VOLR2=',THICK+BLDTRH*VRAIN2
!
                IF (PRAUT > 0.) WRITE(6,"(a12,g11.4,1x)") '{} QW0=',QW0
!
                IF (PRACW > 0.) WRITE(6,"(a12,g11.4,1x)") '{} FWR=',FWR
!
                IF (PIACR > 0.) WRITE(6,"(a12,g11.4,1x)") '{} FIR=',FIR
!
                DUM = PIMLT + PICND - PREVP - PIEVP
                IF (DUM > 0. .or. DWVi /= 0.)                           &
     &            WRITE(6,"(4(a12,g11.4,1x))")                          &
     & '{} TFACTOR=',TFACTOR,'DYNVIS=',DYNVIS,                          &
     &   'THERM_CON=',THERM_COND,'DIFFUS=',DIFFUS
!
                IF (PREVP < 0.) WRITE(6,"(4(a12,g11.4,1x))")            &
     & '{} RFACTOR=',RFACTOR,'ABW=',ABW,'VENTR=',VENTR,'CREVP=',CREVP,  &
     & '{} DWVr=',DWVr,'DENOMW=',DENOMW
!
                IF (PIDEP /= 0. .AND. DWVi /= 0.)                       &
     &            WRITE(6,"(4(a12,g11.4,1x))")                          &
     & '{} DWVi=',DWVi,'DENOMI=',DENOMI,'PIDEP_max=',PIDEP_max,         &
     &   'SFACTOR=',SFACTOR,                                            &
     & '{} ABI=',ABI,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),        &
     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                             &
     & '{} VENTIS=',VENTIS,'DIDEP=',DIDEP
!
                IF (PIDEP > 0. .AND. PCOND /= 0.)                       &
     &            WRITE(6,"(4(a12,g11.4,1x))")                          &
     & '{} DENOMW=',DENOMW,'DENOMWI=',DENOMWI,'DENOMF=',DENOMF,         &
     &    'DUM2=',PCOND-PIACW
!
                IF (FWS > 0.) WRITE(6,"(4(a12,g11.4,1x))") '{} FWS=',FWS
!
                DUM = PIMLT + PICND - PIEVP
                IF (DUM >  0.) WRITE(6,"(4(a12,g11.4,1x))")             &
     & '{} SFACTOR=',SFACTOR,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),&
     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                             &
     & '{} AIEVP=',AIEVP,'DIEVP=',DIEVP,'QSW0=',QSW0,'DWV0=',DWV0
   !
              ENDIF
!
!----------------------------------------------------------------------
!-------------- Water budget statistics & maximum values --------------
!----------------------------------------------------------------------
!
              IF (PRINT_diag) THEN
                ITdx = MAX( ITLO, MIN( INT(Tnew-T0C), ITHI ) )
                IF (QInew > EPSQ) NSTATS(ITdx,1) = NSTATS(ITdx,1)+1
                IF (QInew > EPSQ .AND.  QRnew+QWnew > EPSQ)             &
     &            NSTATS(ITdx,2) = NSTATS(ITdx,2)+1
                IF (QWnew > EPSQ) NSTATS(ITdx,3) = NSTATS(ITdx,3)+1 
                IF (QRnew > EPSQ) NSTATS(ITdx,4) = NSTATS(ITdx,4)+1
  !
                QMAX(ITdx,1)  = MAX(QMAX(ITdx,1), QInew)
                QMAX(ITdx,2)  = MAX(QMAX(ITdx,2), QWnew)
                QMAX(ITdx,3)  = MAX(QMAX(ITdx,3), QRnew)
                QMAX(ITdx,4)  = MAX(QMAX(ITdx,4), ASNOWnew)
                QMAX(ITdx,5)  = MAX(QMAX(ITdx,5), ARAINnew)
                QTOT(ITdx,1)  = QTOT(ITdx,1)+QInew*THICK
                QTOT(ITdx,2)  = QTOT(ITdx,2)+QWnew*THICK
                QTOT(ITdx,3)  = QTOT(ITdx,3)+QRnew*THICK
  !
                QTOT(ITdx,4)  = QTOT(ITdx,4)+PCOND*THICK
                QTOT(ITdx,5)  = QTOT(ITdx,5)+PICND*THICK
                QTOT(ITdx,6)  = QTOT(ITdx,6)+PIEVP*THICK
                QTOT(ITdx,7)  = QTOT(ITdx,7)+PIDEP*THICK
                QTOT(ITdx,8)  = QTOT(ITdx,8)+PREVP*THICK
                QTOT(ITdx,9)  = QTOT(ITdx,9)+PRAUT*THICK
                QTOT(ITdx,10) = QTOT(ITdx,10)+PRACW*THICK
                QTOT(ITdx,11) = QTOT(ITdx,11)+PIMLT*THICK
                QTOT(ITdx,12) = QTOT(ITdx,12)+PIACW*THICK
                QTOT(ITdx,13) = QTOT(ITdx,13)+PIACWI*THICK
                QTOT(ITdx,14) = QTOT(ITdx,14)+PIACWR*THICK
                QTOT(ITdx,15) = QTOT(ITdx,15)+PIACR*THICK
  !
                QTOT(ITdx,16) = QTOT(ITdx,16)+(WVnew-WV)*THICK
                QTOT(ITdx,17) = QTOT(ITdx,17)+(QWnew-QW)*THICK
                QTOT(ITdx,18) = QTOT(ITdx,18)+(QInew-QI)*THICK
                QTOT(ITdx,19) = QTOT(ITdx,19)+(QRnew-QR)*THICK
                QTOT(ITdx,20) = QTOT(ITdx,20)+(ARAINnew-ARAIN)
                QTOT(ITdx,21) = QTOT(ITdx,21)+(ASNOWnew-ASNOW)
                IF (QInew > 0.)                                         &
     &            QTOT(ITdx,22) = QTOT(ITdx,22)+QInew*THICK/RimeF
  !
              ENDIF
!
!----------------------------------------------------------------------
!------------------------- Update arrays ------------------------------
!----------------------------------------------------------------------
!
              T_col(L)     = Tnew                        ! temperature
!
!             QV_col(L)    = max(EPSQ, WVnew/(1.+WVnew)) ! specific humidity
              QV_col(L)    = max(0.0, WVnew           )  ! specific humidity
              WC_col(L)    = max(0.0, WCnew)             ! total condensate mixing ratio
              QI_col(L)    = max(0.0, QInew)             ! ice mixing ratio
              QR_col(L)    = max(0.0, QRnew)             ! rain mixing ratio
              QW_col(L)    = max(0.0, QWnew)             ! cloud water mixing ratio
              RimeF_col(L) = RimeF                       ! rime factor
              ASNOW        = ASNOWnew                    ! accumulated snow
              ARAIN        = ARAINnew                    ! accumulated rain
!
!#######################################################################
!
            ENDIF   ! End of IF (.NOT. CLEAR) THEN
          ENDIF     ! End of IF (QV_col(L) <= EPSQ .AND. WC_col(L) <= EPSQ) THEN
!
        ENDDO       ! ##### End "L" loop through model levels #####
!
        ARAING = ARAING + ARAIN
        ASNOWG = ASNOWG + ASNOW
      enddo              ! do for ntimes=1,mic_step
!
!#######################################################################
!
!-----------------------------------------------------------------------
!--------------------------- Return to GSMDRIVE -----------------------
!-----------------------------------------------------------------------
!
      CONTAINS
!     END  SUBROUTINE GSMCOLUMN
!
!#######################################################################
!--------- Produces accurate calculation of cloud condensation ---------
!#######################################################################
!
      REAL FUNCTION CONDENSE (PP, QW, RHgrd, TK, WV)
!
      implicit none
!
!---------------------------------------------------------------------------------
!------ The Asai (1965) algorithm takes into consideration the release of ------
!------ latent heat in increasing the temperature & in increasing the     ------
!------ saturation mixing ratio (following the Clausius-Clapeyron eqn.).  ------
!---------------------------------------------------------------------------------
!
      real pp, qw, rhgrd, tk, wv
      INTEGER, PARAMETER :: HIGH_PRES=kind_phys
!     INTEGER, PARAMETER :: HIGH_PRES=Selected_Real_Kind(15)
      REAL (KIND=HIGH_PRES), PARAMETER ::                               &
     & RHLIMIT=.001, RHLIMIT1=-RHLIMIT
      REAL, PARAMETER :: RCP=1./CP, RCPRV=RCP/RV
      REAL (KIND=HIGH_PRES) :: COND, SSAT, WCdum, tsq
      real wvdum, tdum, xlv, xlv1, xlv2, ws, dwv, esw, rfac
!
!-----------------------------------------------------------------------
!
!--- LV (T) is from Bolton (JAS, 1980)
!
!     XLV=3.148E6-2370.*TK
!     XLV1=XLV*RCP
!     XLV2=XLV*XLV*RCPRV
!
      Tdum     = TK
      WVdum    = WV
      WCdum    = QW
      ESW      = min(PP, FPVSL(Tdum))          ! Saturation vapor press w/r/t water
!     WS       = RHgrd*EPS*ESW/(PP-ESW)        ! Saturation mixing ratio w/r/t water
      WS       = RHgrd*EPS*ESW/(PP+epsm1*ESW)  ! Saturation specific hum w/r/t water
      DWV      = WVdum - WS                    ! Deficit grid-scale specific humidity
      SSAT     = DWV / WS                      ! Supersaturation ratio
      CONDENSE = 0.
      rfac     = 0.5                           ! converges faster with 0.5
      DO WHILE ((SSAT < RHLIMIT1 .AND. WCdum > EPSQ)                    &
     &           .OR. SSAT > RHLIMIT)
!
        XLV   = 3.148E6-2370.*Tdum
        XLV1  = XLV*RCP
        XLV2  = XLV*XLV*RCPRV
!
!       COND  = DWV/(1.+XLV2*WS/(Tdum*Tdum))  ! Asai (1965, J. Japan)
        tsq   = Tdum*Tdum
        COND  = rfac*DWV*tsq/(tsq+XLV2*WS)    ! Asai (1965, J. Japan)
!       COND  =      DWV*tsq/(tsq+XLV2*WS)    ! Asai (1965, J. Japan)
        COND  = MAX(COND, -WCdum)             ! Limit cloud water evaporation
        Tdum  = Tdum+XLV1*COND                ! Updated temperature
        WVdum = WVdum-COND                    ! Updated water vapor mixing ratio
        WCdum = WCdum+COND                    ! Updated cloud water mixing ratio
        CONDENSE = CONDENSE + COND            ! Total cloud water condensation
        ESW  = min(PP, FPVSL(Tdum))           ! Updated saturation vapor press w/r/t water
!       WS   = RHgrd*EPS*ESW/(PP-ESW)         ! Updated saturation mixing ratio w/r/t water
        WS   = RHgrd*EPS*ESW/(PP+epsm1*ESW)   ! Updated saturation mixing ratio w/r/t water
        DWV  = WVdum-WS                       ! Deficit grid-scale water vapor mixing ratio
        SSAT = DWV / WS                       ! Grid-scale supersaturation ratio
        rfac = 1.0
      ENDDO

      END FUNCTION CONDENSE
!
!#######################################################################
!---------------- Calculate ice deposition at T<T_ICE ------------------
!#######################################################################
!
      REAL FUNCTION DEPOSIT (PP, RHgrd, Tdum, WVdum)
!
      implicit none
!
!--- Also uses the Asai (1965) algorithm, but uses a different target
!      vapor pressure for the adjustment
!
      REAL PP, RHgrd, Tdum, WVdum
      INTEGER, PARAMETER :: HIGH_PRES=kind_phys
!     INTEGER, PARAMETER :: HIGH_PRES=Selected_Real_Kind(15)
      REAL (KIND=HIGH_PRES), PARAMETER :: RHLIMIT=.001,                 & 
     & RHLIMIT1=-RHLIMIT
      REAL, PARAMETER :: RCP=1./CP, RCPRV=RCP/RV, XLS=HVAP+HFUS         &
     &,                  XLS1=XLS*RCP, XLS2=XLS*XLS*RCPRV
      REAL (KIND=HIGH_PRES) :: DEP, SSAT
      real esi, ws, dwv
!
!-----------------------------------------------------------------------
!
      ESI=min(PP, FPVSI(Tdum))                  ! Saturation vapor press w/r/t ice
!     WS=RHgrd*EPS*ESI/(PP-ESI)                 ! Saturation mixing ratio
      WS=RHgrd*EPS*ESI/(PP+epsm1*ESI)           ! Saturation mixing ratio
      DWV=WVdum-WS                              ! Deficit grid-scale water vapor mixing ratio
      SSAT=DWV/WS                               ! Supersaturation ratio
      DEPOSIT=0.
      DO WHILE (SSAT > RHLIMIT .OR. SSAT < RHLIMIT1)
   !
   !--- Note that XLVS2=LS*LV/(CP*RV)=LV*WS/(RV*T*T)*(LS/CP*DEP1),
   !     where WS is the saturation mixing ratio following Clausius-
   !     Clapeyron (see Asai,1965; Young,1993,p.405)
   !
        DEP=DWV/(1.+XLS2*WS/(Tdum*Tdum))        ! Asai (1965, J. Japan)
        Tdum=Tdum+XLS1*DEP                      ! Updated temperature
        WVdum=WVdum-DEP                         ! Updated ice mixing ratio
        DEPOSIT=DEPOSIT+DEP                     ! Total ice deposition
        ESI=min(PP, FPVSI(Tdum))                ! Updated saturation vapor press w/r/t ice
!       WS=RHgrd*EPS*ESI/(PP-ESI)               ! Updated saturation mixing ratio w/r/t ice
        WS=RHgrd*EPS*ESI/(PP+epsm1*ESI)         ! Updated saturation mixing ratio w/r/t ice
        DWV=WVdum-WS                            ! Deficit grid-scale water vapor mixing ratio
        SSAT=DWV/WS                             ! Grid-scale supersaturation ratio
      ENDDO
      END FUNCTION DEPOSIT
!
      END  SUBROUTINE GSMCOLUMN


      SUBROUTINE rsipath(im, ix, ix2, levs, prsl, prsi, t, q, clw       &
     &,                  f_ice, f_rain, f_rime, flgmin                  &
     &,                  cwatp, cicep, rainp, snowp                     &
     &,                  recwat, rerain, resnow, lprnt, ipr)
!
      implicit none
!
!--------------------CLOUD----------------------------------------------
      integer im, ix, ix2, levs, ipr
      real    prsl(ix,levs), prsi(ix,levs+1), t(ix,levs), q(ix,levs)    &
     &,       clw(ix2,levs), f_ice(ix2,levs), f_rain(ix2,levs)          &
     &,       f_rime(ix2,levs)                                          &
     &,       cwatp(ix,levs), rainp(ix,levs),  cicep(ix,levs)           &
     &,       snowp(ix,levs), recwat(ix,levs), resnow(ix,levs)          &
     &,       rerain(ix,levs)
      real    flgmin
      real    frice, frrain, qcice, qcwat,  qrain, qsnow,   qtot, sden  &
     &,       cpath, rho,    dsnow, flarge, rimef, xsimass, nlice       &
     &,       tc,    recw1,  drain, xli,    dum,   NLImax, pfac, pp     &
     &,       snofac, tem
!
      real, parameter :: cexp=1./3.
      integer i, l, indexs
      logical lprnt
!

      RECW1 = 620.3505 / TNW**CEXP         ! cloud droplet effective radius

      do l=1,levs
        do i=1,im
                                           !--- HYDROMETEOR'S OPTICAL PATH
           CWATP(I,L) = 0.
           CICEP(I,L) = 0.
           RAINP(I,L) = 0.
           SNOWP(I,L) = 0.
                                           !--- HYDROMETEOR'S EFFECTIVE RADIUS
           RECWAT(I,L) = RECWmin
           RERAIN(I,L) = RERAINmin
           RESNOW(I,L) = RESNOWmin
        ENDDO
      ENDDO

      do l=1,levs
        DO I=1,im

   !        Assume parameterized condensate is 
   !         all water for T>=-10C,
   !         all ice for T<=-30C,
   !         and a linear mixture at -10C > T > -30C
   !
   !    * Determine hydrometeor composition of total condensate (QTOT)
   !
!         pp   = prsl(i,l) * 1000.0
          pp   = prsl(i,l) / prsi(i,levs+1)
!         pfac = max(0.25, sqrt(sqrt(min(1.0, pp*0.000025))))
!         pfac = max(0.5, sqrt(sqrt(min(1.0, pp*0.000025))))
!         pfac = max(0.5, sqrt(sqrt(min(1.0, pp*0.00002))))
!         pfac = max(0.25, sqrt(sqrt(min(1.0, pp*0.00001))))
!         pfac = max(0.25, sqrt(sqrt(min(1.0, pp))))
!         pfac = max(0.1, sqrt(min(1.0, pp*0.00001)))
!         pfac = max(0.5, sqrt(sqrt(min(1.0, pp*0.000033))))
!         pfac = max(0.5, sqrt(sqrt(min(1.0, pp*0.00004))))
!go       pfac = max(0.5, (sqrt(min(1.0, pp*0.000025))))
          pfac = 1.0
          TC   = T(I,L) - t0c
          QTOT = clw(I,L)
          IF (QTOT > EPSQ) THEN
             QCWAT=0.
             QCICE=0.
             QRAIN=0.
             QSNOW=0.
             FRice  = max(0.0, min(1.0, F_ice(I,L)))
             FRrain = max(0.0, min(1.0, F_rain(I,L)))
             IF(TC <= Thom) then
                QCICE = QTOT
             ELSE
                QCICE = FRice  * QTOT
                QCWAT = QTOT   - QCICE
                QRAIN = FRrain * QCWAT
                QCWAT = QCWAT  - QRAIN
             ENDIF
    !
    !--- Air density (RHO), model mass thickness (CPATH)
    !
             RHO   = PRSL(I,L)/(RD*T(I,L)*(1.+EPS1*Q(I,L)))
             CPATH = (PRSI(I,L+1)-PRSI(I,L))*(1000000.0/grav)

    !! CLOUD WATER
    !
    !--- Effective radius (RECWAT) & total water path (CWATP)
    !    Assume monodisperse distribution of droplets (no factor of 1.5)
    !
             IF(QCWAT > 0.) THEN
                RECWAT(I,L) = MAX(RECWmin, RECW1*(RHO*QCWAT)**CEXP)
                CWATP(I,L)  = CPATH*QCWAT         ! cloud water path
!               tem         = 5.0*(1+max(0.0,min(1.0,-0.05*tc)))
!               RECWAT(I,L) = max(RECWAT(I,L), tem)
             ENDIF

    !! RAIN
    !
    !--- Effective radius (RERAIN) & total water path (RAINP)
    !--- Factor of 1.5 accounts for r**3/r**2 moments for exponentially
    !    distributed drops in effective radius calculations
    !    (from M.D. Chou's code provided to Y.-T. Hou)
    !
             IF(QRAIN > 0.) THEN
                DRAIN       = CN0r0*sqrt(sqrt((RHO*QRAIN)))
                RERAIN(I,L) = 1.5*MAX(XMRmin, MIN(DRAIN, XMRmax))
                RAINP(I,L)  = CPATH*QRAIN         ! rain water path
             ENDIF

    !! SNOW (large ice) & CLOUD ICE
    !
    !--- Effective radius (RESNOW) & total ice path (SNOWP)
    !--- Total ice path (CICEP) for cloud ice 
    !--- Factor of 1.5 accounts for r**3/r**2 moments for exponentially
    !    distributed ice particles in effective radius calculations 
    !
    !--- Separation of cloud ice & "snow" uses algorithm from
    !    subroutine GSMCOLUMN
    !
             IF(QCICE > 0.) THEN
    !
    !--- Mean particle size following Houze et al. (JAS, 1979, p. 160), 
    !    converted from Fig. 5 plot of LAMDAs.  An analogous set of
    !    relationships also shown by Fig. 8 of Ryan (BAMS, 1996, p. 66),
    !    but with a variety of different relationships that parallel the
    !    Houze curves.
    !
!               DUM=MAX(0.05, MIN(1., EXP(.0536*TC)) )
                DUM=MAX(0.05, MIN(1., EXP(.0564*TC)) )
                INDEXS=MIN(MDImax, MAX(MDImin, INT(XMImax*DUM) ) )
!               indexs=max(INDEXSmin, indexs)
!               NLImax=5.E3/sqrt(DUM)        !- Ver3
                DUM=MAX(FLGmin*pfac, DUM)
!               DUM=MAX(FLGmin, DUM)
!               NLImax=20.E3            !- Ver3
!               NLImax=50.E3            !- Ver3 => comment this line out
                NLImax=10.E3/sqrt(DUM)        !- Ver3
!               NLImax=5.E3/sqrt(DUM)        !- Ver3
!               NLImax=6.E3/sqrt(DUM)        !- Ver3
!               NLImax=7.5E3/sqrt(DUM)       !- Ver3
!               NLImax=20.E3/DUM        !- Ver3
!               NLImax=20.E3/max(0.2,DUM)        !- Ver3
!               NLImax=2.0E3/max(0.1,DUM)        !- Ver3
!               NLImax=2.5E3/max(0.1,DUM)        !- Ver3
!               NLImax=10.E3/max(0.2,DUM)        !- Ver3
!               NLImax=4.E3/max(0.2,DUM)        !- Ver3
!Moorthi        DSNOW  = XMImax*EXP(.0536*TC)
!Moorthi        INDEXS = MAX(INDEXSmin, MIN(MDImax, INT(DSNOW)))
    !
    !--- Assumed number fraction of large ice to total (large & small)
    !    ice particles, which is based on a general impression of the
    !    literature.
    !
    !    Small ice are assumed to have a mean diameter of 50 microns.
    !
                IF(TC >= 0.) THEN
                  FLARGE=FLG1P0
                ELSE
                  FLARGE = dum
                ENDIF
!------------------------Commented by Moorthi -----------------------------
!               ELSEIF (TC >= -25.) THEN
!
!--- Note that absence of cloud water (QCWAT) is used as a quick
!    substitute for calculating water subsaturation as in GSMCOLUMN
!
!                  IF(QCWAT <= 0. .OR. TC < -8.
!    &                            .OR. TC > -3.)THEN
!                     FLARGE=FLG0P2
!                  ELSE
      
!--- Parameterize effects of rime splintering by increasing
!    number of small ice particles
!
!                     FLARGE=FLG0P1
!                  ENDIF
!               ELSEIF (TC <= -50.) THEN
!                  FLARGE=.01
!               ELSE
!                  FLARGE=.2*EXP(.1198*(TC+25))
!               ENDIF
!____________________________________________________________________________

                RimeF=MAX(1., F_RIME(I,L))
                XSIMASS=MASSI(MDImin)*(1.-FLARGE)/FLARGE
!     if (lprnt) print *,' rimef=',rimef,' xsimass=',xsimass
!    &,' indexs=',indexs,' massi=',massi(indexs),' flarge=',flarge
                NLICE=RHO*QCICE/(XSIMASS+RimeF*MASSI(INDEXS))
    !
    !--- From subroutine GSMCOLUMN:
    !--- Minimum number concentration for large ice of NLImin=10/m**3 
    !    at T>=0C.  Done in order to prevent unrealistically small 
    !    melting rates and tiny amounts of snow from falling to 
    !    unrealistically warm temperatures.
    !
                IF(TC >= 0.) THEN
                   NLICE=MAX(NLImin, NLICE)
                ELSEIF (NLICE > NLImax) THEN
      !
      !--- Ferrier 6/13/01:  Prevent excess accumulation of ice
      !
                   XLI=(RHO*QCICE/NLImax-XSIMASS)/RimeF
 
                   IF(XLI <= MASSI(450) ) THEN
                      DSNOW=9.5885E5*XLI**.42066
                   ELSE
                      DSNOW=3.9751E6*XLI**.49870
                   ENDIF
 
                   INDEXS=MIN(MDImax, MAX(INDEXS, INT(DSNOW)))
                   NLICE=RHO*QCICE/(XSIMASS+RimeF*MASSI(INDEXS))
                ENDIF

!               if (tc > -20.0 .and. indexs >= indexsmin) then
!                 snofac = max(0.0, min(1.0, exp(1.0*(tc+20.0))))
!               if (indexs >= indexsmin) then
!               if (tc > -20.0 .or. indexs >= indexsmin) then
!               if (tc > -40.0) then
!               if (tc >= -40.0 .or. prsl(i,l) > 50.0) then
!!              if (tc >= -20.0) then
!               if (tc >= -20.0 .or. prsl(i,l) > 50.0) then
!               if ((tc >= -20.0 .or.
!    &              prsi(i,levs+1)-prsi(i,l) < 30.0)
                if (prsi(i,levs+1)-prsi(i,l) < 40.0                     &
!               if (prsi(i,levs+1)-prsi(i,l) < 70.0
     &              .and. indexs >= indexsmin) then
!    &              prsi(i,levs)-prsl(i,l) < 20.0) then
!    &              prsi(i,levs)-prsl(i,l) < 30.0) then
!    &              prsi(i,levs)-prsl(i,l) < 40.0) then
!                 snofac = max(0.0, min(1.0, 0.05*(tc+40.0)))
!                 snofac = max(0.0, min(1.0, 0.1*(tc+25.0)))
!                 snofac = max(0.0, min(1.0, 0.0667*(tc+25.0)))
!               if (indexs > indexsmin) then
                  QSNOW     = MIN(QCICE, NLICE*RimeF*MASSI(INDEXS)/RHO)
!    &                      * snofac
                endif
!               qsnow       = qcice
                QCICE       = MAX(0., QCICE-QSNOW)
!               qsnow       = 0.0
                CICEP (I,L) = CPATH*QCICE          ! cloud ice path
                RESNOW(I,L) = 1.5*FLOAT(INDEXS)
                SDEN        = SDENS(INDEXS)/RimeF  ! 1/snow density
                SNOWP (I,L) = CPATH*QSNOW*SDEN     ! snow path / snow density
!               SNOWP (I,L) = CPATH*QSNOW          ! snow path / snow density
!               if (lprnt .and. i == ipr) then
!                 print *,' L=',L,' snowp=',snowp(i,l),' cpath=',cpath
!    &,' qsnow=',qsnow,' sden=',sden,' rimef=',rimef,' indexs=',indexs
!    &,' sdens=',sdens(indexs),' resnow=',resnow(i,l)
!    &,' qcice=',qcice,' cicep=',cicep(i,l)
!               endif

                
             ENDIF                                 ! END QCICE BLOCK
           ENDIF                                   ! QTOT IF BLOCK

        ENDDO
      ENDDO
!
      END SUBROUTINE rsipath



!-----------------------------------
      subroutine rsipath2                                               &
!...................................

!  ---  inputs:
     &     ( plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime,        &
     &       IM, LEVS, iflip, flgmin,                                   &
!  ---  outputs:
     &       cwatp, cicep, rainp, snowp, recwat, rerain, resnow, snden  &
     &     )

! =================   subprogram documentation block   ================ !
!                                                                       !
! abstract:  this program is a modified version of ferrier's original   !
!   "rsipath" subprogram.  it computes layer's cloud liquid, ice, rain, !
!   and snow water condensate path and the partical effective radius    !
!   for liquid droplet, rain drop, and snow flake.                      !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IM,LEVS) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IM,LEVS+1):model level pressure in mb (100Pa)                !
!   tlyr  (IM,LEVS) : model layer mean temperature in k                 !
!   qlyr  (IM,LEVS) : layer specific humidity in gm/gm                  !
!   qcwat (IM,LEVS) : layer cloud liquid water condensate amount        !
!   qcice (IM,LEVS) : layer cloud ice water condensate amount           !
!   qrain (IM,LEVS) : layer rain drop water amount                      !
!   rrime (IM,LEVS) : mass ratio of total to unrimed ice ( >= 1 )       !
!   IM              : horizontal dimention                              !
!   LEVS            : vertical layer dimensions                         !
!   iflip           : control flag for in/out vertical indexing         !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   flgmin          : Minimum large ice fraction                        !
!   lprnt           : logical check print control flag                  !
!                                                                       !
! output variables:                                                     !
!   cwatp (IM,LEVS) : layer cloud liquid water path                     !
!   cicep (IM,LEVS) : layer cloud ice water path                        !
!   rainp (IM,LEVS) : layer rain water path                             !
!   snowp (IM,LEVS) : layer snow water path                             !
!   recwat(IM,LEVS) : layer cloud eff radius for liqid water (micron)   !
!   rerain(IM,LEVS) : layer rain water effective radius      (micron)   !
!   resnow(IM,LEVS) : layer snow flake effective radius      (micron)   !
!   snden (IM,LEVS) : 1/snow density                                    !
!                                                                       !
!                                                                       !
! usage:     call rsipath2                                              !
!                                                                       !
! subroutines called:  none                                             !
!                                                                       !
! program history log:                                                  !
!      xx-xx-2001   b. ferrier     - original program                   !
!      xx-xx-2004   s. moorthi     - modified for use in gfs model      !
!      05-20-2004   y. hou         - modified, added vertical index flag!
!                     to reduce data flipping, and rearrange code to    !
!                     be comformable with radiation part programs.      !
!                                                                       !
!  ====================    end of description    =====================  !
!

      implicit none

!  ---  constant parameter:
      real, parameter :: CEXP= 1.0/3.0

!  ---  inputs:
      real, dimension(:,:), intent(in) ::                               &
     &       plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime

      integer, intent(in) :: IM, LEVS, iflip
      real, dimension(:),   intent(in) :: flgmin
!     logical, intent(in) :: lprnt

!  ---  output:
      real, dimension(:,:), intent(out) ::                              &
     &       cwatp, cicep, rainp, snowp, recwat, rerain, resnow, snden

!  ---  locals:
!     real,    dimension(IM,LEVS) :: delp, pp1, pp2

      real    :: recw1, dsnow, qsnow, qqcice, flarge, xsimass, pfac,    &
     &           nlice, xli, nlimax, dum, tem,                          &
     &           rho, cpath, rc, totcnd, tc

      integer :: i, k, indexs, ksfc, k1
!
!===>  ...  begin here
!
      recw1 = 620.3505 / TNW**CEXP         ! cloud droplet effective radius

      do k = 1, LEVS
        do i = 1, IM
                                           !--- hydrometeor's optical path
           cwatp(i,k) = 0.0
           cicep(i,k) = 0.0
           rainp(i,k) = 0.0
           snowp(i,k) = 0.0
           snden(i,k) = 0.0
                                           !--- hydrometeor's effective radius
           recwat(i,k) = RECWmin
           rerain(i,k) = RERAINmin
           resnow(i,k) = RESNOWmin
        enddo
      enddo

!  ---  set up pressure related arrays, convert unit from mb to cb (10Pa)
!       cause the rest part uses cb in computation

      if (iflip == 0) then        ! data from toa to sfc
        ksfc = levs + 1
        k1   = 0
      else                        ! data from sfc to top
        ksfc = 1
        k1   = 1
      endif                       ! end_if_iflip
!
      do k = 1, LEVS
        do i = 1, IM
          totcnd = qcwat(i,k) + qcice(i,k) + qrain(i,k)
          qsnow = 0.0
          if(totcnd > EPSQ) then

!  ---  air density (rho), model mass thickness (cpath), temperature in c (tc)

            rho   = 0.1 * plyr(i,k)                                     &
     &            / (RD* tlyr(i,k) * (1.0 + EPS1*qlyr(i,k)))
            cpath = abs(plvl(i,k+1) - plvl(i,k)) * (100000.0 / GRAV)
            tc    = tlyr(i,k) - T0C

!! cloud water
!
!  ---  effective radius (recwat) & total water path (cwatp):
!       assume monodisperse distribution of droplets (no factor of 1.5)

            if (qcwat(i,k) > 0.0) then
              recwat(i,k) = max(RECWmin,recw1*(rho*qcwat(i,k))**CEXP)
              cwatp (i,k) = cpath * qcwat(i,k)           ! cloud water path
!             tem         = 5.0*(1.0 + max(0.0, min(1.0,-0.05*tc)))
!             recwat(i,k) = max(recwat(i,k), tem)
            endif

!! rain
!
!  ---  effective radius (rerain) & total water path (rainp):
!       factor of 1.5 accounts for r**3/r**2 moments for exponentially
!       distributed drops in effective radius calculations
!       (from m.d. chou's code provided to y.-t. hou)

            if (qrain(i,k) > 0.0) then
              tem         = CN0r0 * sqrt(sqrt(rho*qrain(i,k)))
              rerain(i,k) = 1.5 * max(XMRmin, min(XMRmax, tem))
              rainp (i,k) = cpath * qrain(i,k)           ! rain water path
            endif

!! snow (large ice) & cloud ice
!
!  ---  effective radius (resnow) & total ice path (snowp) for snow, and
!       total ice path (cicep) for cloud ice:
!       factor of 1.5 accounts for r**3/r**2 moments for exponentially
!       distributed ice particles in effective radius calculations
!       separation of cloud ice & "snow" uses algorithm from subroutine gsmcolumn

!           pfac = max(0.5, sqrt(sqrt(min(1.0, pp1(i,k)*0.00004))))
!go         pfac = max(0.5, (sqrt(min(1.0, pp1(i,k)*0.000025))))
            pfac = 1.0

            if (qcice(i,k) > 0.0) then

!  ---  mean particle size following houze et al. (jas, 1979, p. 160),
!       converted from fig. 5 plot of lamdas.  an analogous set of
!       relationships also shown by fig. 8 of ryan (bams, 1996, p. 66),
!       but with a variety of different relationships that parallel
!       the houze curves.

!             dum = max(0.05, min(1.0, exp(0.0536*tc) ))
              dum = max(0.05, min(1.0, exp(0.0564*tc) ))
              indexs = min(MDImax, max(MDImin, int(XMImax*dum) ))
              DUM=MAX(FLGmin(i)*pfac, DUM)

!  ---  assumed number fraction of large ice to total (large & small) ice
!       particles, which is based on a general impression of the literature.
!       small ice are assumed to have a mean diameter of 50 microns.

              if (tc >= 0.0) then
                flarge = FLG1P0
              else
                flarge = dum
!               flarge = max(FLGmin*pfac, dum)
              endif
!------------------------commented by moorthi -----------------------------
!             elseif (tc >= -25.0) then
!
!  ---  note that absence of cloud water (qcwat) is used as a quick
!       substitute for calculating water subsaturation as in gsmcolumn
!
!               if (qcwat(i,k) <= 0.0 .or. tc < -8.0                 &
!    &                                .or. tc > -3.0) then
!                 flarge = FLG0P2
!               else
!
!  ---  parameterize effects of rime splintering by increasing
!       number of small ice particles
!
!                 flarge = FLG0P1
!               endif
!             elseif (tc <= -50.0) then
!               flarge = 0.01
!             else
!               flarge = 0.2 * exp(0.1198*(tc+25.0))
!             endif
!____________________________________________________________________________

              xsimass = MASSI(MDImin) * (1.0 - flarge) / flarge
!             nlimax = 20.0e3                                      !- ver3
!             NLImax=50.E3                 !- Ver3 => comment this line out
              NLImax=10.E3/sqrt(DUM)       !- Ver3
!             NLImax=5.E3/sqrt(DUM)        !- Ver3
!             NLImax=6.E3/sqrt(DUM)        !- Ver3
!             NLImax=7.5E3/sqrt(DUM)       !- Ver3

!             indexs = min(MDImax, max(MDImin, int(XMImax*dum) ))
!moorthi      dsnow  = XMImax * exp(0.0536*tc)
!moorthi      indexs = max(INDEXSmin, min(MDImax, int(dsnow)))

!             if (lprnt) print *,' rrime=',rrime,' xsimass=',xsimass,   &
!    &       ' indexs=',indexs,' massi=',massi(indexs),' flarge=',flarge

              tem = rho * qcice(i,k)
              nlice = tem / (xsimass +rrime(i,k)*MASSI(indexs))

!  ---  from subroutine gsmcolumn:
!       minimum number concentration for large ice of NLImin=10/m**3
!       at t>=0c.  done in order to prevent unrealistically small
!       melting rates and tiny amounts of snow from falling to
!       unrealistically warm temperatures.

              if (tc >= 0.0) then

                nlice = max(NLImin, nlice)

              elseif (nlice > nlimax) then

!  ---  ferrier 6/13/01:  prevent excess accumulation of ice

                xli = (tem/nlimax - xsimass) / rrime(i,k)

                if (xli <= MASSI(450) ) then
                  dsnow = 9.5885e5 * xli**0.42066
                else
                  dsnow = 3.9751e6 * xli** 0.49870
                endif

                indexs = min(MDImax, max(indexs, int(dsnow)))
                nlice = tem / (xsimass + rrime(i,k)*MASSI(indexs))

              endif                               ! end if_tc block

!             if (abs(plvl(i,ksfc)-plvl(i,k+k1)) < 300.0                &
!             if (abs(plvl(i,ksfc)-plvl(i,k+k1)) < 400.0                &
!             if (plvl(i,k+k1) > 600.0                                  &
!    &                            .and. indexs >= INDEXSmin) then
!             if (tc > -20.0 .and. indexs >= indexsmin) then
              if (plvl(i,ksfc) > 850.0 .and.                            &
!    &            plvl(i,k+k1) > 600.0 .and. indexs >= indexsmin) then
     &            plvl(i,k+k1) > 700.0 .and. indexs >= indexsmin) then ! 20060516
!!            if (plvl(i,ksfc) > 800.0 .and.                            &
!!   &            plvl(i,k+k1) > 700.0 .and. indexs >= indexsmin) then
!             if (plvl(i,ksfc) > 700.0 .and.                            &
!    &            plvl(i,k+k1) > 600.0 .and. indexs >= indexsmin) then
                qsnow = min( qcice(i,k),                                &
     &                       nlice*rrime(i,k)*MASSI(indexs)/rho )
              endif

              qqcice      = max(0.0, qcice(i,k)-qsnow)
              cicep (i,k) = cpath * qqcice          ! cloud ice path
              resnow(i,k) = 1.5 * float(indexs)
              snden (i,k) = SDENS(indexs) / rrime(i,k)   ! 1/snow density
              snowp (i,k) = cpath*qsnow             ! snow path
!             snowp (i,k) = cpath*qsnow*snden(i,k)  ! snow path / snow density

!             if (lprnt .and. i == ipr) then
!             if (i == 2) then
!               print *,' L=',k,' snowp=',snowp(i,k),' cpath=',cpath,   &
!    &         ' qsnow=',qsnow,' sden=',snden(i,k),' rrime=',rrime(i,k),&
!    &         ' indexs=',indexs,' sdens=',sdens(indexs),' resnow=',    &
!    &           resnow(i,k),' qcice=',qqcice,' cicep=',cicep(i,k)
!           endif

            endif                                 ! end if_qcice block
          endif                                   ! end if_totcnd block

        enddo
      enddo
!
!...................................
      end subroutine rsipath2
!-----------------------------------

      end MODULE module_microphysics

