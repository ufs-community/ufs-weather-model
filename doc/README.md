NEMSFV3GFS Build and Test Instructions
======================================

This page documents the build system and test suite for the NEMSFV3
app.  This is a NEMS app that runs the standalone (uncoupled) Finite
Volume Model (FV3). 

<a name="supbuild"></a>Supported Builds and Platforms
-----------------------------------------------------

The NEMSFV3 app currently works on NOAA WCOSS Cray and theia.
No other platforms are tested or supported at this time.  
There are plans to port to the WCOSS phase1 and phase2, theia machine.

The app can be built in the following ways via the the script 
compile.sh under tests directory. The simple regression test
suite is set up for senity check when the changes are coming in.

<a name="nemsfv3gfs"></a>NEMS top level driver and fv3_cap
------------------------------------------------------------
NEMSFV3GFS is using NEMS as top level driver, the fv3_cap connects NEMS and 
main fv3 code. The basic calling sequence is listed below:

1)The main program for NEMSfv3gfs is: trunk/NEMS/src/MAIN_NEMS.F90, 
which sets up nems grid component, reads date/time in configure file and 
sets up central clock for nems grid component.

2)The initialiation, run and finalzation of NEMS grid component are 
in trunk/NEMS/src/module_NEMS_GRID_COMP.F90, it calls earth component.

3)earth component  (trunk/NEMS/src/module_EARTH_GRID_COMP.F90) has 
generic setting for the subgrid component defined by ifdef

#ifdef FRONT_FV3
      use FRONT_FV3,        only: FV3_SS   => SetServices
#endif

, and read in nems.configure to decide the run sequence. 
For FV3, the SetServices is defined in fv3_cap.F90.

4)fv3_cap connects nems with fv3. The SetService defines fv3 initialization, 
run and finalize phases. 
4a)In fv3_cap, the initialization also has several phases,
which is required for coupling with other component. 
The fv3 model set up is in one of the initialization phases called InitializeAdvertise, 
In this subroutine, fv3 subroutines are called to initialize fms, mpp, and diag_manager, etc. 
4b)fv3 run phase is in subroutine ModelAdvance. The forecast integration loop is in ModelAdvance,
fv3 subroutines are called to run fv3 dyn core (update_atmos_model_dynamics), physics
(update_atmos_radiation_physics) etc. 
4c)fv3 finalization phase is in subroutine atmos_model_finalize. This subroutine
writes the time stamp for restart and finalizes the mpp, diag_manager and fms.

GFDL is preparing documentation for fv3 code.  For documentation besides NEMS and fv3 cap, 
please contact GFDL. 


