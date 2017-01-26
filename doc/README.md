NEMSFV3 Build and Test Instructions
======================================

This page documents the build system and test suite for the NEMSFV3
app.  This is a NEMS app that runs the standalone (uncoupled) Finite
Volume Model (FV3). 

<a name="supbuild"></a>Supported Builds and Platforms
-----------------------------------------------------

The NEMSFV3 app currently only works on NOAA WCOSS Cray.
No other platforms are tested or supported at this time.  
There are plans to port to the WCOSS phase1 and phase2, theia machine.

The app can be built in the following ways via the the script 
compile_nems.sh under tests directory. The simple regression test
suite is set up for senity check when the changes are coming in.

