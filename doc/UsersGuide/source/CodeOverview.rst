.. _CodeOverview:

*************************
Code Overview
*************************

===================================================
UFS Weather Model Hierarchical Repository Structure
===================================================

The ufs-weather-model repository supports the short- and medium-range UFS applications. It contains atmosphere and wave components and some infrastructure components. Each of these components has its own repository. All the repositories are currently located in GitHub with public access to the broad community. :numref:`Table %s <Repo_Structure>` describes the list of repositories that comprises the ufs-weather-model.

.. _Repo_Structure:

.. list-table:: *List of Repositories that comprise the ufs-weather-model*
  :widths: 50 50
  :header-rows: 1

  * - Repository Description
    - Authoritative repository URL
  * - Umbrella repository for the UFS Weather Model
    - https://github.com/ufs-community/ufs-weather-model
  * - Infrastructure: NOAA Environmental Modeling System
    - https://github.com/NOAA-EMC/NEMS
  * - Framework to connect the CCPP library to a host model
    - https://github.com/NCAR/ccpp-framework
  * - CCPP library of physical parameterizations
    - https://github.com/NCAR/ccpp-physics
  * - Umbrella repository for the physics and dynamics of the atmospheric model
    - https://github.com/NOAA-EMC/fv3atm
  * - FV3 dynamical core
    - https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere
  * - Stochastic physics pattern generator
    - https://github.com/noaa-psd/stochastic_physics
  * - Modular Ocean Model (MOM6)
    - https://github.com/NOAA-EMC/MOM6
  * - HYbrid Coordinate Ocean Model (HYCOM)
    - https://github.com/NOAA-EMC/HYCOM-src
  * - Los Alamos sea ice model (CICE6)
    - https://github.com/NOAA-EMC/CICE
  * - NOAA/NCEP WAVEWATCH III Model (WW3)
    - https://github.com/NOAA-EMC/WW3
  * - The Goddard Chemistry Aerosol Radiation and Transport (GOCART)
    - https://github.com/GEOS-ESM/GOCART 
  * - NUOPC Community Mediator for Earth Prediction Systems (CMEPS)
    - https://github.com/NOAA-EMC/CMEPS
  * - Community Data Models for Earth Prediction Systems (CDEPS)
    - https://github.com/NOAA-EMC/CDEPS

In the table, the left column contains a description of each repository, and the right column shows the component repositories which are pointing to (or will point to) the authoritative repositories. The ufs-weather-model currently uses git submodule to manage the sub-components.

The umbrella repository for the UFS Weather Model is named ufs-weather-model.  Under this repository reside a number of submodules that are nested in specific directories under the parent repository’s working directory.  When the ufs-weather-model repository is cloned, the *.gitmodules* file creates the following directories:

.. code-block:: console

   ufs-weather-model/
   ├── FV3                                     https://github.com/NOAA-EMC/fv3atm
   │   ├── atmos_cubed_sphere                  https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere
   │   ├── ccpp
   │   │   ├── framework                       https://github.com/NCAR/ccpp-framework
   │   │   ├── physics                         https://github.com/NCAR/ccpp-physics
   ├── NEMS                                    https://github.com/NOAA-EMC/NEMS
   ├── stochastic_physics                      https://github.com/noaa-psd/stochastic_physics
   ├── MOM6-interface
   │    └── MOM6                               https://github.com/NOAA-EMC/MOM6
   ├── HYCOM-interface
   │    └── HYCOM                              https://github.com/NOAA-EMC/HYCOM-src
   ├── CICE-interface
   │    └── CICE                               https://github.com/NOAA-EMC/CICE
   ├── WW3                                     https://github.com/NOAA-EMC/WW3
   ├── GOCART                                  https://github.com/GEOS-ESM/GOCART
   ├── CMEPS-interface
   │    └── CMEPS                              https://github.com/NOAA-EMC/CMEPS
   ├── CDEPS-interface
   │    └── CDEPS                              https://github.com/NOAA-EMC/CDEPS
   
===================
Directory Structure
===================

When the ufs-weather-model is cloned, the basic directory structure will be similar to the example below. Files and some directories have been removed for brevity.


.. code-block:: console

   ufs-weather-model/
   ├── cmake                -------- cmake configuration files
   ├── doc                  -------- User Guide files
   ├── FV3                  -------- UFSAtm atmosphere model
   │   ├── atmos_cubed_sphere ------ FV3 dynamic core
   │   │   ├── docs
   │   │   ├── driver
   │   │   ├── model
   │   │   └── tools
   │   ├── ccpp             -------- Common Community Physics Package
   │   │   ├── config
   │   │   ├── driver
   │   │   ├── framework    -------- CCPP framework
   │   │   ├── physics      -------- CCPP compliant physics schemes
   │   │   └── suites       -------- CCPP physics suite definition files (SDFs)
   │   ├── cpl              -------- Coupling field data structures
   │   ├── io               -------- UFSAtm write grid comp code
   │   └── stochastic_physics ------ Wrapper for stochastic physics
   │
   ├── NEMS                 -------- NOAA Earth Modeling System framework
   ├── stochastic_physics   -------- stochastic physics pattern generator
   ├── MOM6-interface
   │    └── MOM6
   │        ├── src                              --------- MOM6 ocean model
   │        └── config_source/drivers/nuopc_cap  --------- NUOPC MOM6 cap
   ├── HYCOM-interface
   │    └── HYCOM                                --------- HYCOM ocean model
   │        └── NUOPC                            --------- NUOPC HYCOM cap
   ├── CICE-interface
   ├── CICE-interface
   │    └── CICE                                 --------- CICE6 sea ice model
   │        ├── icepack                          --------- Sea ice column physics
   │        └── cicecore/drivers/nuopc/cmeps     --------- NUOPC CICE6 cap
   ├── WW3
   │    └── model                                --------- WW3 model
   │        └── esmf                             --------- NUOPC WW3 cap
   ├── GOCART
   │    └── ESMF                                 --------- GOCART model
   ├── CDEPS-interface
   │     └── CDEPS
   │         ├── datm                            --------- CDEPS DATM
   │         └── docn                            --------- CDEPS DOCN
   ├── modulefiles          -------- system module files for supported HPC systems
   ├── tests                -------- regression test infrastructure
   │   └── parm
   │   └── tests
   │   └── fv3_conf   

The physics subdirectory in the *gfsphysics* directory  is not used or supported
as part of this release (all physics is available through the :term:`CCPP` using
the repository described in :numref:`Table %s <Repo_Structure>`).
