.. _CodeOverview:

*************************
Code Overview
*************************

===================================================
UFS Weather Model Hierarchical Repository Structure
===================================================

The UFS Weather Model (:term:`WM`) repository supports the :term:`UFS` short- and medium-range weather applications (:term:`SRW` / :term:`MRW` Apps). The WM repository contains atmosphere, ocean, sea ice, and wave components, as well as some infrastructure components. Each of these subcomponents has its own repository. All the repositories are currently located in GitHub with public access to the broad community. :numref:`Table %s <Repo_Structure>` describes the list of repositories that comprises the UFS WM.

.. _Repo_Structure:

.. list-table:: *List of Repositories that comprise the ufs-weather-model*
  :widths: 50 50
  :header-rows: 1

  * - Repository Description
    - Authoritative repository URL
  * - Umbrella repository for the UFS Weather Model
    - https://github.com/ufs-community/ufs-weather-model
  * - Framework to connect the CCPP library to a host model
    - https://github.com/NCAR/ccpp-framework
  * - CCPP library of physical parameterizations
    - https://github.com/NCAR/ccpp-physics
  * - Umbrella repository for the physics and dynamics of the atmospheric model (FV3) 
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
  * - Air Quality Model (AQM)
    - https://github.com/NOAA-EMC/AQM

In the table, the left column contains a description of each repository, and the right column shows the component repositories which are pointing to (or will point to) the authoritative repositories. The UFS WM currently uses Git submodules to manage the sub-components.
   
===================
Directory Structure
===================

The umbrella repository for the UFS Weather Model is named ``ufs-weather-model``.  Under this repository reside a number of submodules that are nested in specific directories under the parent repository's working directory. When the ``ufs-weather-model`` repository is cloned, the basic directory structure will be similar to the example below. Files and some directories have been removed for brevity. Directories in parentheses will appear only after a submodule update (``git submodule update --init --recursive``). 

.. code-block:: console

   ufs-weather-model/
   ├── build.sh                 -------- script for building the WM
   ├── cmake                    -------- cmake configuration files
   ├── CMakeLists.txt         
   ├── CMakeModules           
   ├── doc                      -------- User Guide files
   ├── driver                 
   ├── FV3                      -------- UFSAtm atmosphere model
   │   ├── (atmos_cubed_sphere) -------- FV3 dynamical core
   │   │   ├── (docs)
   │   │   ├── (driver)
   │   │   ├── (model)
   │   │   └── (tools)
   │   ├── (ccpp)               -------- Common Community Physics Package
   │   │   ├── (config)
   │   │   ├── (driver)
   │   │   ├── (framework)      -------- CCPP framework
   │   │   ├── (physics)        -------- CCPP compliant physics schemes
   │   │   └── (suites)         -------- CCPP physics suite definition files (SDFs)
   │   ├── (cpl)                -------- Coupling field data structures
   │   ├── (io)                 -------- UFSAtm write grid comp code
   │   └── (stochastic_physics) -------- Wrapper for stochastic physics
   │
   ├── stochastic_physics       -------- stochastic physics pattern generator
   ├── AQM
   │     └── (src)
   │         ├── (model)
   │            └── (CMAQ)                         --------- EPA AQ Model
   ├── CICE-interface
   │    └── CICE                                   --------- CICE6 sea ice model
   │        ├── (icepack)                          --------- Sea ice column physics
   │        └── (cicecore/drivers/nuopc/cmeps)     --------- NUOPC CICE6 cap
   ├── GOCART
   │    └── (ESMF)                                 --------- GOCART model
   ├── HYCOM-interface
   │    └── HYCOM                                  --------- HYCOM ocean model
   │        └── (NUOPC)                            --------- NUOPC HYCOM cap
   ├── MOM6-interface
   │    └── MOM6
   │        ├── (src)                              --------- MOM6 ocean model
   │        └── (config_source/drivers/nuopc_cap)  --------- NUOPC MOM6 cap
   ├── WW3
   │    └── (model)                                --------- WW3 model
   │        └── (esmf)                             --------- NUOPC WW3 cap
   ├── CDEPS-interface
   │     └── CDEPS
   │         ├── (datm)                            --------- CDEPS DATM
   │         └── (docn)                            --------- CDEPS DOCN
   ├── CMEPS-interface
   │    └── CMEPS
   │         └── (cesm)                            --------- CMEPS CESM
   ├── modulefiles          -------- system module files for supported HPC systems
   └── tests                -------- regression test infrastructure
       └── parm
       └── tests
       └── fv3_conf   

The physics subdirectory in the ``gfsphysics`` directory  is not used or supported
as part of this release (all physics is available through the :term:`CCPP` using
the repository described in :numref:`Table %s <Repo_Structure>`).

.. COMMENT: Should we omit CMakeLists.txt, CMakeModules, driver (which I added) or add a description?
.. COMMENT: I don't see a "gfsphysics" directory... Can we remove it or be more specific about where it is? There are two CCPP repos in the table referenced above... Framework and Physics.