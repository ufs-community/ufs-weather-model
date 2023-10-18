.. _CodeOverview:

*************************
Technical Overview
*************************

.. _SupportedPlatforms:

=====================================================================
Supported Platforms and Compilers for Running the UFS Weather Model
=====================================================================

Four levels of support have been defined for :term:`UFS` applications, and the UFS Weather Model (:term:`WM`) operates under this paradigm: 

* **Level 1** *(Preconfigured)*: Prerequisite software libraries are pre-built and available in a central location; code builds and runs; full testing of model.
* **Level 2** *(Configurable)*: Prerequisite libraries are not available in a centralized location but are expected to install successfully; code builds and runs; full testing of model.
* **Level 3** *(Limited-test platforms)*: Libraries and code build on these systems, but there is limited testing with running the model.
* **Level 4** *(Build-only platforms)*: Libraries and code build, but running the model is not tested.

Level 1 Systems
==================
Preconfigured (Level 1) systems for the UFS WM already have the required external libraries available in a central location via :term:`spack-stack`. The WM is expected to build and run out-of-the-box on these systems, and users can download the WM code without first installing prerequisite software. Additionally, regression test data is already available on these systems. In general, users must have access to these Level 1 systems in order to use them.

Currently, Level 1 (or Tier-1) platforms for regression testing are: 

   * WCOSS2 (Intel)
   * Gaea (Intel)
   * Hera (Intel/GNU compilers)
   * Jet (Intel)
   * Orion (Intel)
   * Hercules (Intel/GNU compilers)
   * AWS Docker container (Intel)

More information is available in the `UFS WM wiki <https://github.com/ufs-community/ufs-weather-model/wiki/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`__. 

Level 2-4 Systems
===================

On non-Level 1 platforms, users must install the required libraries before building the UFS WM. Additionally, users must stage the required data in order to run regression tests. Once the prerequisite libraries are installed, and the data has been staged, the WM should build and run successfully. However, users may need to perform additional troubleshooting on Level 3 or 4 systems since little or no testing is conducted on these systems.

Currently, Level 2 platforms for regression testing are:

   * S4 (Intel)

===================================================
UFS Weather Model Hierarchical Repository Structure
===================================================

The UFS :term:`WM` repository supports the :term:`UFS` short- and medium-range weather applications (:term:`SRW` / :term:`MRW` Apps). The WM repository contains atmosphere, ocean, sea ice, land, and wave components, as well as some infrastructure components. Each of these subcomponents has its own repository. All the repositories are currently located in GitHub with public access to the broader community. :numref:`Table %s <Repo_Structure>` describes the list of repositories that comprise the UFS WM.

.. _Repo_Structure:

.. list-table:: *List of Repositories that comprise the ufs-weather-model*
  :widths: 50 50
  :header-rows: 1

  * - Repository Description
    - Authoritative repository URL
  * - Umbrella repository for the UFS Weather Model
    - https://github.com/ufs-community/ufs-weather-model
  * - Framework to connect the :term:`CCPP` library to a host model
    - https://github.com/NCAR/ccpp-framework
  * - CCPP library of physical parameterizations
    - https://github.com/NCAR/ccpp-physics
  * - Umbrella repository for the physics and dynamics of the atmospheric model (FV3) 
    - https://github.com/NOAA-EMC/fv3atm
  * - :term:`FV3` dynamical core
    - https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere
  * - Stochastic physics pattern generator
    - https://github.com/noaa-psd/stochastic_physics
  * - Modular Ocean Model (:term:`MOM6`)
    - https://github.com/NOAA-EMC/MOM6
  * - HYbrid Coordinate Ocean Model (:term:`HYCOM`)
    - https://github.com/NOAA-EMC/HYCOM-src
  * - Los Alamos sea ice model (:term:`CICE6`)
    - https://github.com/NOAA-EMC/CICE
  * - NOAA/NCEP WAVEWATCH III Model (:term:`WW3`)
    - https://github.com/NOAA-EMC/WW3
  * - The Goddard Chemistry Aerosol Radiation and Transport (:term:`GOCART`)
    - https://github.com/GEOS-ESM/GOCART 
  * - NUOPC Community Mediator for Earth Prediction Systems (:term:`CMEPS`)
    - https://github.com/NOAA-EMC/CMEPS
  * - Community Data Models for Earth Prediction Systems (:term:`CDEPS`)
    - https://github.com/NOAA-EMC/CDEPS
  * - Air Quality Model (:term:`AQM`)
    - https://github.com/NOAA-EMC/AQM
  * - Noah-MP Land Surface Model (Noah-MP)
    - https://github.com/NOAA-EMC/noahmp

In the table, the left-hand column contains a description of each repository, and the 
right-hand column shows the GitHub location of the authoritative component repositories. 
The UFS WM currently uses Git submodules to manage these subcomponents.
   
===================
Directory Structure
===================

The umbrella repository for the UFS WM is named ``ufs-weather-model``. Under this repository reside a number of submodules that are nested in specific directories under the parent repository's working directory. When the ``ufs-weather-model`` repository is cloned, the basic directory structure will be similar to the example below. Files and some directories have been removed for brevity. Directories in parentheses will appear only after a recursive clone or submodule update (``git submodule update --init --recursive``). 

.. code-block:: console

   ufs-weather-model
    ├── AQM
    │     └── (src)
    │         ├── (model)
    │            └── (CMAQ)                      -------- EPA Air Quality Model
    ├── build.sh                                 -------- script for building the WM
    ├── CDEPS-interface
    │     └── CDEPS
    │         ├── (datm)                         -------- CDEPS DATM
    │         └── (docn)                         -------- CDEPS DOCN
    ├── CICE-interface
    │    └── CICE                                -------- CICE6 sea ice model
    │        ├── (icepack)                       -------- Sea ice column physics
    │        └── (cicecore/drivers/nuopc/cmeps)  -------- NUOPC CICE6 cap
    ├── cmake                                    -------- cmake configuration files
    ├── CMakeLists.txt         
    ├── CMakeModules           
    ├── CMEPS-interface
    │    └── CMEPS
    │         └── (cesm)                         -------- CMEPS CESM
    ├── doc                                      -------- User Guide files
    ├── driver                 
    ├── FV3                                      -------- UFSAtm atmosphere model
    │   ├── (atmos_cubed_sphere)                 -------- FV3 dynamical core
    │   │   ├── (docs)
    │   │   ├── (driver)
    │   │   ├── (model)
    │   │   └── (tools)
    │   ├── (ccpp)                               -------- Common Community Physics Package
    │   │   ├── (config)
    │   │   ├── (driver)
    │   │   ├── (framework)                      -------- CCPP framework
    │   │   ├── (physics)                        -------- CCPP-compliant physics schemes
    │   │   └── (suites)                         -------- CCPP physics suite definition files (SDFs)
    │   ├── (cpl)                                -------- Coupling field data structures
    │   ├── (io)                                 -------- UFSAtm write grid comp code
    │   └── (stochastic_physics)                 -------- Wrapper for stochastic physics
    ├── GOCART
    │    └── (ESMF)                              -------- GOCART model
    ├── HYCOM-interface
    │    └── HYCOM                               -------- HYCOM ocean model
    │        └── (NUOPC)                         -------- NUOPC HYCOM cap
    ├── LICENSE.md
    ├── modulefiles                              -------- system module files for supported HPC systems
    ├── MOM6-interface
    │    └── MOM6
    │        ├── (src)                           -------- MOM6 ocean model
    │        └── (config_source/drivers/nuopc_cap)  -------- NUOPC MOM6 cap
    ├── NOAHMP-interface
    │    └── noahmp
    │        ├── (cmake)                         -------- Noah-MP land model
    │        ├── (drivers/nuopc)                 -------- NUOPC Noah-MP cap
    │        ├── (parameters)
    │        └── (src)
    ├── README.md
    ├── stochastic_physics                       -------- stochastic physics pattern generator
    ├── tests                                    -------- regression test infrastructure
    │   └── parm
    │   └── tests
    │   └── fv3_conf   
    └── WW3
         └── (model)                             -------- WW3 model
             └── (src)                           -------- NUOPC WW3 caps