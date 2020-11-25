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
  * - Infrastructure: Flexible Modeling System
    - https://github.com/NOAA-GFDL/FMS
  * - Infrastructure: NOAA Environmental Modeling System
    - https://github.com/NOAA-EMC/NEMS
  * - Infrastructure: Utilities
    - https://github.com/NOAA-EMC/NCEPLIBS-pyprodutil
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

In the table, the left column contains a description of each repository, and the right column shows the component repositories which are pointing to (or will point to) the authoritative repositories. The ufs-weather-model currently uses git submodule to manage the sub-components.

The umbrella repository for the UFS Weather Model is named ufs-weather-model.  Under this repository reside a number of submodules that are nested in specific directories under the parent repository’s working directory.  When the ufs-weather-model repository is cloned, the *.gitmodules* file creates the following directories:

.. code-block:: console

   ufs-weather-model/
   ├── FMS                                     https://github.com/NOAA-GFDL/FMS
   ├── FV3                                     https://github.com/NOAA-EMC/fv3atm
   │   ├── atmos_cubed_sphere                  https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere
   │   ├── ccpp
   │   │   ├── framework                       https://github.com/NCAR/ccpp-framework
   │   │   ├── physics                         https://github.com/NCAR/ccpp-physics
   ├── NEMS                                    https://github.com/NOAA-EMC/NEMS
   │   └── tests/produtil/NCEPLIBS-pyprodutil  https://github.com/NOAA-EMC/NCEPLIBS-pyprodutil
   ├── stochastic_physics                      https://github.com/noaa-psd/stochastic_physics

===================
Directory Structure
===================

When the ufs-weather-model is cloned, the basic directory structure will be similar to the example below. Files and some directories have been removed for brevity.

.. code-block:: console

   ufs-weather-model/
   ├── cmake               --------- cmake configuration files
   ├── compsets            --------- configurations used by some regression tests
   ├── conf                --------- compile options for Tier 1 and 2 platforms
   ├── doc                 --------- READMEs with build, reg-test hints
   ├── FMS                 --------- The Flexible Modeling System (FMS),a software framework
   ├── FV3                 --------- FV3 atmosphere model
   │   ├── atmos_cubed_sphere   ---- FV3 dynamic core
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
   │   ├── gfsphysics
   │   │   ├── CCPP_layer
   │   │   ├── GFS_layer
   │   │   └── physics     --------- unused - IPD version of physics codes
   │   ├── io              --------- FV3 write grid comp code
   │   ├── ipd             --------- unused - IPD driver/interfaces
   |   ├── stochastic_physics  ----- Cmakefile for stochastic physics code
   ├── log                 --------- log files from NEMS compset regression tests
   ├── modulefiles         --------- system module files for supported HPC systems
   ├── NEMS                --------- NOAA Earth Modeling System framework
   │   ├── exe
   │   ├── src
   │   └── test
   ├── parm                --------- regression test configurations
   ├── stochastic_physics   -------- stochastic physics pattern generator
   ├── tests               --------- regression test scripts

The physics subdirectory in the *gfsphysics* directory  is not used or supported
as part of this release (all physics is available through the :term:`CCPP` using
the repository described in :numref:`Table %s <Repo_Structure>`).
