.. _BuildingAndRunning:

******************************************
Building and Running the UFS Weather Model
******************************************

===================================
Supported Platforms & Compilers
===================================
Before running the Weather Model (:term:`WM`), users should determine which of the 
:ref:`levels of support <SupportedPlatforms>` 
is applicable to their system. Generally, Level 1 & 2 systems are restricted to those with access 
through NOAA and its affiliates. These systems are named (e.g., Hera, Orion, Cheyenne). 
Level 3 & 4 systems include certain personal computers or non-NOAA-affiliated HPC systems. 
The prerequisite software libraries for building the WM already exist in a centralized location on Level 1/preconfigured 
systems, so users may skip directly to :ref:`getting the data <GetData>` and downloading the code. 
On other systems, users will need to build the prerequisite libraries using :term:`spack-stack` or :term:`HPC-Stack`. 

=======================
Prerequisite Libraries
=======================

The UFS WM requires a number of libraries.
The WM uses two categories of libraries, which are available as a bundle via 
:term:`spack-stack` or :term:`HPC-Stack`:

   #. :term:`NCEP` libraries (:term:`NCEPLIBS`): These are libraries developed for use with NOAA weather models.
      Most have an NCEPLIBS prefix in the repository (e.g., NCEPLIBS-bacio). Select tools from the UFS
      Utilities repository (:term:`UFS_UTILS`) are also included in this category. 

   #. Third-party libraries (:term:`NCEPLIBS-external`): These are libraries that were developed externally to
      the UFS Weather Model. They are general software packages that are also used by other community models. 
      Building these libraries is optional if users can point to existing builds of these libraries on their system
      instead. 

.. note::
   Currently, spack-stack is the software stack validated by the UFS WM for running 
   :term:`regression tests <RT>`. Spack-stack is a Spack-based method for installing UFS 
   prerequisite software libraries. UFS applications and components are also shifting to 
   spack-stack from HPC-Stack but are at various stages of this transition. 
   Although users can still build and use HPC-Stack, the UFS WM no longer uses HPC-Stack 
   for validation, and support for this option is being deprecated. 

----------------
Common Modules
----------------

As of May 19, 2023, the UFS WM Regression Tests (:term:`RTs <RT>`) on Level 1 systems use the following common modules: 

.. code-block:: console

   bacio/2.4.1
   crtm/2.4.0
   esmf/8.3.0b09
   fms/2022.04
   g2/3.4.5
   g2tmpl/1.10.2
   gftl-shared/v1.5.0
   hdf5/1.10.6
   ip/3.3.3
   jasper/2.0.25
   libpng/1.6.37
   mapl/2.22.0-esmf-8.3.0b09
   netcdf/4.7.4
   pio/2.5.7
   sp/2.3.3
   w3emc/2.9.2
   zlib/1.2.11

The most updated list of common modules can be viewed in ``ufs_common.lua`` 
`here <https://github.com/ufs-community/ufs-weather-model/blob/develop/modulefiles/ufs_common.lua>`__.

.. attention::
   Documentation is available for installing `spack-stack <https://spack-stack.readthedocs.io/en/latest/>`__
   and `HPC-Stack <https://hpc-stack.readthedocs.io/en/latest/>`__, respectively. 
   One of these software stacks (or the libraries they contain) must be installed before running the UFS Weather Model. 

.. _GetData:

============
Get Data
============

The WM RTs require input files to run. 
These include static datasets, files that depend on grid resolution and 
initial/boundary conditions, and model configuration files. On Level 1 and 2 systems, 
the data required to run the WM RTs are already available in the following locations: 

.. _DataLocations:
.. table:: Data Locations for Level 1 & 2 Systems

   +--------------+-----------------------------------------------------+
   | Machine      | File location                                       |
   +==============+=====================================================+
   | Cheyenne     | /glade/scratch/epicufsrt/GMTB/ufs-weather-model/RT  |
   +--------------+-----------------------------------------------------+
   | Gaea         | /lustre/f2/pdata/ncep_shared/emc.nemspara/RT        |
   +--------------+-----------------------------------------------------+
   | Hera         | /scratch1/NCEPDEV/nems/emc.nemspara/RT              |
   +--------------+-----------------------------------------------------+
   | Jet          | /mnt/lfs4/HFIP/hfv3gfs/role.epic/RT                 |
   +--------------+-----------------------------------------------------+
   | Orion        | /work/noaa/nems/emc.nemspara/RT                     |
   +--------------+-----------------------------------------------------+
   | S4           | /data/prod/emc.nemspara/RT                          |
   +--------------+-----------------------------------------------------+ 
   | WCOSS2       | /lfs/h2/emc/nems/noscrub/emc.nems/RT                |
   +--------------+-----------------------------------------------------+ 

For Level 3-4 systems, the data must be added to the user's system. 
Publicly available RT data is available in the `UFS WM Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`__. 
Data for running RTs off of the develop branch is available for the most recent 60 days. 
To view the data, users can visit https://noaa-ufs-regtests-pds.s3.amazonaws.com/index.html. 
To download data, users must select the data they want from the bucket and either download it in their browser or via a ``wget`` command. 
For example, to get the data for ``control_p8`` (specifically the May 17, 2023 ``develop`` branch version of the WM), run: 

.. code-block:: console

   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/atmf000.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/atmf021.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/atmf024.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/GFSFLX.GrbF00
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/GFSFLX.GrbF21
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/GFSFLX.GrbF24
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/GFSPRS.GrbF00
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/GFSPRS.GrbF21
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/GFSPRS.GrbF24
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/sfcf000.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/sfcf021.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/sfcf024.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.coupler.res
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_core.res.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_core.res.tile1.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_core.res.tile2.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_core.res.tile3.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_core.res.tile4.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_core.res.tile5.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_core.res.tile6.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_srf_wnd.res.tile1.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_srf_wnd.res.tile2.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_srf_wnd.res.tile3.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_srf_wnd.res.tile4.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_srf_wnd.res.tile5.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_srf_wnd.res.tile6.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_tracer.res.tile1.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_tracer.res.tile2.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_tracer.res.tile3.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_tracer.res.tile4.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_tracer.res.tile5.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.fv_tracer.res.tile6.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.phy_data.tile1.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.phy_data.tile2.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.phy_data.tile3.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.phy_data.tile4.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.phy_data.tile5.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.phy_data.tile6.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.sfc_data.tile1.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.sfc_data.tile2.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.sfc_data.tile3.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.sfc_data.tile4.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.sfc_data.tile5.nc
   wget https://noaa-ufs-regtests-pds.s3.amazonaws.com/develop-20230517/INTEL/control_p8/RESTART/20210323.060000.sfc_data.tile6.nc

Detailed information on input files can be found in :numref:`Chapter %s <InputsOutputs>`. 

.. _DownloadingWMCode:

==================================
Downloading the Weather Model Code
==================================

To clone the develop branch of the ``ufs-weather-model`` repository and update its submodules, execute the following commands:

.. code-block:: console

  git clone --recursive https://github.com/ufs-community/ufs-weather-model.git ufs-weather-model
  cd ufs-weather-model

Compiling the model will take place within the ``ufs-weather-model`` directory created by this command.

==========================
Building the Weather Model
==========================

----------------------------
Loading the Required Modules
----------------------------

The process for loading modules is fairly straightforward on NOAA :ref:`Level 1 Systems <SupportedPlatforms>`. 
Users may need to make adjustments when running on other systems. 


On NOAA Level 1 & 2 Systems
-----------------------------

Modulefiles for :ref:`preconfigured platforms <SupportedPlatforms>` are located in 
``modulefiles/ufs_<platform>.<compiler>``. For example, to load the modules from the 
``ufs-weather-model`` directory on Hera:

.. code-block:: console

    module use modulefiles
    module load ufs_hera.intel

Note that loading this module file will also set the CMake environment variables shown in
:numref:`Table %s <CMakeEnv>`.

.. _CMakeEnv:

.. table:: *CMake environment variables required to configure the build for the Weather Model*

   +-------------------------+----------------------------------------------+----------------------+
   | **EnvironmentVariable** | **Description**                              | **Hera Intel Value** |
   +=========================+==============================================+======================+
   |  CMAKE_C_COMPILER       | Name of C compiler                           | mpiicc               |
   +-------------------------+----------------------------------------------+----------------------+
   |  CMAKE_CXX_COMPILER     | Name of C++ compiler                         | mpiicpc              |
   +-------------------------+----------------------------------------------+----------------------+
   |  CMAKE_Fortran_COMPILER | Name of Fortran compiler                     | mpiifort             |
   +-------------------------+----------------------------------------------+----------------------+
   |  CMAKE_Platform         | String containing platform and compiler name | hera.intel           |
   +-------------------------+----------------------------------------------+----------------------+

On Other Systems
-------------------

If you are not running on one of the pre-configured platforms, you will need to set the environment variables
manually. For example, in a bash shell, a command in the following form will set the C compiler environment variable:

.. code-block:: console

   export CMAKE_C_COMPILER=</path/to/C/compiler>

.. COMMENT: Update after Zach's PR is merged. 

------------------------------------------------------------------------
Setting the ``CMAKE_FLAGS`` and ``CCPP_SUITES`` Environment Variables
------------------------------------------------------------------------

The UFS Weather Model can be built in one of several configurations (see :numref:`Table %s <UFS-configurations>` for common options). 
The ``CMAKE_FLAGS`` environment variable specifies which configuration to build using the ``-DAPP`` and ``-DCCPP_SUITES`` variables.
Users set which components to build using ``-DAPP``. Users select the :term:`CCPP` suite(s) by setting the 
``CCPP_SUITES`` environment variable at build time in order to have one or more CCPP physics suites available at runtime. 
Multiple suites can be set. Additional variables, such as ``-D32BIT=ON``, 
can be set if the user chooses. These options are documented in :numref:`Section %s <other-build-options>`. 
The following examples assume a bash shell.

ATM Configurations
---------------------

.. _atm:

**Standalone ATM**

For the ``ufs-weather-model ATM`` configuration (standalone :term:`ATM`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16"

.. _atmw:

**ATMW**

For the ``ufs-weather-model ATMW`` configuration (standalone ATM coupled to :term:`WW3`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATMW -DCCPP_SUITES=FV3_GFS_v16"

.. _atmaero:

**ATMAERO**

For the ``ufs-weather-model ATMAERO`` configuration (standalone ATM coupled to :term:`GOCART`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATMAERO -DCCPP_SUITES=FV3_GFS_v17_p8"

.. _atmaq:

**ATMAQ**

For the ``ufs-weather-model ATMAQ`` configuration (standalone ATM coupled to :term:`CMAQ`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATMAQ -DCCPP_SUITES=FV3_GFS_v15p2"

.. _atml:

**ATML**

For the ``ufs-weather-model ATML`` configuration (standalone ATM coupled to :term:`LND`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATML -DCCPP_SUITES=FV3_GFS_v17_p8"

S2S Configurations 
----------------------

.. _s2s:

**S2S**

For the ``ufs-weather-model S2S`` configuration (coupled atm/ice/ocean):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2S -DCCPP_SUITES=FV3_GFS_v17_coupled_p8"

To turn on debugging flags, add ``-DDEBUG=ON`` flag after ``-DAPP=S2S``. Users can allow verbose build messages by running: 

.. code-block:: console

    export BUILD_VERBOSE=1

To receive atmosphere-ocean fluxes from the CMEPS :term:`mediator`, add the argument ``-DCMEPS_AOFLUX=ON``.
For example:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2S -DCCPP_SUITES=FV3_GFS_v17_coupled_p8_sfcocn -DCMEPS_AOFLUX=ON"

.. _s2sa:

**S2SA**

For the ``ufs-weather-model S2SA`` configuration (atm/ice/ocean/aerosols):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2SA -DCCPP_SUITES=FV3_GFS_2017_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled,FV3_GFS_v16_coupled_noahmp"

..
   CHECK: DAPP flag and physics suites

.. _s2sw:

**S2SW**

For the ``ufs-weather-model S2SW`` configuration (atm/ice/ocean/wave):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2SW -DCCPP_SUITES=FV3_GFS_v17_coupled_p8"

.. _s2swa:

**S2SWA**

For the ``ufs-weather-model S2SWA`` configuration (atm/ice/ocean/wave/aerosols):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2SWA -DCCPP_SUITES=FV3_GFS_v17_coupled_p8,FV3_GFS_cpld_rasmgshocnsstnoahmp_ugwp"

.. _ng-godas:

NG-GODAS Configuration
------------------------

For the ``ufs-weather-model NG-GODAS`` configuration (atm/ocean/ice/data assimilation): 

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=NG-GODAS"

.. COMMENT: Check! --> In rt.conf, no CCPP suite is set. Is there a default one?

HAFS Configurations
----------------------

.. _hafs:

**HAFS**

For the ``ufs-weather-model HAFS`` configuration (atm/ocean) in 32 bit:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFS -D32BIT=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf_nonsst,FV3_HAFS_v0_gfdlmp_tedmf"

.. _hafsw:

**HAFSW**

For the ``ufs-weather-model HAFSW`` configuration (atm/ocean/wave) in 32-bit with moving nest:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFSW -D32BIT=ON -DMOVING_NEST=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf,FV3_HAFS_v0_gfdlmp_tedmf_nonsst,FV3_HAFS_v0_thompson_tedmf_gfdlsf"

.. _hafs-all:

**HAFS-ALL**

For the ``ufs-weather-model HAFS-ALL`` configuration (data/atm/ocean/wave) in 32 bit:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFS-ALL -D32BIT=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf,FV3_HAFS_v0_gfdlmp_tedmf_nonsst"

LND Configuration
----------------------

.. _lnd:

**LND**

For the ``ufs-weather-model LND`` configuration (datm/land):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=LND"

------------------
Building the Model
------------------

.. COMMENT: Is the "Building the Model" section necessary? Can users just run the RT without?

The UFS Weather Model uses the CMake build system. There is a build script called ``build.sh`` in the
top-level directory of the WM repository that configures the build environment and runs the ``make``
command. This script also checks that all necessary environment variables have been set.

If any of the environment variables have not been set, the ``build.sh`` script will exit with a message similar to:

.. code-block:: console

   ./build.sh: line 11: CMAKE_Platform: Please set the CMAKE_Platform environment variable, e.g. [macosx.gnu|linux.gnu|linux.intel|hera.intel|...]

The WM can be built by running the following command from the ``ufs-weather-model`` directory:

.. code-block:: console

   ./build.sh

Once ``build.sh`` is finished, users should see the executable, named ``ufs_model``, in the ``ufs-weather-model/build/`` directory.
If users prefer to build in a different directory, specify the ``BUILD_DIR`` environment variable. For example: ``export BUILD_DIR=test_cpld``
will build in the ``ufs-weather-model/test_cpld`` directory instead.

Expert help is available through `GitHub Discussions <https://github.com/ufs-community/ufs-weather-model/discussions/categories/q-a>`__. Users may post questions there for help with difficulties related to the UFS WM.

.. _run-wm:

=================
Running the Model
=================

.. attention::
   Although the following discussions are general, users may not be able to execute the script successfully "as is" unless they are on a 
   `Tier-1 platform <https://github.com/ufs-community/ ufs-weather-model/wiki/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`__.

.. _UsingRegressionTest:

--------------------------------
Using the Regression Test Script
--------------------------------

Users can run a number of preconfigured regression test cases from the ``rt.conf`` file 
using the regression test script ``rt.sh`` in the ``tests`` directory. 
``rt.sh`` is the top-level script that calls lower-level scripts to build specified 
WM configurations, set up environments, and run tests. 
Users must edit the ``rt.conf`` file to indicate which tests/configurations to run. 

.. _rt.conf:

The ``rt.conf`` File
------------------------

Each line in the PSV (Pipe-separated values) file, ``rt.conf``, contains four columns of information. 
The first column specifies whether to build a test (``COMPILE``) or run a test (``RUN``). 
The second column specifies either configuration information for building a test or 
the name of a test to run.
Thus, the second column in a ``COMPILE`` line will list the application to build (e.g., ``-DAPP=S2S``), 
the CCPP suite to use (e.g., ``-DCCPP_SUITES=FV3_GFS_2017_coupled``), and additional build options 
(e.g., ``-DDEBUG=ON``) as needed. On a ``RUN`` line, the second column will contain a test name 
(e.g., ``control_p8``). The test name should match the name of one of the test files in the 
``tests/tests`` directory or, if the user is adding a new test, the name of the new test file. 
The third column of ``rt.conf`` relates to the platform; 
if blank, the test can run on any WM Tier-1 platform. 
The fourth column deals with baseline creation 
(see information on ``-c`` option :ref:`below <cmd-line-opts>` for more), 
and ``fv3`` means that the test will be included during baseline creation.

The order of lines in ``rt.conf`` matters
since ``rt.sh`` processes them sequentially; a ``RUN`` line should be preceeded
by a ``COMPILE`` line that builds the model used in the test. The following
``rt.conf`` file excerpt builds the standalone ATM model with GFS_v16 physics 
in 32-bit mode and then runs the ``control`` test:

.. code-block:: console

    COMPILE | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16 -D32BIT=ON | | fv3
    RUN     | control                                        | | fv3

The ``rt.conf`` file includes a large number of tests. If the user wants to run
only specific tests, s/he can either (1) comment out the tests to be skipped (using the ``#`` prefix)
or (2) create a new file (e.g., ``my_rt.conf``), add the tests, and execute ``./rt.sh -l my_rt.conf``.

On NOAA RDHPCS
------------------

On `Tier-1 platforms <https://github.com/ufs-community/ufs-weather-model/wiki
/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`__, users can run 
regression tests by editing the ``rt.conf`` file and executing:

.. code-block:: console

    ./rt.sh -l rt.conf

Users may need to add additional command line arguments or change information in the ``rt.sh`` file as well. 
This information is provided in :numref:`Section %s <rt.sh>` below. 

On Other Systems
------------------

Users on non-NOAA systems will need to make adjustments to several files in the 
``tests`` directory before running ``rt.sh``, including:
  
   * ``rt.sh``
   * ``run_test.sh``
   * ``detect_machine.sh``
   * ``default_vars.sh``
   * ``fv3_conf/fv3_slurm.IN_*``
   * ``fv3_conf/compile_slurm.IN_*``
   * ``compile.sh``
   * ``module-setup.sh``

.. _rt.sh:

The ``rt.sh`` File
---------------------

This section contains additional information on command line options and troubleshooting for the ``rt.sh`` file. 

.. _cmd-line-opts:

Optional Arguments
^^^^^^^^^^^^^^^^^^^^^

To display detailed information on how to use ``rt.sh``, users can simply run ``./rt.sh``, which will output the following options: 

.. code-block:: console

   ./rt.sh -c | -e | -h | -k | -w | -d | -l <file> | -m | -n <name> | -r 
      -c  create new baseline results
      -e  use ecFlow workflow manager
      -h  display this help 
      -k  keep run directory after rt.sh is completed
      -l  runs test specified in <file>
      -m  compare against new baseline results
      -n  run single test <name>
      -r  use Rocoto workflow manager
      -w  for weekly_test, skip comparing baseline results
      -d  delete run direcotries that are not used by other tests

.. COMMENT: An -n option is discussed below. Why is this not printed when running ./rt.sh? 

When running a large number (10's or 100's) of tests, the ``-e`` or ``-r`` options can significantly
decrease testing time by using a workflow manager (ecFlow or Rocoto, respectively) to queue the jobs 
according to dependencies and run them concurrently. 
The ``-n`` option can be used to run a single test; for example, ``./rt.sh -n control`` 
will build the ATM model and run the ``control`` test. 
The ``-c`` option is used to create a baseline. New baselines are needed when code changes lead 
to result changes and therefore deviate from existing baselines on a bit-for-bit basis.

To run ``rt.sh`` using a custom configuration file and the Rocoto workflow manager, 
create the configuration file (e.g. ``my_rt.conf``) based on the desired tests in 
``rt.conf``, and run:

.. code-block:: console

   ./rt.sh -r -l my_rt.conf

adding additional arguments as desired. 

To run a single test, users can try the following command instead of creating a ``my_rt.conf`` file:

.. code-block:: console

   ./rt.sh -r -k -n control_p8

Troubleshooting
^^^^^^^^^^^^^^^^^^

Users may need to adjust certain information in the ``rt.sh`` file, such as 
the *Machine* and *Account* variables (``$MACHINE_ID`` and ``$ACCNR``), for the tests to run 
correctly. If there is a problem with these or other variables (e.g., file paths), the output should indicate where: 

.. code-block:: console
   :emphasize-lines: 5,6

   + echo 'Machine: ' hera.intel '    Account: ' nems
   Machine:  hera.intel     Account:  nems
   + mkdir -p /scratch1/NCEPDEV/stmp4/First.Last
   mkdir: cannot create directory ‘/scratch1/NCEPDEV/stmp4/First.Last’: Permission denied
   ++ echo 'rt.sh error on line 370'
   rt.sh error on line 370

Then, users can adjust the information in ``rt.sh`` accordingly. 

.. _log-files:

Log Files
------------

The regression test generates a number of log files. The summary log file
``RegressionTests_<machine>.<compiler>.log`` in the ``tests`` directory compares
the results of the test against the baseline for a given platform and
reports the outcome: 

   * ``'Missing file'`` results when the expected files from the simulation are not found and typically occurs when the simulation did not run to completion; 
   * ``'OK'`` means that the simulation results are bit-for-bit identical to those of the baseline; 
   * ``'NOT OK'`` when the results are **not** bit-for-bit identical; and 
   * ``'Missing baseline'`` when there is no baseline data to compare against.

More detailed log files are located in the ``tests/log_<machine>.<compiler>/`` directory.
The run directory path, which corresponds to the value of ``RUNDIR`` in the ``run_<test-name>`` file, 
is particularly useful. ``$RUNDIR`` is a self-contained (i.e., sandboxed) 
directory with the executable file, initial conditions, model configuration files, 
environment setup scripts and a batch job submission script. The user can run the test 
by navigating into ``$RUNDIR`` and invoking the command:

.. code-block:: console

    sbatch job_card

This can be particularly useful for debugging and testing code changes. Note that
``$RUNDIR`` is automatically deleted at the end of a successful regression test;
specifying the ``-k`` option retains the ``$RUNDIR``, e.g. ``./rt.sh -l rt.conf -k``.

Inside the ``$RUNDIR`` directory are a number of model configuration files (``input.nml``, 
``model_configure``, ``nems.configure``) and other application
dependent files (e.g., ``ice_in`` for the Subseasonal-to-Seasonal Application).
These model configuration files are
generated by ``rt.sh`` from the template files in the ``tests/parm`` directory.
Specific values used to fill in the template files are test-dependent and
are set in two stages. First, default values are specified in ``tests/default_vars.sh``, and
the default values are overriden if necessary by values specified in a test file
``tests/tests/<test-name>``. For example, the variable ``DT_ATMOS`` is initially assigned 1800 
in the function ``export_fv3`` of the script ``default_vars.sh``, but the test file 
``tests/tests/control`` overrides this setting by reassigning 720 to the variable.

The files ``fv3_run`` and ``job_card`` also reside in the ``$RUNDIR`` directory. 
These files are generated from the template files in the ``tests/fv3_conf``
directory. ``job_card`` is a platform-specific batch job submission script, while 
``fv3_run`` prepares the initial conditions for the test by copying relevant data from the
input data directory of a given platform to the ``$RUNDIR`` directory.
:numref:`Table %s <RTSubDirs>` summarizes the subdirectories discussed above.

.. _RTSubDirs:

.. table:: *Regression Test Subdirectories*

   +-----------------+--------------------------------------------------------------------------------------+
   | **Name**        | **Description**                                                                      |
   +=================+======================================================================================+
   | tests/          | Regression test root directory. Contains rt-related scripts and the summary log file |
   +-----------------+--------------------------------------------------------------------------------------+
   | tests/tests/    | Contains specific test files                                                         |
   +-----------------+--------------------------------------------------------------------------------------+
   | tests/parm/     | Contains templates for model configuration files                                     |
   +-----------------+--------------------------------------------------------------------------------------+
   | tests/fv3_conf/ | Contains templates for setting up initial conditions and a batch job                 |
   +-----------------+--------------------------------------------------------------------------------------+
   | tests/log_*/    | Contains fine-grained log files                                                      |
   +-----------------+--------------------------------------------------------------------------------------+


.. _new-test:

Creating a New Test
----------------------

When a developer needs to create a new test for his/her implementation, the
first step would be to identify a test in the ``tests/tests`` directory that can
be used as a basis and to examine the variables defined in the test file. As
mentioned above, some of the variables may be overrides for those defined in
``default_vars.sh``. Others may be new variables that are needed specifically
for that test. Default variables and their values are defined in the ``export_fv3``
function of the ``default_vars.sh`` script for ATM configurations, the ``export_cpl``
function for S2S configurations, and the ``export_datm`` function for the NG-GODAS configuration.
Also, the names of template files for model configuration and initial conditions
can be identified via variables ``INPUT_NML``, ``NEMS_CONFIGURE`` and ``FV3_RUN`` 
by running ``grep -n INPUT_NML *`` inside the ``tests`` and ``tests/tests`` directories.

.. COMMENT: Is NEMS_CONFIGURE still in there?

.. _UsingOpnReqTest:

---------------------------------------------
Using the Operational Requirement Test Script
---------------------------------------------
The operational requirement test script ``opnReqTest`` in the ``tests`` directory can be used to run
tests in place of ``rt.sh``. Given the name of a test, ``opnReqTest`` carries out a suite of test cases.
Each test case addresses an aspect of the requirements that new operational implementations
must satisfy. These requirements are shown in :numref:`Table %s <OperationalRequirement>`.
For the following discussions on opnReqTest, the user should note the distinction between
``'test name'`` and ``'test case'``. Examples of test names are ``control``, ``cpld_control``
and ``regional_control`` which are all found in the ``tests/tests`` directory, whereas
test case refers to any one of the operational requirements: ``thr``, ``mpi``, ``dcp``, ``rst``, ``bit`` and ``dbg``.

.. _OperationalRequirement:

.. table:: *Operational Requirements*

  +----------+-------------------------------------------------------------------------------+
  | **Case** | **Description**                                                               |
  +==========+===============================================================================+
  | thr      | Varying the number of threads produces the same results                       |
  +----------+-------------------------------------------------------------------------------+
  | mpi      | Varying the number of MPI tasks produces the same results                     |
  +----------+-------------------------------------------------------------------------------+
  | dcp      | Varying the decomposition (i.e. tile layout of FV3) produces the same results |
  +----------+-------------------------------------------------------------------------------+
  | rst      | Restarting produces the same results                                          |
  +----------+-------------------------------------------------------------------------------+
  | bit      | Model can be compiled in double/single precision and run to completion        |
  +----------+-------------------------------------------------------------------------------+
  | dbg      | Model can be compiled and run to completion in debug mode                     |
  +----------+-------------------------------------------------------------------------------+

The operational requirement testing uses the same testing framework as the regression
tests, so it is recommened that the user first read :numref:`Section %s <UsingRegressionTest>`. 
All the files in the subdirectories shown in :numref:`Table %s <RTSubDirs>` are relevant to the
operational requirement test. The only difference is that the ``opnReqTest`` script replaces ``rt.sh``.
The ``tests/opnReqTests`` directory contains
opnReqTest-specific lower-level scripts used to set up run configurations.

On `Tier-1 platforms <https://github.com/ufs-community/ ufs-weather-model/wiki
/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`_, tests can
be run by invoking

.. code-block:: console

    ./opnReqTest -n <test-name>

For example, ``./opnReqTest -n control`` performs all six test cases
listed in :numref:`Table %s <OperationalRequirement>` for the ``control``
test. At the end of the run, a log file ``OpnReqTests_<machine>.<compiler>.log``
is generated in the ``tests`` directory, which informs the user whether each test case
passed or failed. The user can choose to run a specific test case by invoking

.. code-block:: console

    ./opnReqTest -n <test-name> -c <test-case>

where ``<test-case>`` is one or
more comma-separated values selected from ``thr``, ``mpi``, ``dcp``, ``rst``,
``bit``, ``dbg``. For example, ``./opnReqTest -n control -c thr,rst`` runs the
``control`` test and checks the reproducibility of threading and restart.


The user can see different command line options available to ``opnReqTest`` by
executing ``./opnReqTest -h``, which produces the following results:

.. code-block:: console
 
   Usage: opnReqTest -n <test-name> [ -c <test-case> ] [-b] [-d] [-e] [-k] [-h] [-x] [-z]

      -n  specify <test-name>

      -c  specify <test-case>
            defaults to all test-cases: thr,mpi,dcp,rst,bit,dbg,fhz
            comma-separated list of any combination of std,thr,mpi,dcp,rst,bit,dbg,fhz
            
      -b  test reproducibility for bit; compare against baseline
      -d  test reproducibility for dbg; compare against baseline
      -s  test reproducibility for std; compare against baseline
      -e  use ecFlow workflow manager
      -k  keep run directory
      -h  display this help and exit
      -x  skip compile
      -z  skip run

Frequently used options are ``-e`` to use the ecFlow
workflow manager, and ``-k`` to keep the ``$RUNDIR``. Not that the Rocoto workflow manager 
is not used operationally and therefore is not an option. 

As discussed in :numref:`Section %s <log-files>`, the variables and
values used to configure model parameters and to set up initial conditions in the
``$RUNDIR`` directory are set up in two stages. First, ``tests/default_vars.sh``
define default values; then a specific test file in the ``tests/tests`` subdirectory
either overrides the default values or creates new variables if required by the test.
The regression test treats the different test cases shown in
:numref:`Table %s <OperationalRequirement>` as different tests. Therefore, each
test case requires a test file in the ``tests/tests`` subdirectory. Examples include
``control_2threads``, ``control_decomp``, ``control_restart`` and ``control_debug``,
which are just variations of the ``control`` test to check various reproducibilities.
There are two potential issues with this approach. First, if several different
variations of a given test were created and included in the ``rt.conf`` file,
there would be too many tests to run. Second, if a new test is added by the user, s/he
will also have to create these variations. The idea behind the operational requirement test is to
automatically configure and run these variations, or test cases, given a test file.
For example, ``./opnReqTest -n control`` will run all six test cases in
:numref:`Table %s <OperationalRequirement>` based on a single ``control`` test file.
Similarly, if the user adds a new test ``new_test``, then ``./opnReqTest -n new_test`` will
run all test cases. This is done by the operational requirement test script ``opnReqTest`` by adding a third
stage of variable overrides. The related scripts can be found in the ``tests/opnReqTests``
directory.
