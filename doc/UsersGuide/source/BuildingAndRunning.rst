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
through NOAA and its affiliates. These systems are named (e.g., Hera, Orion, Derecho). 
Level 3 & 4 systems include certain personal computers or non-NOAA-affiliated HPC systems. 
The prerequisite software libraries for building the WM already exist in a centralized location on Level 1/preconfigured 
systems, so users may skip directly to :ref:`getting the data <GetData>` and downloading the code. 
On other systems, users will need to build the prerequisite libraries using :term:`spack-stack`.

=======================
Prerequisite Libraries
=======================

The UFS WM requires a number of libraries.
The WM uses two categories of libraries, which are available as a bundle via 
:term:`spack-stack`:

   #. :term:`NCEP` libraries (:term:`NCEPLIBS`): These are libraries developed for use with NOAA weather models.
      Most have an NCEPLIBS prefix in the repository (e.g., NCEPLIBS-bacio). Select tools from the UFS
      Utilities repository (:term:`UFS_UTILS`) are also included in this category. 

   #. Third-party libraries (:term:`NCEPLIBS-external`): These are libraries that were developed externally to
      the UFS Weather Model. They are general software packages that are also used by other community models. 
      Building these libraries is optional if users can point to existing builds of these libraries on their system
      instead. 

----------------
Common Modules
----------------

As of February 24, 2025, the UFS WM Regression Tests (:term:`RTs <RT>`) on Level 1 systems use the following common modules: 

.. code-block:: console

   bacio/2.4.1
   crtm/2.4.0
   esmf/8.6.0
   fms/2024.01
   g2/3.5.1
   g2tmpl/1.13.0
   gftl-shared/1.6.1
   hdf5/1.14.0
   ip/4.3.0
   jasper/2.0.32
   libpng/1.6.37
   mapl/2.40.3-esmf-8.6.0
   netcdf-c/4.9.2
   netcdf-fortran/4.6.1
   parallelio/2.5.10
   scotch/7.0.4
   sp/2.5.0
   w3emc/2.10.0
   zlib/1.2.13

The most updated list of common modules can be viewed in ``ufs_common.lua`` 
:wm-repo:`here <blob/develop/modulefiles/ufs_common.lua>`.

.. attention::
   Documentation is available for installing `spack-stack <https://spack-stack.readthedocs.io/en/latest/>`_. 
   Spack-stack (or the libraries it contains) must be installed before running the UFS Weather Model. 

.. _GetData:

============
Get Data
============

The WM RTs require input files to run. 
These include static datasets, files that depend on grid resolution and 
initial/boundary conditions, and model configuration files. On Level 1 and 2 systems, 
the data required to run the WM RTs are already available at the following ``DISKNM`` locations: 

.. _DataLocations:

.. list-table:: Data Locations (``$DISKNM``) for Level 1 & 2 Systems
   :widths: 20 50
   :header-rows: 1

   * - Machine
     - File location
   * - Derecho
     - /glade/derecho/scratch/epicufsrt/ufs-weather-model/RT/
   * - Gaea-C6
     - /gpfs/f6/bil-fire8/world-shared/role.epic/UFS-WM_RT
   * - Hera
     - /scratch2/NAGAPE/epic/UFS-WM_RT
   * - Hercules
     - /work/noaa/epic/hercules/UFS-WM_RT
   * - NOAA Cloud (Level 2)
     - /contrib/ufs-weather-model/RT
   * - Orion
     - /work/noaa/epic/UFS-WM_RT
   * - S4 (Level 2)
     - /data/prod/emc.nemspara/RT
   * - WCOSS2
     - /lfs/h2/emc/nems/noscrub/emc.nems/RT

Within ``DISKNM``, input data for the UFS WM is located at the following locations: 

  * **INPUTDATA_ROOT**: ``${DISKNM}/NEMSfv3gfs/input-data-20240501``
  * **INPUTDATA_ROOT_WW3** ``${INPUTDATA_ROOT}/WW3_input_data_20250212``
  * **INPUTDATA_ROOT_BMIC**: ``${DISKNM}/NEMSfv3gfs/BM_IC-20220207``
  * **INPUTDATA_LM4**: ``${INPUTDATA_ROOT}/LM4_input_data``

For Level 3-4 systems, the data must be added to the user's system. 
Publicly available data is available in the `UFS WM Data Bucket <https://registry.opendata.aws/noaa-ufs-regtests/>`_. 
Baseline data for the ``develop`` branch is available for the most recent 60 days. 
The regression testing script (``rt.sh``) has certain default data directories (i.e., ``INPUTDATA_*``) that users may need to change when working on Level 3-4 systems. 
The corresponding data is publicly available in the data bucket. To view the data, users can visit https://noaa-ufs-regtests-pds.s3.amazonaws.com/index.html. 
Users can download the data and update the ``rt.sh`` script to point to the appropriate locations in order to run RTs on their own system: 
  
* ``INPUTDATA_ROOT``: https://noaa-ufs-regtests-pds.s3.amazonaws.com/index.html#input-data-20240501/
* ``INPUTDATA_ROOT_WW3`` https://noaa-ufs-regtests-pds.s3.amazonaws.com/index.html#input-data-20240501/WW3_input_data_20240214/
* ``INPUTDATA_ROOT_BMIC``: https://noaa-ufs-regtests-pds.s3.amazonaws.com/index.html#BM_IC-20220207/
* ``INPUTDATA_LM4``: https://noaa-ufs-regtests-pds.s3.amazonaws.com/index.html#LM4_input_data

To download data, users must select the files they want from the bucket and download them either in their browser, via a ``wget`` command, or through the AWS CLI. 

Detailed information on input files can be found in :numref:`Chapter %s <InputsOutputs>`. 

.. _DownloadingWMCode:

==================================
Downloading the Weather Model Code
==================================

To clone the develop branch of the ``ufs-weather-model`` repository and update its submodules, execute the following commands:

.. code-block:: console

  git clone --recursive https://github.com/ufs-community/ufs-weather-model.git
  cd ufs-weather-model

Compiling the model will take place within the ``ufs-weather-model`` directory created by the clone command.

.. _build-wm:

==========================
Building the Weather Model
==========================

.. note:: 

   The most straightforward way to run the UFS WM is to use the regression testing (RT) framework. The RT framework will load modulefiles, build (compile) the desired WM configuration, and run the test(s). Users can create new tests or modify existing tests to correspond to the WM configuration(s) they wish to run. This section is provided for those who do not want to use the RT framework to run the WM. However, most users should skip to :numref:`Section %s <rt-config>` to learn more about RT configuration or :numref:`Section %s <run-wm>` to build/run the WM with the RT framework. 

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

.. _atmf:

**ATMF**

For the ``ufs-weather-model ATMF`` configuration (standalone ATM coupled to :term:`UFS Fire`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATMF -DCCPP_SUITES=FV3_HRRR -D32BIT=ON"

.. _atm_ds2s:

**ATM_DS2S**

For the ``ufs-weather-model ATM_DS2S`` configuration (:term:`ATM`/:term:`DOCN`/:term:`DICE`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM_DS2S  -DCCPP_SUITES=FV3_GFS_v17_coupled_p8_ugwpv1"


.. _atm_ds2s-pcice:

**ATM_DS2S-PCICE**

For the ``ufs-weather-model ATM_DS2S-PCICE`` configuration (:term:`ATM`/:term:`DOCN`/:term:`CICE6` [prescribed ice mode]):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM_DS2S-PCICE -DCCPP_SUITES=FV3_GFS_v17_coupled_p8"


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

.. _s2swal:

**S2SWAL**

For the ``ufs-weather-model S2SWAL`` configuration (atm/ice/ocean/wave/aerosols/land):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2SWAL -DCCPP_SUITES=FV3_GFS_v17_coupled_p8,FV3_GFS_v17_coupled_p8_ugwpv1"


.. _ng-godas:

NG-GODAS Configuration
------------------------

For the ``ufs-weather-model NG-GODAS`` configuration (atm/ocean/ice/data assimilation): 

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=NG-GODAS"

HAFS Configurations
----------------------

.. _hafs:

**HAFS**

For the ``ufs-weather-model HAFS`` configuration (atm/ocean) in 32 bit:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFS -D32BIT=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf_nonsst,FV3_HAFS_v0_gfdlmp_tedmf"

.. _hafsw:

**HAFSW**

For the ``ufs-weather-model HAFSW`` configuration (atm/:term:`HYCOM`/wave) in 32-bit with moving nest:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFSW -D32BIT=ON -DMOVING_NEST=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf,FV3_HAFS_v0_gfdlmp_tedmf_nonsst,FV3_HAFS_v0_thompson_tedmf_gfdlsf"

.. _hafs-mom6w:

**HAFS-MOM6W**

For the ``ufs-weather-model HAFS-MOM6`` configuration (atm/:term:`MOM6`/wave) in 32-bit with moving nest:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFS-MOM6W -DREGIONAL_MOM6=ON -DCDEPS_INLINE=ON -DMOVING_NEST=ON -DCCPP_SUITES=FV3_HAFS_v1_gfdlmp_tedmf,FV3_HAFS_v1_gfdlmp_tedmf_nonsst,FV3_HAFS_v1_thompson,FV3_HAFS_v1_thompson_nonsst -D32BIT=ON"

.. _hafs-all:

**HAFS-ALL**

For the ``ufs-weather-model HAFS-ALL`` configuration (data/atm/ocean/wave) in 32 bit:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFS-ALL -D32BIT=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf,FV3_HAFS_v0_gfdlmp_tedmf_nonsst"

Land Configurations
----------------------

.. _lnd:

**LND**

For the ``ufs-weather-model LND`` configuration (:term:`DATM`/land [:term:`NOAHMP`]):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=LND"

.. _lnd-lm4:

**LM4**

For the ``ufs-weather-model LND-LM4`` configuration (:term:`DATM`/land [:term:`LM4`]):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=LND-LM4"

------------------
Building the Model
------------------

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

.. _rt-config:

====================
Test Configuration
====================

.. note:: 
   
   This section explains how forecasts are configured using the regression test (RT) framework. For a full list of 
   supported RT configurations, view the `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`_ file
   or visit the :wm-repo:`tests/tests <tree/develop/tests/tests>` directory. 


The UFS Weather Model (WM) can be run in any of several configurations, from a single-component atmospheric 
model to a fully coupled model with multiple earth system components (e.g., atmosphere, ocean, sea-ice, land, and 
mediator). Each RT test configuration file (located in the :wm-repo:`tests/tests directory <tree/develop/tests/tests>`) 
sets default variables by calling functions from :wm-repo:`tests/default_vars.sh <blob/develop/tests/default_vars.sh>`. 
Then, the test file sets test-specific variables. These values will override 
the defaults. 

---------------------
``default_vars.sh`` 
---------------------

``default_vars.sh`` first sets a series of machine-specific variables. It also contains several functions that set defaults for different types of tests. :numref:`Table %s <def-funcs>` describes what each function does. 

.. _def-funcs:

.. list-table:: ``default_vars.sh`` functions
   :widths: 10 70
   :header-rows: 1
   
   * - Function Name
     - Description
   * - export_fv3_v16
     - Set variables to the FV3 default values for GFS v16 cases. This section will be removed once support for GFSv16 is officially depricated.
   * - export_fv3
     - Set variables to the FV3 default values.
   * - export_tiled
     - Set default values for tiled grid namelist.
   * - export_ugwpv1 
     - Set default values for the Unified Gravity Wave Drag Physics v1. 
   * - export_cice6
     - Set default values for the CICE6 model namelist and ``mx100``. 
   * - export_mom6 
     - Set default values for the MOM6 model namelist and ``mx100``. 
   * - export_ww3
     - Set default values for the WW3 global model. 
   * - export_fire_behavior
     - Set default values for the Fire Behavior model. 
   * - export_cmeps
     - Set default values for the coupled 5-component tests using CMEPS.
   * - export_cpl
     - Set default values for *coupled* / S2S configurations. 
   * - export_35d_run
     - Set default values for EMC's weekly coupled benchmark 35d tests (see `rt_35d.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt_35d.conf>`__). 
   * - export_datm_cdeps
     - Set default values for configurations that use the data atmosphere (:term:`DATM`) component. 
   * - export_hafs_datm_cdeps
     - Set default values for HAFS configurations that use the data atmosphere (DATM) component. 
   * - export_hafs_docn_cdeps
     - Set default values for HAFS configurations that use the data ocean (:term:`DOCN`) component. 
   * - export_hafs_regional
     - Set default values for regional HAFS configurations. 
   * - export_hafs
     - Set default values for HAFS configurations. 
   * - export_hrrr 
     - Set default values for HRRR test configurations. 
   * - export_hrrr_conus13km
     - Set default values for hrrr_conus13km test configurations. 
   * - export_rap_common
     - Set default values that are common to RAP and RRFS v1 test configurations. 
   * - export_rap
     - Set default values for RAP test configurations. 
   * - export_rrfs_v1
     - Set default values for RRFS v1 test configurations.
   
Multiple ``default_vars.sh`` functions may be called in a given test, usually starting with the most general function and ending with the most specific. Values set in one function will be overridden when the same values are set in a subsequent function. 

------------
Test Files
------------

Individual test files typically start with an ``export TEST_DESCR`` statement describing the test, followed by an ``export CNTL_DIR`` statement indicating the name of the directory that contains the baselines for the experiment. Next, an ``export LIST_FILES`` statement indicates which files the test expects to output from the model run. This list often includes RESTART files. After the LIST_FILES statement, the tests typically call functions from ``default_vars.sh`` to set default values. 

For example, the ``hafs_regional_atm_ocn_wav`` test file lists the files that it will output and then calls three ``export_*`` functions from ``default_vars.sh``, moving from the most general to the most specific:

.. code-block:: console

   export LIST_FILES="atmf006.nc \
                   sfcf006.nc \
                   archv.2019_241_06.a \
                   archs.2019_241_06.a \
                   20190829.060000.out_grd.ww3 \
                   20190829.060000.out_pnt.ww3 \
                   ufs.hafs.ww3.r.2019-08-29-21600.nc \
                   ufs.hafs.cpl.r.2019-08-29-21600.nc"

   export_fv3
   export_hafs
   export_hafs_regional

Lastly, the :wm-repo:`test configuration file <blob/develop/tests/tests/hafs_regional_atm_ocn_wav>` sets any test-specific variables for the experiment. These variables will override the default values from ``default_vars.sh``. In the excerpt below, ``...`` indicates omitted lines: 

.. code-block:: console

   export HAFS=true
   export FHMAX=6
   export RESTART_N=${FHMAX}
   export DT_ATMOS=180
   export IDEFLATE=1
   export OUTPUT_FH='3 -1'
   export OUTPUT_FILE="'netcdf' 'netcdf'"
   export SDAY=29
   export SHOUR=00
   export SMONTH=08
   export SYEAR=2019

   ...

   export CDEPS_DOCN=false
   export OCEAN_START_DTG=43340.00000

   export atm_model=fv3
   export ocn_model=hycom
   export wav_model=ww3
   OCN_tasks=60
   WAV_tasks=60
   export coupling_interval_sec=360
   export MESH_ATM=unset

   export FIELD_TABLE=field_table_hafs
   export DIAG_TABLE=diag_table_hafs_template
   export INPUT_NML=input_regional_hafs.nml.IN
   export MODEL_CONFIGURE=model_configure_hafs.IN
   export UFS_CONFIGURE=ufs.configure.hafs_atm_ocn_wav.IN
   export FV3_RUN="hafs_fv3_run.IN hycom_hat10_run.IN hafs_ww3_run.IN"

   if [[ $MACHINE_ID = orion ]]; then
   WLCLK=40
   fi
   ...

.. _new-test:

--------------------
Creating New Tests
--------------------

Users are welcome to modify current tests for their own use or create new tests to facilitate their own research. 
When creating a test, users will need to add a row for the test in ``rt.conf`` or in their own custom file. 
See :numref:`Section %s <rt-conf>` for more information. 

Typically, when a developer needs to create a new test for his/her implementation, the
first step would be to identify a test in the ``tests/tests`` directory that can
be used as a basis and to examine the variables defined in the test file. 
The names of appropriate template files for model configuration and initial conditions
can be identified via variables ``INPUT_NML``, ``UFS_CONFIGURE``, ``MODEL_CONFIGURE`` and ``FV3_RUN`` 
by running ``grep -n INPUT_NML *`` inside the ``tests`` and ``tests/tests`` directories.

.. _rt-conf:

-----------------------
The ``rt.conf`` File
-----------------------

The ``rt.conf`` file is a pipe-separated values (PSV) file grouped into sections of tests with a ``COMPILE`` line followed by several ``RUN`` lines. The ``COMPILE`` line contains information needed to compile the tests, while the ``RUN`` lines contain information on specific tests. 
``COMPILE`` lines have 6 columns:

   #. ``COMPILE`` indicator
   #. **Compile name** -- a category of test to compile
   #. **Compiler** to use in build (``intel`` or ``gnu``)
   #. **CMAKE Options** -- Provides all CMAKE options for the build. This typically includes the ``-DAPP`` and ``-DCCPP_SUITES`` flags; these flags set which components to build and which physics suites will be available at runtime. Additional options are documented in :numref:`Section %s <other-build-options>`, but users can examine the :wm-repo:`CMakeLists.txt <blob/develop/CMakeLists.txt>` file for the most up-to-date list of options. 
   #. **Machines** to run on (``-`` is used to ignore specified machines, ``+`` is used to run only on specified machines). For example: 
      
      * ``+ hera orion gaea``: Compile will only run on Hera, Orion, and Gaea machines
      * ``- wcoss2 acorn``: Compile will NOT be run on WCOSS2 or Acorn

   #. ``fv3``: Set as fv3. Previously, this was used to run a test without compiling code (e.g., if FV3 was already present). 

After each compile line is one or more ``RUN`` lines. ``RUN`` lines have five columns. The build resulting from the ``COMPILE`` line above the ``RUN`` line will be used to run the tests. 

   #. ``RUN`` indicator
   #. **Test name** -- indicates which test in the :wm-repo:`tests/tests <tree/develop/tests/tests>` directory should be sourced.
   #. **Machines** to run on (``+``) or ignore (``-``).
   #. **Baseline Creation** -- controls whether the run creates its own baseline or uses the baseline from a different (control) test (see information on ``-c`` option :ref:`below <cmd-line-opts>` for more).
   #. **Comparison Test** -- Test name to compare baselines with if not itself.

The order of lines in ``rt.conf`` matters since ``rt.sh`` processes them sequentially; a ``RUN`` line should be preceeded
by a ``COMPILE`` line that builds the model used in the test. The following
``rt.conf`` file excerpt builds the standalone ATM model with GFS_v16 physics 
in 32-bit mode and then runs the ``control`` test:

.. code-block:: console

   COMPILE | s2swa_32bit_pdlib  | intel | -DAPP=S2SWA -D32BIT=ON -DCCPP_SUITES=FV3_GFS_v17_coupled_p8_ugwpv1 -DPDLIB=ON | - noaacloud | fv3 |
   RUN | cpld_control_gfsv17                               | - noaacloud                          | baseline |
   RUN | cpld_control_gfsv17_iau                           | - noaacloud                          | baseline | cpld_control_gfsv17
   RUN | cpld_restart_gfsv17                               | - noaacloud                          |          | cpld_control_gfsv17
   RUN | cpld_mpi_gfsv17                                   | - noaacloud                          |          |

The ``rt.conf`` file includes a large number of tests. If the user wants to run
only specific tests, s/he can either (1) comment out the tests to be skipped (using the ``#`` prefix)
or (2) create a new file (e.g., ``my_rt.conf``), add the tests, and execute ``./rt.sh -l my_rt.conf``.


.. _run-wm:

=================
Running the Model
=================

.. attention::
   Although the following discussions are general, users may not be able to execute the script successfully "as is" unless they are on a 
   :wm-wiki:`Tier-1 platform <Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`.

.. _UsingRegressionTest:

--------------------------------
Using the Regression Test Script
--------------------------------

Users can run a number of preconfigured regression test cases from the ``rt.conf`` file 
using the regression test script ``rt.sh`` in the ``tests`` directory. 
``rt.sh`` is the top-level script that calls lower-level scripts to build specified 
WM configurations, set up environments, and run tests. 
Users should edit the ``rt.conf`` file to indicate which tests/configurations to run or create their own configuration file (e.g., ``my_tests.conf``) with the subset of tests they want to run. 

On NOAA RDHPCS
------------------

On :wm-wiki:`Tier-1 platforms <Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`, users can run 
regression tests by editing the ``rt.conf`` file and executing:

.. code-block:: console

    ./rt.sh -a <account> -l rt.conf

where ``<account>`` is to the account/project number where users submit their batch jobs. 
Users may need to add additional command line arguments or change information in the ``rt.sh`` file as well. 
This information is provided in :numref:`Section %s <rt.sh>` below. 

.. _other-systems:

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

   ./rt.sh -a <account> | -b <file> | -c | -d | -e | -h | -k | -l <file> | -m | -n <name> | -o | -r | -v | -w
      -a  <account> to use on for HPC queue
      -b  create new baselines only for tests listed in <file>
      -c  create new baseline results
      -d  delete run directories that are not used by other tests
      -e  use ecFlow workflow manager
      -h  display this help
      -k  keep run directory after rt.sh is completed
      -l  runs test specified in <file>
      -m  compare against new baseline results
      -n  run single test <name>
      -o  compile only, skip tests
      -r  use Rocoto workflow manager
      -v  verbose output
      -w  for weekly_test, skip comparing baseline results

When running a large number (10's or 100's) of tests, the ``-e`` or ``-r`` options can significantly
decrease testing time by using a workflow manager (ecFlow or Rocoto, respectively) to queue the jobs 
according to dependencies and run them concurrently. 
The ``-n`` option can be used to run a single test; for example, ``./rt.sh -a epic -n "control_c48 intel"`` 
will build the ATM model and run the ``control_c48`` test with an Intel compiler using the "epic" account 
(users should substitute an account where they can charge computational resources).
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

   ./rt.sh -r -k -n "control_p8 <compiler>"

where ``<compiler>`` is ``gnu`` or ``intel``. 

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
``model_configure``, ``ufs.configure``) and other application
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
tests, so it is recommended that the user first read :numref:`Section %s <UsingRegressionTest>`. 
All the files in the subdirectories shown in :numref:`Table %s <RTSubDirs>` are relevant to the
operational requirement test. The only difference is that the ``opnReqTest`` script replaces ``rt.sh``.
The ``tests/opnReqTests`` directory contains
opnReqTest-specific lower-level scripts used to set up run configurations.

On :wm-wiki:`Tier-1 platforms <Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`, tests can
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
 
   Usage: opnReqTest -n <test-name> -a <account> [ -c <test-case> ] [-b] [-d] [-e] [-k] [-h] [-x] [-z]
  
      -a  specify HPC <account> to use for batch job
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
workflow manager, and ``-k`` to keep the ``$RUNDIR``. The Rocoto workflow manager 
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
