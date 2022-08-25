.. _BuildingAndRunning:

******************************************
Building and Running the UFS Weather Model
******************************************

===================================
Supported Platforms & Compilers
===================================
Before running the Weather Model (WM), users should determine which of the 
`four levels of support <https://github.com/ufs-community/ufs/wiki/Supported-Platforms-and-Compilers>`__ 
is applicable to their system. Generally, Level 1 & 2 systems are restricted to those with access 
through NOAA and its affiliates. These systems are named (e.g., Hera, Orion, Cheyenne). 
Level 3 & 4 systems include certain personal computers or non-NOAA-affiliated HPC systems. 
The prerequisite software libraries for building the WM already exist on Level 1/preconfigured 
systems, so users may skip directly :ref:`downloading the code <DownloadingWMCode>`. 
On other systems, users will need to build the prerequisite libraries using 
:term:`HPC-Stack` or :term:`spack-stack`. 

..
   COMMENT: Update link w/supported platforms and compilers!

======================
Prerequisite Libraries
======================

The UFS Weather Model (WM) requires a number of libraries for it to compile.
The WM uses two categories of libraries, which are available as a bundle via 
:term:`HPC-Stack` or :term:`spack-stack`:

   #. :term:`NCEP` libraries (NCEPLIBS): These are libraries developed for use with NOAA weather models.
      Most have an NCEPLIBS prefix in the repository (e.g., NCEPLIBS-bacio). Select tools from the UFS
      Utilities repository (:term:`UFS_UTILS`) are also included in this category. 
      A list of the bundled libraries tested with this WM release is available in the top-level ``README`` of the 
      `NCEPLIBS repository <https://github.com/NOAA-EMC/NCEPLIBS/tree/ufs-v2.0.0>`__ (**be sure to look at 
      the tag in that repository that matches the tag on the most recent WM release**).

   #. Third-party libraries (NCEPLIBS-external): These are libraries that were developed external to
      the UFS Weather Model. They are general software packages that are also used by other community models. 
      Building these libraries is optional if users can point to existing builds of these libraries on their system
      instead. A list of the external libraries tested with this WM release is in the top-level ``README``
      of the `NCEPLIBS-external <https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v2.0.0>`__ repository. Again, be
      sure to look at the tag in that repository that matches the tag on this WM release.

.. note::
   Documentation is available for installing `HPC-Stack <https://hpc-stack.readthedocs.io/en/latest/>`__ 
   and `spack-stack <https://spack-stack.readthedocs.io/en/latest/>`__, respectively. 
   One of these software stacks (or the libraries they contain) must be installed before running the Weather Model. 

For users who *do* need to build the prerequisite libraries, it is a good idea to check the platform- and compiler-specific
``README`` files in the ``doc`` directory of the `NCEPLIBS-external repository <https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v2.0.0>`_
first, to see if their system or one similar to it is included. These files have detailed
instructions for building NCEPLIBS-external, NCEPLIBS, and the UFS Weather Model. They may be all the
documentation you need. Be sure to use the tag that corresponds to this version of the WM, and define a
WORK directory path before you get started.

..
   COMMENT: What is meant by a WORK directory path?

If your platform is not included in these platform- and compiler-specific ``README`` files, there is a more
generic set of instructions in the ``README`` file at the top level of the `NCEPLIBS-external repository
<https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v2.0.0>`__ and at the top level of the `NCEPLIBS repository
<https://github.com/NOAA-EMC/NCEPLIBS/tree/ufs-v2.0.0>`__. It may still be a good idea to look at some of the platform-
and compiler-specific ``README`` files as a guide. Again, be sure to use the tag that corresponds to this version of the WM.

The top-level ``README`` in the NCEPLIBS-external repository includes a troubleshooting section that may be helpful.

You can also get expert help through a `user support forum <https://forums.ufscommunity.org/forum/build-dependencies>`__
set up specifically for issues related to build dependencies.

..
   COMMENT: Will this be deprecated? When?

.. _DownloadingWMCode:

==================================
Downloading the Weather Model Code
==================================

To clone the develop branch of the ``ufs-weather-model`` repository and update its submodules, execute the following commands:

.. code-block:: console

  git clone https://github.com/ufs-community/ufs-weather-model.git ufs-weather-model
  cd ufs-weather-model
  git submodule update --init --recursive

Compiling the model will take place within the ``ufs-weather-model`` directory you just created.

==========================
Building the Weather Model
==========================

----------------------------
Loading the Required Modules
----------------------------

Modulefiles for `pre-configured platforms <https://github.com/ufs-community/ufs/wiki/Supported-Platforms-and-Compilers>`_
are located in ``modulefiles/ufs_<platform>.<compiler>``. For example, to load the modules from the ``ufs-weather-model``
directory on Hera:

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

If you are not running on one of the pre-configured platforms, you will need to set the environment variables
manually. For example, in a bash shell, a command in the following form will set the C compiler environment variable:

.. code-block:: console

   export CMAKE_C_COMPILER=</path/to/C/compiler>


------------------------------------------------------------------------
Setting the ``CMAKE_FLAGS`` and ``CCPP_SUITES`` Environment Variables
------------------------------------------------------------------------

The UFS Weather Model can be built in one of twelve configurations (cf. :numref:`Table %s <UFS-configurations>`). 
The ``CMAKE_FLAGS`` environment variable specifies which configuration to build.
Additionally, users must select the :term:`CCPP` suite(s) by setting the ``CCPP_SUITES`` environment variable at
build time in order to have one or more CCPP physics suites available at runtime. Multiple suites can be set. 
Additional environment variables, such as ``-D32BIT=ON``, can be set if the user chooses. These options are documented 
in :numref:`Section %s <other-build-options>`. 
The following examples assume a bash shell.

ATM Configurations
---------------------

For the ``ufs-weather-model ATM`` configuration (standalone :term:`ATM`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16"

For the ``ufs-weather-model ATMW`` configuration (standalone ATM coupled to :term:`WW3`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATMW -DCCPP_SUITES=FV3_GFS_2017_coupled,FV3_GFS_v16_coupled"

.. CHECK above!!

For the ``ufs-weather-model ATMAERO`` configuration (standalone ATM coupled to :term:`GOCART`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATMAERO -DCCPP_SUITES=FV3_GFS_v16"

For the ``ufs-weather-model ATMAQ`` configuration (standalone ATM coupled to :term:`CMAQ`):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATMAQ -DCCPP_SUITES=FV3_GFS_v15p2"


S2S Configurations 
----------------------

For the ``ufs-weather-model S2S`` configuration (coupled atm/ice/ocean):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2S -DCCPP_SUITES=FV3_GFS_2017_coupled,FV3_GFS_2017_satmedmf_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled,FV3_GFS_v16_couplednsst"

.. Which ocean model is it coupled to? All? 

To turn on debugging flags, add ``-DDEBUG=ON`` flag after ``-DAPP=S2S``. Users can allow verbose build messages by running: 

.. code-block:: console

    export BUILD_VERBOSE=1

For the ``ufs-weather-model S2S`` configuration (atm/ice/ocean) with activating CCPP host model under CMEPS and receiving atmosphere-ocean fluxes from mediator:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2S -DCCPP_SUITES=FV3_GFS_v17_coupled_p8_sfcocn -DCMEPS_AOFLUX=ON"

..
   COMMENT: Need some clarification on what the above code does with CCPP/CMEPS... not clear from description. 

For the ``ufs-weather-model S2SW`` configuration (atm/ice/ocean/wave):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2SW -DCCPP_SUITES=FV3_GFS_2017_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled,FV3_GFS_v16_coupled_noahmp"

For the ``ufs-weather-model S2SA`` configuration (atm/ice/ocean/aerosols):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2SA -DCCPP_SUITES=FV3_GFS_2017_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled,FV3_GFS_v16_coupled_noahmp"

..
   CHECK: DAPP flag and physics suites

For the ``ufs-weather-model S2SWA`` configuration (atm/ice/ocean/wave/aerosols):


.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2SWA -DCCPP_SUITES=FV3_GFS_2017_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled,FV3_GFS_v16_coupled_noahmp"

..
   CHECK: physics suites

NG-GODAS Configuration
------------------------

For the ``ufs-weather-model NG-GODAS`` configuration (atm/ocean/ice/data assimilation): 


.. code-block:: console

    export CMAKE_FLAGS="-DAPP=NG-GODAS -DCCPP_SUITES=FV3_GFS_2017_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled"
..
   COMMENT: NG-GODAS --> Coupled CDEPS-DATM-MOM6-CICE6-CMEPS
   What is the DAPP argument? And the physics suites?

HAFS Configurations
----------------------

For the ``ufs-weather-model HAFS`` configuration (atm/ocean) in 32 bit:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFS -D32BIT=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf_nonsst,FV3_HAFS_v0_gfdlmp_tedmf,FV3_HAFS_v0_hwrf_thompson,FV3_HAFS_v0_hwrf"

For the ``ufs-weather-model HAFSW`` configuration (atm/ocean/wave) in 32 bit:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFSW -D32BIT=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf_nonsst,FV3_HAFS_v0_gfdlmp_tedmf,FV3_HAFS_v0_hwrf_thompson,FV3_HAFS_v0_hwrf"

For the ``ufs-weather-model HAFS-ALL`` configuration (data/atm/ocean/wave) in 32 bit:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=HAFS-ALL -D32BIT=ON -DCCPP_SUITES=FV3_HAFS_v0_gfdlmp_tedmf_nonsst,FV3_HAFS_v0_gfdlmp_tedmf,FV3_HAFS_v0_hwrf_thompson,FV3_HAFS_v0_hwrf"

------------------
Building the Model
------------------
The UFS Weather Model uses the CMake build system.  There is a build script called ``build.sh`` in the
top-level directory of the WM repository that configures the build environment and runs the ``make``
command.  This script also checks that all necessary environment variables have been set.

If any of the environment variables have not been set, the ``build.sh`` script will exit with a message similar to:

.. code-block:: console

   ./build.sh: line 11: CMAKE_Platform: Please set the CMAKE_Platform environment variable, e.g. [macosx.gnu|linux.gnu|linux.intel|hera.intel|...]

The WM can be built by running the following command from the `ufs-weather-model` directory:

.. code-block:: console

   ./build.sh

Once ``build.sh`` is finished, you should see the executable, named ``ufs_model``, in the `ufs-weather-model/build/` directory.
If it is desired to build in a different directory, specify the ``BUILD_DIR`` environment variable: e.g. ``export BUILD_DIR=test_cpld``
will build in the `ufs-weather-model/test_cpld` directory instead.

Expert help is available through a `user support forum <https://forums.ufscommunity.org/forum/ufs-weather-model>`__
set up specifically for issues related to the Weather Model.

=================
Running the Model
=================

.. _UsingRegressionTest:

--------------------------------
Using the Regression Test Script
--------------------------------
The regression test script ``rt.sh`` in the tests/ directory can be
used to run a number of preconfigured test cases. It is the top-level script
that calls lower-level scripts to build, set up environments and run tests.
On `Tier-1 platforms <https://github.com/ufs-community/ ufs-weather-model/wiki
/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`_, it can
be as simple as editing the ``rt.conf`` file and subsequently executing

.. code-block:: console

    ./rt.sh -l rt.conf

Following discussions are general, but the user may not be able to successfully
execute the script as is unless s/he is on one of the Tier-1 platforms.

Each line in the PSV (Pipe-separated values) file ``rt.conf`` is used to either
build or run. The ``COMPILE`` line specifies the application to build (e.g.
``APP=S2S``), CCPP suite to use (e.g. ``SUITES=FV3_GFS_2017_coupled``), and
additional build options (e.g. ``DEBUG=Y``) as necessary. The ``RUN`` line
specifies the name of a test to run. The test name should match the name of one
of the test files in the tests/tests/ directory or, if the user is adding a new
test, the name of the new test file. The order of lines in ``rt.conf`` matters
since ``rt.sh`` processes them sequentially; a ``RUN`` line should be proceeded
by a ``COMPILE`` line that builds the model used in the test. The following example
``rt.conf`` file builds the standalone ATM model in 32 bit and then runs the
``control`` test:

.. code-block:: console

    COMPILE | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16 -D32BIT=ON | | fv3
    RUN     | control                                        | | fv3

The third column of ``rt.conf`` relates to the platform; if left blank, the test
runs on all Tier-1 platforms. The fourth column deals with baseline creation (more
on this later) and ``fv3`` means the test will be included during baseline creation.
The ``rt.conf`` file includes a large number of tests. If the user wants to run
only a specific test, s/he can either comment out (using the ``#`` prefix) the
tests to be skipped, or create a new file, e.g. ``my_rt.conf``, then execute
``./rt.sh -l my_rt.conf``.

The regression test generates a number of log files. The summary log file
``RegressionTests_<machine>.<compiler>.log`` in the tests/ directory compares
the results of the test against the baseline specific to a given platform and
reports the outcome (hence, the 'regression' test): 'Missing file' results when
the expected files from the simulation are not found, and typically occurs
when the simulation did not run to completion; 'OK' means that the simulation
results are bit-for-bit identical to those of the baseline; 'NOT OK' when
the results are not bit-for-bit identical; and 'Missing baseline' when there
is no baseline data to compare against.

More detailed log files are found in the tests/log_<machine>.<compiler>/ directory.
In particular, the user may find useful the run directory path provided as the
value of ``RUNDIR`` variable in the ``run_<test-name>`` file. ``$RUNDIR`` is a
self-contained (i.e. sandboxed) directory with the executable file, initial
conditions, model configuration files, environment setup scripts and a batch job
submission script. The user can run the test by cd'ing into ``$RUNDIR`` and
invoking the command

.. code-block:: console

    sbatch job_card

This can be particularly useful for debugging and testing code changes. Note that
``$RUNDIR`` is automatically deleted at the end of a successful regression test;
specifying the ``-k`` option retains the ``$RUNDIR``, e.g. ``./rt.sh -l rt.conf -k``.

Found inside the ``$RUNDIR`` directory are a number of model configuration files:
``input.nml``, ``model_configure``, ``nems.configure``, and other application
dependent files, e.g. ``ice_in`` for Subseasonal-to-Seasonal application.
These model configuration files are
generated by ``rt.sh`` from the template files in the tests/parm/ directory.
Specific values used to fill in the template files depend on the test being run, and
are set in two stages: default values are specified in ``tests/default_vars.sh`` and
the default values are overriden if necessary by those specified in a test file
``tests/tests/<test-name>``. For example, the variable ``DT_ATMOS``, which is
substituted into the template file ``model_configure.IN`` to generate
``model_configure``, is initially assigned 1800 in the function ``export_fv3`` of the
script ``default_vars.sh``, but the test file ``tests/tests/control`` overrides by
reassigning 720 to the variable.

Also found inside the ``$RUNDIR`` directory are the files ``fv3_run`` and
``job_card``, which are generated from the template files in the tests/fv3_conf/
directory. The latter is a platform-specific batch job submission script, while
the former prepares the initial conditions by copying relevant data from the
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

There are a number of command line options available to the ``rt.sh`` script.
The user can execute ``./rt.sh`` to see information on these options. A couple
of them are discussed here. When running a large number (10's or 100's) of
tests, the ``-e`` option to use the ecFlow workflow manager can significantly
decrease the testing time by queuing the jobs according to dependencies and
running them concurrently. The ``-n`` option can be used to run a single test;
for example, ``./rt.sh -n control`` will build the ATM model and run the
``control`` test. The ``-c`` option is used to create baseline. New
baslines are needed when code changes lead to result changes, and therefore
deviate from existing baselines on a bit-for-bit basis.

When a developer needs to create a new test for his/her implementation, the
first step would be to identify a test in the tests/tests/ directory that can
be used as a basis and to examine the variables defined in the test file. As
mentioned above, some of the variables may be overrides for those defined in
``default_vars.sh``; others may be new variables that are needed specifically
for the test. Default variables and their values are defined in the ``export_fv3``
function of the ``default_vars.sh`` script for ATM application, ``export_cpl``
function for S2S application and ``export_datm`` function for GODAS application.
Also, the names of template files for model configuration and initial conditions
can be identified via variables ``INPUT_NML``, ``NEMS_CONFIGURE`` and ``FV3_RUN``;
for example, by trying ``grep -n INPUT_NML *`` inside the tests/ and tests/tests/
directories.

.. _UsingOpnReqTest:

---------------------------------------------
Using the operational requirement test script
---------------------------------------------
The operational requirement test script ``opnReqTest`` in the tests/ directory can also be used to run
tests. Given the name of a test, ``opnReqTest`` carries out a suite of test cases.
Each test case addresses an aspect of the requirements new implementations
should satisfy, which are shown in :numref:`Table %s <OperationalRequirement>`.
For the following discussions on opnReqTest, the user should note the distinction between
'test name' and 'test case': examples of test name are ``control``, ``cpld_control``
and ``regional_control`` which are all found in the /tests/tests/ directory, whereas
test case refers to any one of ``thr``, ``mpi``, ``dcp``, ``rst``, ``bit`` and ``dbg``.

.. _OperationalRequirement:

.. table:: *Operational Requirements*

  +----------+------------------------------------------------------------------------+
  | **Case** | **Description**                                                        |
  +==========+========================================================================+
  | thr      | Varying the number of threads produces the same results                |
  +----------+------------------------------------------------------------------------+
  | mpi      | Varying the number of MPI tasks reproduces                             |
  +----------+------------------------------------------------------------------------+
  | dcp      | Varying the decomposition (i.e. tile layout of FV3) reproduces         |
  +----------+------------------------------------------------------------------------+
  | rst      | Restarting reproduces                                                  |
  +----------+------------------------------------------------------------------------+
  | bit      | Model can be compiled in double/single precision and run to completion |
  +----------+------------------------------------------------------------------------+
  | dbg      | Model can be compiled and run to completion in debug mode              |
  +----------+------------------------------------------------------------------------+

The operational requirement test uses the same testing framework used by the regression
test, and therefore it is recommened that the user first read
:numref:`Section %s <UsingRegressionTest>`. All the files in
the subdirectories shown in :numref:`Table %s <RTSubDirs>` are relavant to the
operational requirement test except that the ``opnReqTest`` script replaces ``rt.sh``.
The /tests/opnReqTests/ directory contains
opnReqTest-specific lower-level scripts used to set up run configurations.

On `Tier-1 platforms <https://github.com/ufs-community/ ufs-weather-model/wiki
/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`_, tests can
be run by invoking

.. code-block:: console

    ./opnReqTest -n <test-name>

For example, ``./opnReqTest -n control`` performs all six test cases
listed in :numref:`Table %s <OperationalRequirement>` for ``control``
test. At the end of the run, a log file ``OpnReqTests_<machine>.<compiler>.log``
is generated in tests/ directory, which informs the user whether each test case
passed or failed. The user can choose to run a specific test case by invoking

.. code-block:: console

    ./opnReqTest -n <test-name> -c <test-case>

where ``<test-case>`` is one or
more comma-separated values selected from ``thr``, ``mpi``, ``dcp``, ``rst``,
``bit``, ``dbg``. For example, ``./opnReqTest -n control -c thr,rst`` runs the
``control`` test and checks the reproducibility of threading and restart.
The user can see different command line options available to ``opnReqTest`` by
executing ``./opnReqTest -h``; frequently used options are ``-e`` to use the ecFlow
workflow manager, and ``-k`` to keep the ``$RUNDIR``. In the following,
comparisons are made between the regression and operational requirement tests on how they handle
different reproducibility tests.

As discussed in :numref:`Section %s <UsingRegressionTest>`, the variables and
values used to configure model parameters and to set up initial conditions in the
``$RUNDIR`` directory are set up in two stages: first, ``tests/default_vars.sh``
define default values; then a specific test file in the tests/tests/ subdirectory
either overrides the default values or creates new variables if required by the test.
The regression test treats the different test cases shown in
:numref:`Table %s <OperationalRequirement>` as different tests. Therefore, each
test case requires a test file in the tests/tests/ subdirectory; examples are
``control_2threads``, ``control_decomp``, ``control_restart`` and ``control_debug``,
which are just variations of ``control`` test to check various reproducibilities.
There are two potential issues with this approach. First, if several different
variations of a given test were to be created and included in the ``rt.conf`` file,
there are too many tests to run. Second, if a new test is added by the user, s/he
will also have to create these variations. The idea behind the operational requirement test is to
automatically configure and run these variations, or test cases, given a test file.
For example, ``./opnReqTest -n control`` will run all six test cases in
:numref:`Table %s <OperationalRequirement>` based on a single ``control`` test file.
Similarly, if the user adds a new test ``new_test``, then ``./opnReqTest -n new_test`` will
run all test cases. This is done by the operational requirement test script ``opnReqTest`` by adding a third
stage of variable overrides, and the related scripts can be found in the tests/opnReqTests/
directory.
