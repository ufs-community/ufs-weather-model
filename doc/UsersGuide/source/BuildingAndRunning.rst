.. _BuildingAndRunning:

******************************************
Building and Running the UFS Weather Model
******************************************

======================
Prerequisite Libraries
======================

The UFS Weather Model requires a number of libraries for it to compile.
There are two categories of libraries that are needed:

#. Bundled libraries (NCEPLIBS). These are libraries developed for use with NOAA weather models.
   Most have an NCEPLIBS prefix in the repository, e.g. NCEPLIBS-bacio. Select tools from the UFS
   Utilities repository (UFS-UTILS) are also included in this category. A list of the bundled
   libraries tested with this WM release is in the top-level ``README`` of the `NCEPLIBS repository
   <https://github.com/NOAA-EMC/NCEPLIBS/tree/ufs-v2.0.0>`_ (**be sure to look at the tag in that repository that
   matches the tag on this WM release**).

#. Third-party libraries (NCEPLIBS-external). These are libraries that were developed external to
   the UFS Weather Model. They are general software packages that are also used by other models in
   the community. Building these is optional, since existing builds of these libraries can be pointed
   to instead. A list of the external libraries tested with this WM release is in the top-level ``README``
   of the `NCEPLIBS-external repository <https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v2.0.0>`_. Again, be
   sure to look at the tag in that repository that matches the tag on this WM release.

.. note::
   The libraries in NCEPLIBS-external must be built *before* the libraries in NCEPLIBS.

See this `wiki link <https://github.com/ufs-community/ufs/wiki/Supported-Platforms-and-Compilers>`_ for
an explanation of which platforms and compilers are supported. This will help to determine if you need
to build NCEPLIBS and NCEPLIBS-external or are working on a system that is already pre-configured. On
pre-configured platforms, the libraries are already available.

If you do have to build the libraries, it is a good idea to check the platform- and compiler-specific
``README`` files in the doc/ directory of the `NCEPLIBS-external repository <https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v 2.0.0>`_
as a first step, to see if your system or one similar to it is included. These files have detailed
instructions for building NCEPLIBS-external, NCEPLIBS, and the UFS Weather Model. They may be all the
documentation you need. Be sure to use the tag that corresponds to this version of the WM, and define a
WORK directory path before you get started.

If your platform is not included in these platform- and compiler-specific ``README`` files, there is a more
generic set of instructions in the ``README`` file at the top level of the `NCEPLIBS-external repository
<https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v2.0.0>`_, and at the top level of the `NCEPLIBS repository
<https://github.com/NOAA-EMC/NCEPLIBS/tree/ufs-v2.0.0>`_. It may still be a good idea to look at some of the platform-
and compiler-specific ``README`` files as a guide. Again, be sure to use the tag that corresponds to this version of the WM.

The top-level ``README`` in the NCEPLIBS-external repository includes a troubleshooting section that may be helpful.

You can also get expert help through a `user support forum <https://forums.ufscommunity.org/forum/build-dependencies>`_
set up specifically for issues related to build dependencies.

.. _DownloadingWMCode:

==================================
Downloading the Weather Model Code
==================================

To clone the ufs-weather-model repository, execute the following commands:

.. code-block:: console

  git clone https://github.com/ufs-community/ufs-weather-model.git ufs-weather-model
  cd ufs-weather-model
  git submodule update --init --recursive

Compiling the model will take place within the `ufs-weather-model` directory you just created.

==========================
Building the Weather Model
==========================

-------------------------------------------------------------------------
Setting environment variables for NCEPLIBS, NCEPLIBS-external and CMake
-------------------------------------------------------------------------
You will need to make sure that the WM has the paths to the libraries that it requires. In order to do
that, these environment variables need to be set, as shown in :numref:`Table %s <ReqLibEnvVar>` and
:numref:`Table %s <ReqLibEnvVar2>` for the bash shell.

.. _ReqLibEnvVar:

.. table:: *Bundled libraries (NCEPLIBS) required for the Weather Model*

   +------------------+-----------------------------------------------------------------+
   | **NCEP Library** | **Environment Variables**                                       |
   +==================+=================================================================+
   |  nemsio          || export NEMSIO_INC=<path_to_nemsio_include_dir>                 |
   |                  || export NEMSIO_LIB=<path_to_nemsio_lib_dir>/libnemsio<version>.a|
   +------------------+-----------------------------------------------------------------+
   |  bacio           | export BACIO_LIB4=<path_to_bacio_lib_dir>/libbacio<version>.a   |
   +------------------+-----------------------------------------------------------------+
   |  splib           | export SP_LIBd=<path_to_sp_lib_dir>/libsp<version>_d.a          |
   +------------------+-----------------------------------------------------------------+
   |  w3emc           | export W3EMC_LIBd=<path_to_w3emc_lib_dir>/libw3emc<version>_d.a |
   +------------------+-----------------------------------------------------------------+
   |  w3nco           | export W3NCO_LIBd=<path_to_w3nco_lib_dir>/libw3nco<version>_d.a |
   +------------------+-----------------------------------------------------------------+

|

.. _ReqLibEnvVar2:

.. table:: *Third-party libraries (NCEPLIBS-external) required for the Weather Model*

   +------------------+----------------------------------------------------+
   | **Library**      | **Environment Variables**                          |
   +==================+====================================================+
   |  NetCDF          | export NETCDF=<path_to_netcdf_install_dir>         |
   +------------------+----------------------------------------------------+
   |  ESMF            | export ESMFMKFILE=<path_to_esmfmk_file>/esmf.mk    |
   +------------------+----------------------------------------------------+

The following are a few different ways to set the required environment variables to the correct values.
If you are running on one of the `pre-configured platforms
<https://github.com/ufs-community/ufs/wiki/Supported-Platforms-and-Compilers>`_, you can set them using
modulefiles.  Modulefiles for all supported platforms are located in ``modulefiles/ufs_<platform>.<compiler>``. To
load the modules from the `ufs-weather-model` directory on hera:

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
in a different way.

If you used one of the platform- and compiler-specific ``README`` files in the ``doc/`` directory of NCEPLIBS-external
to build the prerequisite libraries, there is a script in the ``NCEPLIBS-ufs-v2.0.0/bin`` directory called
``setenv_nceplibs.sh`` that will set the NCEPLIBS-external variables for you.

Of course, you can also set the values of these variables yourself if you know where the paths are on your system.

-------------------------------------------------------------
Setting the CMAKE_FLAGS and CCPP_SUITES environment variables
-------------------------------------------------------------

You need to use the ``CMAKE_FLAGS`` environment variable to specify which application to build.
In order to have one or more CCPP physics suites available at runtime, you also need to select those suites at
build time by setting the ``CCPP_SUITES`` environment variable. Multiple suites can be set. Following
examples are for the bash shell.

For the ufs-weather-model ATM app (standalone ATM):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM"
    export CCPP_SUITES="FV3_GFS_v16"

For the ufs-weather-model ATM app (standalone ATM) in 32 bit:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATM -D32BIT=ON"
    export CCPP_SUITES="FV3_GFS_v16"

For the ufs-weather-model ATMW app (standalone ATM with wave):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=ATMW"
    export CCPP_SUITES="FV3_GFS_v16"

For the ufs-weather-model S2S app (atm/ice/ocean):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2S"
    export CCPP_SUITES="FV3_GFS_2017_coupled,FV3_GFS_2017_satmedmf_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled,FV3_GFS_v16_couplednsst"

For the ufs-weather-model S2S app (atm/ice/ocean) with debugging flags turned on, with verbose build messages:

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2S -DDEBUG=ON"
    export CCPP_SUITES="FV3_GFS_2017_coupled,FV3_GFS_2017_satmedmf_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled,FV3_GFS_v16_couplednsst"

For the ufs-weather-model S2SW app (atm/ice/ocean/wave):

.. code-block:: console

    export CMAKE_FLAGS="-DAPP=S2SW"
    export CCPP_SUITES="FV3_GFS_2017_coupled,FV3_GFS_v15p2_coupled,FV3_GFS_v16_coupled,FV3_GFS_v16_coupled_noahmp"

------------------
Building the model
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

Once ``build.sh`` is finished, you should see the executable, named ``ufs_weather_model``, in the top-level directory.

Expert help is available through a `user support forum <https://forums.ufscommunity.org/forum/ufs-weather-model>`_
set up specifically for issues related to the Weather Model.

=================
Running the model
=================

.. _UsingRegressionTest:

--------------------------------
Using the regression test script
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
``rt.conf`` file builds the Subseasonal to Seasonal (S2S) model and then runs the
``cpld_control`` test:

.. code-block:: console

    COMPILE | APP=S2S SUITES=FV3_GFS_2017_coupled | | fv3
    RUN     | cpld_control                        | | fv3

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
invoking the command ``sbatch job_card``. Note that ``$RUNDIR`` is automatically
deleted at the end of a successful regression test; specifying the ``-k`` option
retains the ``$RUNDIR``, e.g. ``./rt.sh -l rt.conf -k``.

Found inside the ``$RUNDIR`` directory are the model configuration files
``data_table``, ``diag_table``, ``ice_in``, ``input.nml``, ``model_configure``
and ``nems.configure``. They are generated by ``rt.sh`` from the template files
in the tests/parm/ directory. Specific values used to fill in the template files
depend on the test being run, and are set in two stages: default values are
specified in ``tests/default_vars.sh`` and the default values are overriden if
necessary by those specified in a test file ``tests/tests/<test-name>``. For
example, the variable ``DT_ATMOS``, which is substituted into the template file
``model_configure.IN`` to generate ``model_configure``, is initially assigned
1800 in the function ``export_fv3`` of the script ``default_vars.sh``, but the
test file ``tests/tests/control`` overrides by reassigning 720 to the variable.

Also found inside the ``$RUNDIR`` directory are the files ``fv3_run`` and
``job_card``, which are generated from the template files in the tests/fv3_conf/
directory. The latter is a platform-specific batch job submission script, while
the former prepares the initial conditions by copying relevant data from the
input data directory of a given platform to the ``$RUNDIR`` directory.
:numref:`Table %s <RTSubDirs>` summarizes the subdirectories discussed above.

.. _RTSubDirs:

.. table:: *Regression test subdirectories*

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
for example, ``./rt.sh -n cpld_control`` will build the S2S model and run the
``cpld_control`` test. The ``-c`` option is used to create baseline. New
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

.. _UsingUnitTest:

--------------------------
Using the unit test script
--------------------------
The unit test script ``utest`` in the tests/ directory can also be used to run
tests. Given the name of a test, ``utest`` carries out a suite of test cases.
Each test case addresses an aspect of the requirements new implementations
should satisfy, which are shown in :numref:`Table %s <ImplementationRequirement>`.
For the following discussions on utest, the user should note the distinction between
'test name' and 'test case': examples of test name are ``control``, ``cpld_control``
and ``regional_control`` which are all found in the /tests/tests/ directory, whereas
test case refers to any one of ``thr``, ``mpi``, ``dcp``, ``rst``, ``bit`` and ``dbg``.

.. _ImplementationRequirement:

.. table:: *Implementation requirements*

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

The unit test uses the same testing framework used by the regression
test, and therefore it is recommened that the user first read the
:numref:`Section on regression test %s <UsingRegressionTest>`. All the files in
the subdirectories shown in :numref:`Table %s <RTSubDirs>` are relavant to the
unit test except that the ``utest`` script replaces ``rt.sh`` and the
``utest.bld`` file replaces ``rt.conf``. The /tests/utests/ directory contains
utest-specific lower-level scripts used to set up run configurations.

On `Tier-1 platforms <https://github.com/ufs-community/ ufs-weather-model/wiki
/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`_, tests can
be run by first modifying the PSV file ``utest.bld`` to specify the build options
and then invoking

.. code-block:: console

    ./utest -n <test-name>

For example, including in the ``utest.bld`` file the following line

.. code-block:: console

    cpld_control | APP=S2S SUITES=FV3_GFS_2017_coupled

and then executing ``./utest -n cpld_control`` performs all six test cases
listed in :numref:`Table %s <ImplementationRequirement>` for ``cpld_control``
test. At the end of the run, a log file ``UnitTests_<machine>.<compiler>.log``
is generated in tests/ directory, which informs the user whether each test case
passed or failed. The user can choose to run a specific test case by invoking

.. code-block:: console

    ./utest -n <test-name> -c <test-case>

where ``<test-case>`` is one or
more comma-separated values selected from ``thr``, ``mpi``, ``dcp``, ``rst``,
``bit``, ``dbg``. For example, ``./utest -n cpld_control -c thr,rst`` runs the
``cpld_control`` test and checks the reproducibility of threading and restart.
The user can see different command line options available to ``utest`` by
executing ``./utest -h``; frequently used options are ``-e`` to use the ecFlow
workflow manager, and ``-k`` to keep the ``$RUNDIR``. In the following,
comparisons are made between the regression and unit tests on how they handle
different reproducibility tests.

As discussed in :numref:`Section %s <UsingRegressionTest>`, the variables and
values used to configure model parameters and to set up initial conditions in the
``$RUNDIR`` directory are set up by the regression test script ``rt.sh`` in two
stages: first, ``tests/default_vars.sh`` define default values; then a specific
test file in the tests/tests/ subdirectory either overrides the default values or
creates new variables if required by the test. The regression test treats
the different test cases shown in :numref:`Table %s <ImplementationRequirement>`
as different tests. Therefore, each test case requires a test file in the tests/tests/
subdirectory; examples are ``control_2threads``, ``control_decomp``,
``control_restart`` and ``control_debug``, which are just variations of ``control``
test to check various reproducibilities. There are two potential issues with this
approach. First, if several different variations of a given test were to be created
and included in the ``rt.conf`` file, there are too many tests to run. Second, if
a new test is added by the user, s/he will also have to create these variations.
The idea behind the unit test is to automatically configure and run these variations,
or test cases, given a test file. For example, ``./utest -n control`` will run all
six test cases in :numref:`Table %s <ImplementationRequirement>` based on a single
``control`` test file. Similarly, if the user adds a new test ``new_test``, then
``./utest -n new_test`` will run all test cases. This is done by the unit test script
``utest`` by adding a third stage of variable overrides, and the related scripts can
be found in the tests/utests/ directory.
