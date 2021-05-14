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

To clone the ufs-weather-model repository for this v2.0.0 release, execute the following commands:

.. code-block:: console

  git clone https://github.com/ufs-community/ufs-weather-model.git ufs-weather-model
  cd ufs-weather-model
  git checkout ufs-v2.0.0
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
   |  nemsio          | export NEMSIO_INC=<path_to_nemsio_include_dir>                  |
   +------------------+-----------------------------------------------------------------+
   |                  | export NEMSIO_LIB=<path_to_nemsio_lib_dir>/libnemsio<version>.a |
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
modulefiles.  Modulefiles for all supported platforms are located in ``modulefiles/<platform>/fv3``. To
load the modules from the `ufs-weather-model` directory on hera:

.. code-block:: console

    cd modulefiles/hera.intel
    module use $(pwd)
    module load fv3
    cd ../..

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

--------------------------------------------
Setting the CCPP_SUITES environment variable
--------------------------------------------

In order to have one or more CCPP physics suites available at runtime, you need to select those suites at
build time by setting the ``CCPP_SUITES`` environment variable. Multiple suites can be set, as shown below
in an example for the bash shell:

.. code-block:: console

    export CCPP_SUITES="FV3_GFS_v15p2,FV3_GFS_v16beta"

If ``CCPP_SUITES`` is not set, the default is set to ``‘FV3_GFS_v15p2’`` in ``build.sh``.

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

--------------------------------
Using the regression test script
--------------------------------
The regression test script ``rt.sh`` in the tests/ directory can be
used to run a number of preconfigured test cases. It is the top-level script
that calls lower-level scripts to build, set up environments and run tests.
On `Tier-1 platforms <https://github.com/ufs-community/ ufs-weather-model/wiki
/Regression-Test-Policy-for-Weather-Model-Platforms-and-Compilers>`_, it can
be as simple as editing the ``rt.conf`` file and subsequently executing
``./rt.sh -l rt.conf``. Following discussions assume that the user has access
to at least one Tier-1 platform.

Each line in the PSV (Pipe-separated values) file ``rt.conf`` is used to either
build or test. The ``COMPILE`` line specifies the application to build (e.g.
``APP=S2S``), CCPP suite to use (e.g. ``SUITES=FV3_GFS_2017_coupled``), and
additional build options (e.g. ``DEBUG=Y``) as necessary. The ``RUN`` line
specifies the name of a test to run. The test name should match one of the test
files in the tests/tests/ directory or, if the user is adding a new test, the name
of the new test file. The order of lines in ``rt.conf`` matters since ``rt.sh``
processes them sequentially; a ``RUN`` line should be proceeded by a ``COMPILE``
line that builds the model used in the test. The following example ``rt.conf``
file builds the Subseasonal to Seasonal (S2S) model and then runs the
``cpld_control`` test.

.. code-block:: console

    COMPILE | APP=S2S SUITES=FV3_GFS_2017_coupled | | fv3
    RUN     | cpld_control                        | | fv3

Regression test generates a number of log files. The summary log file
``RegressionTests_<machine>.<compiler>.log`` in the tests/ directory compares
the results of the test against the baseline specific to a given platform and
reports the outcome (hence, the 'regression' test). More detailed log files are
found in tests/log_<machine> directory. Particularly, the user may find useful
the compile and run directory paths provided as the value of ``RUNDIR``
variable in run file. ``RUNDIR`` is a self-contained directory.


In the following, a more detailed description of ``rt.sh``, its lower-level
scripts, configuration templates for setting simulation parameters, and etc.
is provided. This may be useful when the user wants to create a new test case
to test his/her implementations.

Shown in :numref:`Table %s <LowLevelScripts>` are the lower-level scripts found
in tests/ directory.

.. _LowLevelScripts:

.. table:: *Lower-level scripts*

   +----------------------+-----------------------------------------------------------+
   | **File Name**        | **Description**                                           |
   +======================+===========================================================+
   | edit_inputs.sh       | Sets the default simulation parameters                    |
   +----------------------+-----------------------------------------------------------+
   | default_vars.sh      | Sets the default simulation parameters                    |
   +----------------------+-----------------------------------------------------------+
   | rt_utils.sh          | Sets the default simulation parameters                    |
   +----------------------+-----------------------------------------------------------+
   | detect_machine.sh    | Detects platform and sets relevant parameters             |
   +----------------------+-----------------------------------------------------------+
   | run_compile.sh       | Sets the default simulation parameters                    |
   +----------------------+-----------------------------------------------------------+
   | run_test.sh          | Sets the default simulation parameters                    |
   +----------------------+-----------------------------------------------------------+

Next,

.. _SubDirs:

.. table:: *Subdirectories*

   +---------------------------+-----------------------------------------------------------+
   | **Subdirectory Name**     | **Description**                                           |
   +===========================+===========================================================+
   | fv3_conf                  | Sets the default simulation parameters                    |
   +---------------------------+-----------------------------------------------------------+
   | parm                      | Sets the default simulation parameters                    |
   +---------------------------+-----------------------------------------------------------+
   | tests                     | Sets the default simulation parameters                    |
   +---------------------------+-----------------------------------------------------------+

--------------------------
Using the unit test script
--------------------------
The unit test script ``utest.sh`` in the tests/ directory can be used. Given the name of a test,
it carries out a suite of tests, with each test addressing a single aspect of the requirement
a new implementation should satisfy:

#. Thread (``THR``) reproducibility: varying the number of threads produces the same results
#. MPI process (``MPI``) reproducibility: varying the number of MPI tasks produces the same results
#. Domain decomposition (``DCP``) reproducibility: varying the tile layout of FV3 produces the same results
#. Restart (``RST``) reproducibility: restarting produces the same results
#. 64/32 (``BIT``) reproducibility: changing to double/single precision can be compiled and simulations can be run to completion
#. Debug (``DBG``) reproducibility: can be compiled and simulation can be run to completion in debug mode

The model can be run by executing ``./utest -n <test-name> -c <case-name>``, where ``<test-name>`` is
the name of a test file in tests/tests/ directory, and ``<case-name>`` is one or comma-separated combination
of ``THR``, ``MPI``, ``DCP``, ``RST``, ``BIT``, ``DBG``.
