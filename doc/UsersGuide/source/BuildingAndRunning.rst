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
   <https://github.com/NOAA-EMC/NCEPLIBS/tree/ufs-v1.1.0>`_ (**be sure to look at the tag in that repository that
   matches the tag on this WM release**).

#. Third-party libraries (NCEPLIBS-external). These are libraries that were developed external to
   the UFS Weather Model. They are general software packages that are also used by other models in
   the community. Building these is optional, since existing builds of these libraries can be pointed
   to instead. A list of the external libraries tested with this WM release is in the top-level ``README``
   of the `NCEPLIBS-external repository <https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v1.1.0>`_. Again, be
   sure to look at the tag in that repository that matches the tag on this WM release.

.. note::
   The libraries in NCEPLIBS-external must be built *before* the libraries in NCEPLIBS.

See this `wiki link <https://github.com/ufs-community/ufs/wiki/Supported-Platforms-and-Compilers>`_ for
an explanation of which platforms and compilers are supported. This will help to determine if you need
to build NCEPLIBS and NCEPLIBS-external or are working on a system that is already pre-configured. On
pre-configured platforms, the libraries are already available.

If you do have to build the libraries, it is a good idea to check the platform- and compiler-specific
``README`` files in the doc/ directory of the `NCEPLIBS-external repository <https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v 1.1.0>`_
as a first step, to see if your system or one similar to it is included. These files have detailed
instructions for building NCEPLIBS-external, NCEPLIBS, and the UFS Weather Model. They may be all the
documentation you need. Be sure to use the tag that corresponds to this version of the WM, and define a
WORK directory path before you get started.

If your platform is not included in these platform- and compiler-specific ``README`` files, there is a more
generic set of instructions in the ``README`` file at the top level of the `NCEPLIBS-external repository
<https://github.com/NOAA-EMC/NCEPLIBS-external/tree/ufs-v1.1.0>`_, and at the top level of the `NCEPLIBS repository
<https://github.com/NOAA-EMC/NCEPLIBS/tree/ufs-v1.1.0>`_. It may still be a good idea to look at some of the platform-
and compiler-specific ``README`` files as a guide. Again, be sure to use the tag that corresponds to this version of the WM.

The top-level ``README`` in the NCEPLIBS-external repository includes a troubleshooting section that may be helpful.

You can also get expert help through a `user support forum <https://forums.ufscommunity.org/forum/build-dependencies>`_
set up specifically for issues related to build dependencies.

.. _DownloadingWMCode:

==================================
Downloading the Weather Model Code
==================================

To clone the ufs-weather-model repository for this v1.1.0 release, execute the following commands:

.. code-block:: console

  git clone https://github.com/ufs-community/ufs-weather-model.git ufs-weather-model
  cd ufs-weather-model
  git checkout ufs-v1.1.0
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
to build the prerequisite libraries, there is a script in the ``NCEPLIBS-ufs-v1.1.0/bin`` directory called
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
The `UFS Weather Model wiki <https://github.com/ufs-community/ufs-weather-model/wiki>`_ includes a simple
test case that illustrates how the model can be run.
