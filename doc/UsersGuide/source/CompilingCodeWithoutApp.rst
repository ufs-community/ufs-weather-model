.. _CompilingCodeWithoutApp:
  
*****************************************
Compiling the Code without an Application
*****************************************

UFS weather model uses the cmake build system.  Since UFS weather model requires several external libraries to be available and appropriate environment variables pointing to those libraries to be set, a build script named ``build.sh`` is provided In the top-level directory of the repository to ensure all necessary variables are actually set.

The required libraries and environment variables can be set as shown in :numref:`Table %s <ReqLibEnvVar>` and :numref:`Table %s <ReqLibEnvVar2>` for the bash shell.

1. NCEP libraries:

.. _ReqLibEnvVar: 

.. list-table:: *List of NCEP libraries that comprise the ufs-weather-model*
   :widths: 20 80
   :header-rows: 1

   * - NCEP Library
     - Environment Variables
   * - nemsio
     - export NEMSIO_INC=<path_to_nemsio_include_dir> 
   * -
     - export NEMSIO_LIB=<path_to_nemsio_lib_dir>/libnemsio.a
   * - bacio
     - export BACIO_LIB4=<path_to_bacio_lib_dir>/libbacio.a
   * - splib
     - export SP_LIBd=<path_to_sp_lib_dir>/libsp_d.a
   * - w3emc
     - export W3EMC_LIBd=<path_to_w3emc_lib_dir>/libw3emc_d.a
   * - w3nco
     - export W3NCO_LIBd=<path_to_w3nco_lib_dir>/libw3nco_d.a

2. Third party libraries:

.. _ReqLibEnvVar2: 

.. list-table:: *List of External libraries that comprise the ufs-weather-model*
   :widths: 20 80
   :header-rows: 1

   * - Third Party Library
     - Environment Variables
   * - NetCDF
     - export NETCDF=<path_to_netcdf_install_dir>
   * - ESMF
     - export ESMFMKFILE=<path_to_esmfmk_file>/esmf.mk

You can either manually set those environment variables in your shell, or if you are on one of the supported platforms (Tier 1) the easiest way of setting those variables is by using environment module and load all required modules. Modulefiles for all supported platforms are located in ``modulefiles/<platform>/fv3``. To load those modules, for example on hera, run:

..  code-block:: console

    $ cd modulefiles/hera.intel
    $ module use $(pwd)
    $ module load fv3
    $ cd ../..

Currently the build system supports two compiler families: GNU and Intel. The COMPILE environment variable must be set to either ‘gnu’ or ‘intel’ before running build script. For example, to use Intel compilers when using the bash shell, run:

..  code-block:: console

    $ export COMPILER=intel
  
You can further customize compiler MPI wrappers by setting these four environment variables:

..  code-block:: console

    CMAKE_Platform : if not set the default is linux.${COMPILER}
    CMAKE_C_COMPILER : if not set the default is mpicc
    CMAKE_CXX_COMPILER : if not set the default is mpicxx
    CMAKE_Fortran_COMPILER : if not set the default is mpif90

In order to have one or more CCPP physics suites available at runtime, you need to select those  suites at build time by setting the ``CCPP_SUITES`` environment variable. Multiple suites can be set, as shown below in an example for the bash shell

..  code-block:: console

    $ export CCPP_SUITES=’FV3_GFS_2017_gfdlmp,FV3_GFS_v15’

If ``CCPP_SUITES`` is not set, the default is ``‘FV3_GFS_2017_gfdlmp’``.

After setting all the environment variables, you can build the model using

..  code-block:: console

    $ ./build.sh

Once ``build.sh`` is finished, you should see the executable, named ufs_weather_model, in the top-level directory.





